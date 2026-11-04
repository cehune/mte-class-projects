/* Layer graph with hand-derived forward and backward passes over the assembly kernels. */
#include <dispatch/dispatch.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nn.h"

#define BN_EPS 1e-5f
#define BN_MOMENTUM 0.1f

static dispatch_queue_t Q(void) { return dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0); }

/* Row stride of padded buffers. Rows are rounded up to a multiple of 4, so the kernels never
   need their slow scalar tail, and the extra columns stay zero. */
static int row4(int w) { return (w + 3) & ~3; }

static float *fzalloc(long n) {
    float *p = calloc(n, sizeof(float));
    if (!p) { fprintf(stderr, "out of memory\n"); exit(1); }
    return p;
}

static int add_param(Net *net, long n) {
    Param *p = &net->p[net->n_params];
    p->n = n;
    p->w = fzalloc(n);
    p->g = fzalloc(n);
    p->m = fzalloc(n);
    p->v = fzalloc(n);
    return net->n_params++;
}

/* Matches PyTorch's default kaiming_uniform(a=sqrt(5)): U(-1/sqrt(fan_in), 1/sqrt(fan_in)) */
static void init_uniform(Rng *r, float *w, long n, float bound) {
    for (long i = 0; i < n; i++) w[i] = (2 * rng_uniform(r) - 1) * bound;
}

static Layer *add_layer(Net *net, int type, int co, int ho, int wo) {
    Layer *l = &net->l[net->n_layers];
    memset(l, 0, sizeof(*l));
    l->type = type;
    if (net->n_layers) {
        Layer *prev = l - 1;
        l->ci = prev->co, l->hi = prev->ho, l->wi = prev->wo;
    } else {
        l->ci = 1, l->hi = 28, l->wi = 28;
    }
    l->co = co, l->ho = ho, l->wo = wo;
    l->p0 = l->p1 = -1;
    long size = (long)net->max_batch * co * ho * wo;
    l->y = fzalloc(size);
    l->gy = fzalloc(size);
    net->n_layers++;
    return l;
}

Net *net_build(const Cfg *c, int max_batch, uint64_t seed) {
    Net *net = calloc(1, sizeof(Net));
    net->cfg = *c;
    net->max_batch = max_batch;
    net->rng.s = seed;
    net->mean = 0, net->std = 1;
    net->input = fzalloc((long)max_batch * 784);
    net->ginput = fzalloc((long)max_batch * 784);
    Rng *r = &net->rng;

    int ch = 1, size = 28;
    for (int b = 0; b < c->n_blocks; b++) {
        int out = c->base_filters << b;
        for (int i = 0; i < 2; i++) {
            Layer *l = add_layer(net, L_CONV, out, size, size);
            l->pad = fzalloc((long)max_batch * ch * (size + 2) * (row4(size) + 2));
            l->gpad = fzalloc((long)max_batch * out * (size + 2) * (row4(size) + 2));
            l->p0 = add_param(net, (long)out * ch * 9);
            init_uniform(r, net->p[l->p0].w, net->p[l->p0].n, 1 / sqrtf(ch * 9));
            if (!c->use_batchnorm) {
                l->p1 = add_param(net, out);
                init_uniform(r, net->p[l->p1].w, out, 1 / sqrtf(ch * 9));
            } else {
                Layer *n = add_layer(net, L_BN, out, size, size);
                n->p0 = add_param(net, out);
                n->p1 = add_param(net, out);
                for (int k = 0; k < out; k++) net->p[n->p0].w[k] = 1;
                n->rmean = fzalloc(out), n->rvar = fzalloc(out);
                n->bmean = fzalloc(out), n->bstd = fzalloc(out);
                for (int k = 0; k < out; k++) n->rvar[k] = 1;
            }
            add_layer(net, L_RELU, out, size, size);
            ch = out;
        }
        Layer *p = add_layer(net, L_POOL, ch, size / 2, size / 2);
        p->idx = malloc(sizeof(int32_t) * (long)max_batch * ch * (size / 2) * (size / 2));
        size /= 2;
    }
    int features;
    if (c->gap) {
        add_layer(net, L_GAP, ch, 1, 1);
        features = ch;
    } else {
        features = ch * size * size;
        add_layer(net, L_FLAT, features, 1, 1);
    }
    Layer *d = add_layer(net, L_DROP, features, 1, 1);
    d->p = c->dropout;
    d->mask = fzalloc((long)max_batch * features);
    Layer *fc = add_layer(net, L_FC, 10, 1, 1);
    fc->p0 = add_param(net, 10L * features);
    fc->p1 = add_param(net, 10);
    init_uniform(r, net->p[fc->p0].w, 10L * features, 1 / sqrtf(features));
    init_uniform(r, net->p[fc->p1].w, 10, 1 / sqrtf(features));

    for (int i = 0; i < net->n_layers; i++) {
        Layer *l = &net->l[i];
        l->x = i ? net->l[i - 1].y : net->input;
    }
    return net;
}

long net_param_count(const Net *net) {
    long n = 0;
    for (int i = 0; i < net->n_params; i++) n += net->p[i].n;
    return n;
}

void net_zero_grad(Net *net) {
    for (int i = 0; i < net->n_params; i++) memset(net->p[i].g, 0, net->p[i].n * sizeof(float));
}

/* Copies each H*W image into a buffer with a zero border */
static void pad_images(const float *src, float *dst, long count, int h, int w) {
    long ph = h + 2, pw = row4(w) + 2;
    memset(dst, 0, count * ph * pw * sizeof(float));
    for (long i = 0; i < count; i++)
        for (int y = 0; y < h; y++)
            memcpy(dst + i * ph * pw + (y + 1) * pw + 1, src + i * h * w + y * w, w * sizeof(float));
}

/* Work is split by sample, so each task stays inside one image's features and in cache */
static void conv_forward(Net *net, Layer *l) {
    int N = net->N, C = l->ci, O = l->co, H = l->hi, W = l->wi, W4 = row4(W);
    long HW = H * W, pw = W4 + 2, ph = H + 2;
    const float *w = net->p[l->p0].w;
    const float *bias = l->p1 >= 0 ? net->p[l->p1].w : NULL;
    pad_images(l->x, l->pad, (long)N * C, H, W);
    dispatch_apply(N, Q(), ^(size_t n) {
        float row[H * W4];
        for (int o = 0; o < O; o++) {
            memset(row, 0, sizeof row);
            for (int c = 0; c < C; c++) {
                const float *pc = l->pad + ((long)n * C + c) * ph * pw;
                const float *wc = w + (o * (long)C + c) * 9;
                for (int y = 0; y < H; y++)
                    k_conv3_row(W4, row + y * W4, pc + y * pw + 1, pc + (y + 1) * pw + 1,
                                pc + (y + 2) * pw + 1, wc);
            }
            float *dst = l->y + (long)n * O * HW + o * HW;
            for (int y = 0; y < H; y++) memcpy(dst + y * W, row + y * W4, W * sizeof(float));
            if (bias) k_scale_shift(HW, 1, bias[o], dst, dst);
        }
    });
}

#define GRAD_WORKERS 18

static void conv_backward(Net *net, Layer *l, float *gx) {
    int N = net->N, C = l->ci, O = l->co, H = l->hi, W = l->wi, W4 = row4(W);
    long HW = H * W, pw = W4 + 2, ph = H + 2, nw = net->p[l->p0].n;
    const float *w = net->p[l->p0].w;
    float *wg = net->p[l->p0].g;
    if (!l->wg_local) l->wg_local = fzalloc(nw * GRAD_WORKERS);
    float *locals = l->wg_local;
    float *flip = malloc(nw * sizeof(float)); /* transposed conv uses each kernel rotated 180 degrees */
    for (long i = 0; i < nw; i += 9)
        for (int k = 0; k < 9; k++) flip[i + k] = w[i + 8 - k];
    pad_images(l->gy, l->gpad, (long)N * O, H, W);
    /* Weight gradients are shared across samples, so each worker fills its own copy */
    dispatch_apply(GRAD_WORKERS, Q(), ^(size_t t) {
        float *local = locals + t * nw;
        memset(local, 0, nw * sizeof(float));
        float row[H * W4];
        for (int n = (int)(t * N / GRAD_WORKERS); n < (int)((t + 1) * N / GRAD_WORKERS); n++)
            for (int c = 0; c < C; c++) {
                memset(row, 0, sizeof row);
                const float *px = l->pad + ((long)n * C + c) * ph * pw;
                for (int o = 0; o < O; o++) {
                    const float *pg = l->gpad + ((long)n * O + o) * ph * pw;
                    long at = (o * (long)C + c) * 9;
                    for (int y = 0; y < H; y++) {
                        k_conv3_row(W4, row + y * W4, pg + y * pw + 1, pg + (y + 1) * pw + 1,
                                    pg + (y + 2) * pw + 1, flip + at);
                        /* the padded row is zero past W, so extra columns add nothing */
                        k_conv3_grad_row(W4, local + at, pg + (y + 1) * pw + 1, px + y * pw + 1,
                                         px + (y + 1) * pw + 1, px + (y + 2) * pw + 1);
                    }
                }
                float *dst = gx + (long)n * C * HW + c * HW;
                for (int y = 0; y < H; y++) memcpy(dst + y * W, row + y * W4, W * sizeof(float));
            }
    });
    free(flip);
    for (int t = 0; t < GRAD_WORKERS; t++) k_saxpy(nw, 1, locals + t * nw, wg);
    if (l->p1 >= 0) {
        float *bg = net->p[l->p1].g;
        for (int o = 0; o < O; o++)
            for (int n = 0; n < N; n++) bg[o] += k_ssum(HW, l->gy + (long)n * O * HW + o * HW);
    }
}

static void bn_forward(Net *net, Layer *l) {
    int N = net->N, C = l->ci, train = net->train;
    long HW = (long)l->hi * l->wi, M = N * HW;
    const float *gamma = net->p[l->p0].w, *beta = net->p[l->p1].w;
    dispatch_apply(C, Q(), ^(size_t c) {
        float mean = l->rmean[c], var = l->rvar[c];
        if (train) {
            double s1 = 0, s2 = 0;
            for (int n = 0; n < N; n++) {
                const float *xc = l->x + (long)n * C * HW + c * HW;
                s1 += k_ssum(HW, xc);
                s2 += k_sdot(HW, xc, xc);
            }
            double dm = s1 / M, dv = s2 / M - dm * dm;
            mean = (float)dm;
            var = dv > 0 ? (float)dv : 0;
            l->rmean[c] = (1 - BN_MOMENTUM) * l->rmean[c] + BN_MOMENTUM * mean;
            l->rvar[c] = (1 - BN_MOMENTUM) * l->rvar[c] + BN_MOMENTUM * var;
        }
        float s = sqrtf(var + BN_EPS), a = gamma[c] / s, b = beta[c] - mean * a;
        l->bmean[c] = mean, l->bstd[c] = s;
        for (int n = 0; n < N; n++) {
            long off = (long)n * C * HW + c * HW;
            k_scale_shift(HW, a, b, l->x + off, l->y + off);
        }
    });
}

/* Backward through the batch statistics too, since mean and std depend on every input:
   gx_i = gamma/s (g_i - gbar) - gamma z_i sum(g z) / (M s), with z = (x - mean) / s */
static void bn_backward(Net *net, Layer *l, float *gx) {
    int N = net->N, C = l->ci;
    long HW = (long)l->hi * l->wi, M = N * HW;
    const float *gamma = net->p[l->p0].w;
    float *gg = net->p[l->p0].g, *bg = net->p[l->p1].g;
    dispatch_apply(C, Q(), ^(size_t c) {
        double gsum = 0, gx_dot = 0;
        for (int n = 0; n < N; n++) {
            long off = (long)n * C * HW + c * HW;
            gsum += k_ssum(HW, l->gy + off);
            gx_dot += k_sdot(HW, l->gy + off, l->x + off);
        }
        float s = l->bstd[c], m = l->bmean[c], g0 = gamma[c];
        double gz = (gx_dot - m * gsum) / s;
        float a = g0 / s, gbar = (float)(gsum / M);
        float k = (float)(-g0 * gz / (M * s * s));
        float shift = -a * gbar - k * m;
        for (int n = 0; n < N; n++) {
            long off = (long)n * C * HW + c * HW;
            k_scale_shift(HW, a, shift, l->gy + off, gx + off);
            k_saxpy(HW, k, l->x + off, gx + off);
        }
        gg[c] += (float)gz;
        bg[c] += (float)gsum;
    });
}

static void pool_forward(Net *net, Layer *l) {
    int N = net->N, C = l->ci, H = l->hi, W = l->wi, ho = l->ho, wo = l->wo;
    dispatch_apply((long)N * C, Q(), ^(size_t i) {
        const float *src = l->x + i * H * W;
        for (int y = 0; y < ho; y++)
            k_maxpool2_row(wo, src + 2 * y * W, src + (2 * y + 1) * W,
                           l->y + i * ho * wo + y * wo, l->idx + i * ho * wo + y * wo);
    });
}

static void pool_backward(Net *net, Layer *l, float *gx) {
    int N = net->N, C = l->ci, H = l->hi, W = l->wi, ho = l->ho, wo = l->wo;
    dispatch_apply((long)N * C, Q(), ^(size_t i) {
        float *dst = gx + i * H * W;
        memset(dst, 0, (long)H * W * sizeof(float));
        const float *g = l->gy + i * ho * wo;
        const int32_t *idx = l->idx + i * ho * wo;
        for (int y = 0; y < ho; y++)
            for (int v = 0; v < wo; v++) {
                int which = idx[y * wo + v];
                dst[(2 * y + (which >> 1)) * W + 2 * v + (which & 1)] += g[y * wo + v];
            }
    });
}

static void fc_forward(Net *net, Layer *l) {
    int N = net->N, D = l->ci;
    const float *w = net->p[l->p0].w, *b = net->p[l->p1].w;
    for (int n = 0; n < N; n++)
        for (int k = 0; k < 10; k++)
            l->y[n * 10 + k] = k_sdot(D, w + (long)k * D, l->x + (long)n * D) + b[k];
}

static void fc_backward(Net *net, Layer *l, float *gx) {
    int N = net->N, D = l->ci;
    const float *w = net->p[l->p0].w;
    float *wg = net->p[l->p0].g, *bg = net->p[l->p1].g;
    memset(gx, 0, (long)N * D * sizeof(float));
    for (int n = 0; n < N; n++)
        for (int k = 0; k < 10; k++) {
            float g = l->gy[n * 10 + k];
            k_saxpy(D, g, l->x + (long)n * D, wg + (long)k * D);
            k_saxpy(D, g, w + (long)k * D, gx + (long)n * D);
            bg[k] += g;
        }
}

float *net_forward(Net *net, const float *x, int N, int train) {
    net->N = N;
    net->train = train;
    memcpy(net->input, x, (long)N * 784 * sizeof(float));
    for (int i = 0; i < net->n_layers; i++) {
        Layer *l = &net->l[i];
        long total = (long)N * l->co * l->ho * l->wo;
        switch (l->type) {
        case L_CONV: conv_forward(net, l); break;
        case L_BN: bn_forward(net, l); break;
        case L_RELU: k_relu(total, l->x, l->y); break;
        case L_POOL: pool_forward(net, l); break;
        case L_GAP: {
            long HW = (long)l->hi * l->wi;
            for (long k = 0; k < total; k++) l->y[k] = k_ssum(HW, l->x + k * HW) / HW;
            break;
        }
        case L_FLAT: memcpy(l->y, l->x, total * sizeof(float)); break;
        case L_DROP:
            if (train && l->p > 0) {
                float keep = 1 / (1 - l->p);
                for (long k = 0; k < total; k++) {
                    l->mask[k] = rng_uniform(&net->rng) < l->p ? 0 : keep;
                    l->y[k] = l->x[k] * l->mask[k];
                }
            } else {
                memcpy(l->y, l->x, total * sizeof(float));
            }
            break;
        case L_FC: fc_forward(net, l); break;
        }
    }
    return net->l[net->n_layers - 1].y;
}

void net_loss_backward(Net *net, const int *y) {
    int N = net->N;
    float eps = net->cfg.label_smoothing;
    Layer *last = &net->l[net->n_layers - 1];
    double loss = 0;
    int correct = 0;
    for (int n = 0; n < N; n++) {
        const float *z = last->y + n * 10;
        float *g = last->gy + n * 10;
        float top = z[0];
        int arg = 0;
        for (int k = 1; k < 10; k++)
            if (z[k] > top) top = z[k], arg = k;
        correct += arg == y[n];
        double sum = 0, p[10];
        for (int k = 0; k < 10; k++) sum += exp((double)z[k] - top);
        for (int k = 0; k < 10; k++) {
            p[k] = exp((double)z[k] - top) / sum;
            double target = eps / 10 + (k == y[n] ? 1 - eps : 0);
            loss -= target * log(p[k]);
            g[k] = (float)((p[k] - target) / N);
        }
    }
    net->loss = (float)(loss / N);
    net->correct = correct;

    for (int i = net->n_layers - 1; i >= 0; i--) {
        Layer *l = &net->l[i];
        float *gx = i ? net->l[i - 1].gy : net->ginput;
        long total = (long)N * l->ci * l->hi * l->wi;
        switch (l->type) {
        case L_CONV: conv_backward(net, l, gx); break;
        case L_BN: bn_backward(net, l, gx); break;
        case L_RELU: k_relu_back(total, l->gy, l->x, gx); break;
        case L_POOL: pool_backward(net, l, gx); break;
        case L_GAP: {
            long HW = (long)l->hi * l->wi;
            for (long k = 0; k < (long)N * l->ci; k++)
                for (long j = 0; j < HW; j++) gx[k * HW + j] = l->gy[k] / HW;
            break;
        }
        case L_FLAT: memcpy(gx, l->gy, total * sizeof(float)); break;
        case L_DROP:
            if (net->train && l->p > 0)
                for (long k = 0; k < total; k++) gx[k] = l->gy[k] * l->mask[k];
            else
                memcpy(gx, l->gy, total * sizeof(float));
            break;
        case L_FC: fc_backward(net, l, gx); break;
        }
    }
}

double net_accuracy(Net *net, const Dataset *d, int *confusion) {
    int correct = 0;
    for (int i = 0; i < d->n; i += net->max_batch) {
        int N = d->n - i < net->max_batch ? d->n - i : net->max_batch;
        const float *z = net_forward(net, d->x + (long)i * 784, N, 0);
        for (int n = 0; n < N; n++) {
            int arg = 0;
            for (int k = 1; k < 10; k++)
                if (z[n * 10 + k] > z[n * 10 + arg]) arg = k;
            correct += arg == d->y[i + n];
            if (confusion) confusion[d->y[i + n] * 10 + arg]++;
        }
    }
    return (double)correct / d->n;
}

int net_save(const Net *net, const char *path) {
    FILE *f = fopen(path, "wb");
    if (!f) return -1;
    fwrite(&net->mean, sizeof(float), 1, f);
    fwrite(&net->std, sizeof(float), 1, f);
    for (int i = 0; i < net->n_params; i++) fwrite(net->p[i].w, sizeof(float), net->p[i].n, f);
    for (int i = 0; i < net->n_layers; i++)
        if (net->l[i].type == L_BN) {
            fwrite(net->l[i].rmean, sizeof(float), net->l[i].co, f);
            fwrite(net->l[i].rvar, sizeof(float), net->l[i].co, f);
        }
    return fclose(f);
}

int net_load(Net *net, const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) return -1;
    size_t ok = fread(&net->mean, sizeof(float), 1, f);
    ok += fread(&net->std, sizeof(float), 1, f);
    size_t want = 2;
    for (int i = 0; i < net->n_params; i++) {
        ok += fread(net->p[i].w, sizeof(float), net->p[i].n, f);
        want += net->p[i].n;
    }
    for (int i = 0; i < net->n_layers; i++)
        if (net->l[i].type == L_BN) {
            ok += fread(net->l[i].rmean, sizeof(float), net->l[i].co, f);
            ok += fread(net->l[i].rvar, sizeof(float), net->l[i].co, f);
            want += 2 * net->l[i].co;
        }
    fclose(f);
    return ok == want ? 0 : -1;
}
