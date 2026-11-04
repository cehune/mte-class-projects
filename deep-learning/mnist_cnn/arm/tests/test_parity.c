/* Native training must reproduce PyTorch: logits, loss, gradients, BN stats and parameters
   after several AdamW + OneCycle steps. See make_fixtures.py for the file format. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nn.h"

#define STEPS 5
#define TOTAL_STEPS 20
static int fails = 0;

static void rd(FILE *f, void *p, size_t bytes) {
    if (fread(p, 1, bytes, f) != bytes) { printf("short fixture\n"); exit(1); }
}

static void cmp(const char *what, const float *got, const float *want, long n, float tol) {
    double worst = 0;
    long at = 0;
    for (long i = 0; i < n; i++) {
        double e = fabs((double)got[i] - want[i]);
        if (e > worst) worst = e, at = i;
    }
    printf("  %-22s max error %.2e\n", what, worst);
    if (worst > tol) {
        printf("FAIL %s at %ld: got %g want %g\n", what, at, got[at], want[at]);
        fails++;
    }
}

static void parity(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(1); }
    printf("%s\n", path);
    float h[13];
    int N;
    rd(f, h, sizeof h);
    rd(f, &N, sizeof N);
    Cfg c;
    cfg_default(&c);
    c.n_blocks = (int)h[0], c.base_filters = (int)h[1], c.use_batchnorm = (int)h[2], c.gap = (int)h[3];
    c.batch_size = N, c.weight_decay = h[7], c.label_smoothing = h[8], c.lr = h[9];
    Net *net = net_build(&c, N, 1);

    float *x = malloc((long)N * 784 * sizeof(float));
    int *y = malloc(N * sizeof(int));
    rd(f, x, (long)N * 784 * sizeof(float));
    rd(f, y, N * sizeof(int));
    /* The parameter order differs for BN models only in being weight, then BN, as in PyTorch */
    for (int i = 0; i < net->n_params; i++) rd(f, net->p[i].w, net->p[i].n * sizeof(float));

    float logits[160], loss;
    rd(f, logits, (long)N * 10 * sizeof(float));
    rd(f, &loss, sizeof loss);
    float **grads = malloc(sizeof(float *) * net->n_params);
    for (int i = 0; i < net->n_params; i++) {
        grads[i] = malloc(net->p[i].n * sizeof(float));
        rd(f, grads[i], net->p[i].n * sizeof(float));
    }

    for (int s = 0; s < STEPS; s++) {
        const float *z = net_forward(net, x, N, 1);
        if (s == 0) cmp("logits", z, logits, (long)N * 10, 2e-3f);
        net_zero_grad(net);
        net_loss_backward(net, y);
        if (s == 0) {
            cmp("loss", &net->loss, &loss, 1, 1e-4f);
            for (int i = 0; i < net->n_params; i++) {
                char name[32];
                snprintf(name, sizeof name, "grad param %d", i);
                cmp(name, net->p[i].g, grads[i], net->p[i].n, 2e-3f);
            }
        }
        float lr, b1;
        onecycle_at(c.lr, s, TOTAL_STEPS, &lr, &b1);
        opt_step(net, lr, b1);
    }
    for (int i = 0; i < net->n_params; i++) {
        float *want = malloc(net->p[i].n * sizeof(float));
        rd(f, want, net->p[i].n * sizeof(float));
        char name[32];
        snprintf(name, sizeof name, "param %d after %d steps", i, STEPS);
        cmp(name, net->p[i].w, want, net->p[i].n, 5e-4f);
        free(want);
    }
    for (int i = 0; i < net->n_layers; i++) {
        Layer *l = &net->l[i];
        if (l->type != L_BN) continue;
        float *want = malloc(l->co * sizeof(float));
        rd(f, want, l->co * sizeof(float));
        cmp("BN running mean", l->rmean, want, l->co, 2e-2f);
        rd(f, want, l->co * sizeof(float));
        cmp("BN running var", l->rvar, want, l->co, 2e-2f);
        free(want);
    }
    fclose(f);
}

int main(void) {
    parity("tests/fixtures/bn_flatten.bin");
    parity("tests/fixtures/plain_gap.bin");

    FILE *f = fopen("tests/fixtures/augment.bin", "rb");
    if (!f) { perror("augment fixture"); return 1; }
    float theta[6], image[784], want[784], got[784];
    rd(f, theta, sizeof theta);
    rd(f, image, sizeof image);
    rd(f, want, sizeof want);
    fclose(f);
    augment_apply(image, got, theta, -0.4241f);
    printf("augment\n");
    cmp("grid_sample", got, want, 784, 2e-3f);

    printf("test_parity: %s\n", fails ? "FAILED" : "ok");
    return fails != 0;
}
