/* Finite-difference gradient checks, dropout behaviour, checkpoint round trip, schedule shape */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nn.h"

static int fails = 0;
#define EXPECT(cond, ...) do { if (!(cond)) { printf("FAIL " __VA_ARGS__); printf("\n"); fails++; } } while (0)

static float loss_at(Net *net, const float *x, const int *y, int N, uint64_t seed) {
    net->rng.s = seed; /* the same dropout mask for each evaluation */
    net_forward(net, x, N, 1);
    net_zero_grad(net);
    net_loss_backward(net, y);
    return net->loss;
}

static void fd_check(const char *name, Cfg c) {
    int N = 6;
    Net *net = net_build(&c, N, 11);
    float *x = malloc(sizeof(float) * 784 * N);
    int y[6];
    Rng r = {5};
    for (int i = 0; i < 784 * N; i++) x[i] = 2 * rng_uniform(&r) - 1;
    for (int i = 0; i < N; i++) y[i] = i;

    int checked = 0, bad = 0;
    for (int p = 0; p < net->n_params; p++) {
        Param *q = &net->p[p];
        for (int t = 0; t < 8; t++) {
            long i = rng_u32(&r) % q->n;
            loss_at(net, x, y, N, 99);
            float analytic = q->g[i], saved = q->w[i], h = 2e-3f;
            if (fabsf(analytic) < 1e-3f) continue;
            q->w[i] = saved + h;
            float up = loss_at(net, x, y, N, 99);
            q->w[i] = saved - h;
            float down = loss_at(net, x, y, N, 99);
            q->w[i] = saved;
            float numeric = (up - down) / (2 * h);
            checked++;
            if (fabsf(numeric - analytic) > 0.15f * fabsf(analytic) + 3e-3f) {
                bad++;
                printf("  %s param %d[%ld]: analytic %.5f numeric %.5f\n", name, p, i, analytic, numeric);
            }
        }
    }
    /* ReLU and max pool kinks make a few mismatches possible, so allow a small share */
    EXPECT(checked >= 10 && bad * 10 <= checked, "%s gradient check: %d bad of %d", name, bad, checked);
    printf("  %s: %d of %d coordinates agree\n", name, checked - bad, checked);
}

int main(void) {
    Cfg c;
    cfg_default(&c);
    c.n_blocks = 2, c.base_filters = 4, c.batch_size = 6;
    fd_check("bn+flatten", c);
    c.use_batchnorm = 0, c.gap = 1;
    fd_check("plain+gap", c);
    c.use_batchnorm = 1, c.gap = 0, c.dropout = 0.3f, c.label_smoothing = 0.1f;
    fd_check("bn+dropout+smoothing", c);

    /* Dropout: rate, scale, and off during evaluation */
    Net *net = net_build(&c, 64, 3);
    float *x = calloc(784 * 64, sizeof(float));
    net_forward(net, x, 64, 1);
    Layer *d = NULL;
    for (int i = 0; i < net->n_layers; i++)
        if (net->l[i].type == L_DROP) d = &net->l[i];
    long zero = 0, total = 64L * d->co;
    for (long i = 0; i < total; i++) {
        zero += d->mask[i] == 0;
        EXPECT(d->mask[i] == 0 || fabsf(d->mask[i] - 1 / 0.7f) < 1e-5f, "dropout scale");
    }
    EXPECT(fabs((double)zero / total - 0.3) < 0.02, "dropout rate %f", (double)zero / total);
    net_forward(net, x, 64, 0);
    EXPECT(memcmp(d->x, d->y, total * sizeof(float)) == 0, "dropout off in eval");

    /* Checkpoint round trip */
    Rng r = {8};
    for (int i = 0; i < 784 * 8; i++) x[i] = rng_uniform(&r);
    net_forward(net, x, 8, 1);
    float *before = malloc(sizeof(float) * 80);
    memcpy(before, net_forward(net, x, 8, 0), sizeof(float) * 80);
    net_save(net, "build/ckpt_test.bin");
    Net *other = net_build(&c, 64, 4);
    EXPECT(net_load(other, "build/ckpt_test.bin") == 0, "checkpoint load");
    EXPECT(memcmp(before, net_forward(other, x, 8, 0), sizeof(float) * 80) == 0, "checkpoint outputs differ");

    /* OneCycle endpoints */
    float lr, b1;
    onecycle_at(1e-3f, 0, 1000, &lr, &b1);
    EXPECT(fabsf(lr - 4e-5f) < 1e-7f && fabsf(b1 - 0.95f) < 1e-5f, "onecycle start");
    onecycle_at(1e-3f, 299, 1000, &lr, &b1);
    EXPECT(fabsf(lr - 1e-3f) < 1e-7f && fabsf(b1 - 0.85f) < 1e-5f, "onecycle peak");
    onecycle_at(1e-3f, 999, 1000, &lr, &b1);
    EXPECT(lr < 1e-8f && fabsf(b1 - 0.95f) < 1e-5f, "onecycle end");

    printf("test_model: %s\n", fails ? "FAILED" : "ok");
    return fails != 0;
}
