/* Each assembly kernel against its C reference, for every length from 0 to 67.
   Guard values after each output catch out-of-bounds writes. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nn.h"

#define MAX 67
#define PAD 8
static int fails = 0;
static float a[3][MAX * 2 + 2 * PAD], b[3][MAX * 2 + 2 * PAD], in1[MAX * 2 + 2 * PAD], in2[MAX * 2 + 2 * PAD],
    in3[MAX * 2 + 2 * PAD];

static float rnd(void) { return (rand() % 20001 - 10000) / 2500.0f; }

static void fill(float *x) {
    for (int i = 0; i < MAX * 2 + 2 * PAD; i++) x[i] = rnd();
}

static void reset(void) {
    for (int k = 0; k < 3; k++) {
        fill(a[k]);
        memcpy(b[k], a[k], sizeof a[k]);
    }
    fill(in1), fill(in2), fill(in3);
}

static void check(const char *name, long n, const float *x, const float *y, long count, float tol) {
    for (long i = 0; i < count; i++) {
        float limit = tol * (1 + fabsf(y[i]));
        if (!(fabsf(x[i] - y[i]) <= limit)) {
            printf("FAIL %s n=%ld i=%ld asm=%g ref=%g\n", name, n, i, x[i], y[i]);
            fails++;
            return;
        }
    }
}

#define SAME(name, n, tol) \
    for (int k = 0; k < 3; k++) check(name, n, a[k], b[k], sizeof a[k] / sizeof(float), tol)
#define EXACT 1e-5f
#define SUMS 5e-4f /* NEON adds in a different order, which changes rounding on long sums */

int main(void) {
    srand(577);
    for (int rep = 0; rep < 20; rep++)
        for (long n = 0; n <= MAX; n++) {
            float s = rnd(), t = rnd();
            const float *x = in1 + PAD, *y = in2 + PAD, *z = in3 + PAD;

            reset(); k_saxpy(n, s, x, a[0] + PAD); r_saxpy(n, s, x, b[0] + PAD); SAME("saxpy", n, EXACT);
            reset();
            float d1 = k_sdot(n, x, y), d2 = r_sdot(n, x, y);
            check("sdot", n, &d1, &d2, 1, 1e-4f);
            d1 = k_ssum(n, x), d2 = r_ssum(n, x);
            check("ssum", n, &d1, &d2, 1, 1e-4f);
            reset(); k_scale_shift(n, s, t, x, a[0] + PAD); r_scale_shift(n, s, t, x, b[0] + PAD); SAME("scale_shift", n, EXACT);
            reset(); k_relu(n, x, a[0] + PAD); r_relu(n, x, b[0] + PAD); SAME("relu", n, EXACT);
            reset(); k_relu_back(n, x, y, a[0] + PAD); r_relu_back(n, x, y, b[0] + PAD); SAME("relu_back", n, EXACT);
            reset();
            k_conv3_row(n, a[0] + PAD, x, y, z, in1 + 2 * PAD);
            r_conv3_row(n, b[0] + PAD, x, y, z, in1 + 2 * PAD);
            SAME("conv3_row", n, SUMS);
            reset();
            k_conv3_grad_row(n, a[0] + PAD, x, y, z, in1 + 3 * PAD);
            r_conv3_grad_row(n, b[0] + PAD, x, y, z, in1 + 3 * PAD);
            SAME("conv3_grad_row", n, SUMS);
            reset();
            int32_t i1[MAX + 2], i2[MAX + 2];
            memset(i1, 0x7f, sizeof i1), memset(i2, 0x7f, sizeof i2);
            k_maxpool2_row(n, x, y, a[0] + PAD, i1);
            r_maxpool2_row(n, x, y, b[0] + PAD, i2);
            SAME("maxpool2", n, EXACT);
            if (memcmp(i1, i2, sizeof i1)) { printf("FAIL maxpool idx n=%ld\n", n); fails++; }
            reset();
            for (int i = 0; i < MAX * 2 + 2 * PAD; i++) a[2][i] = b[2][i] = fabsf(a[2][i]);
            k_adamw(n, a[0] + PAD, x, a[1] + PAD, a[2] + PAD, 0.9f, 0.999f, 1e-8f, 0.999f, 0.01f, 1.5f);
            r_adamw(n, b[0] + PAD, x, b[1] + PAD, b[2] + PAD, 0.9f, 0.999f, 1e-8f, 0.999f, 0.01f, 1.5f);
            SAME("adamw", n, EXACT);
        }
    /* Ties: the first maximum wins, as in PyTorch */
    float r0[8] = {1, 1, 0, 2, 3, 0, 1, 1}, r1[8] = {1, 1, 2, 2, 3, 3, 0, 1}, out[4];
    int32_t idx[4], want[4] = {0, 1, 0, 0};
    k_maxpool2_row(4, r0, r1, out, idx);
    if (memcmp(idx, want, sizeof idx)) { printf("FAIL maxpool ties\n"); fails++; }
    printf("test_kernels: %s\n", fails ? "FAILED" : "ok");
    return fails != 0;
}
