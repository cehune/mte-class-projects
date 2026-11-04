/* Plain C versions of the assembly kernels, used by tests and the benchmark. */
#include <math.h>
#include "nn.h"

void r_saxpy(long n, float a, const float *x, float *y) {
    for (long i = 0; i < n; i++) y[i] += a * x[i];
}
float r_sdot(long n, const float *x, const float *y) {
    float s = 0;
    for (long i = 0; i < n; i++) s += x[i] * y[i];
    return s;
}
float r_ssum(long n, const float *x) {
    float s = 0;
    for (long i = 0; i < n; i++) s += x[i];
    return s;
}
void r_scale_shift(long n, float a, float b, const float *src, float *dst) {
    for (long i = 0; i < n; i++) dst[i] = a * src[i] + b;
}
void r_relu(long n, const float *src, float *dst) {
    for (long i = 0; i < n; i++) dst[i] = src[i] > 0 ? src[i] : 0;
}
void r_relu_back(long n, const float *gy, const float *x, float *gx) {
    for (long i = 0; i < n; i++) gx[i] = x[i] > 0 ? gy[i] : 0;
}
void r_conv3_row(long n, float *dst, const float *r0, const float *r1, const float *r2,
                 const float *w) {
    const float *r[3] = {r0, r1, r2};
    for (long i = 0; i < n; i++)
        for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++) dst[i] += w[3 * k + l] * r[k][i - 1 + l];
}
void r_conv3_grad_row(long n, float *wg, const float *g, const float *r0, const float *r1,
                      const float *r2) {
    const float *r[3] = {r0, r1, r2};
    for (long i = 0; i < n; i++)
        for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++) wg[3 * k + l] += g[i] * r[k][i - 1 + l];
}
void r_maxpool2_row(long n_out, const float *r0, const float *r1, float *dst, int32_t *idx) {
    for (long i = 0; i < n_out; i++) {
        float v[4] = {r0[2 * i], r0[2 * i + 1], r1[2 * i], r1[2 * i + 1]};
        int best = 0;
        for (int k = 1; k < 4; k++)
            if (v[k] > v[best]) best = k;
        dst[i] = v[best];
        idx[i] = best;
    }
}
void r_adamw(long n, float *w, const float *g, float *m, float *v, float b1, float b2, float eps,
             float decay, float step, float root_inv) {
    for (long i = 0; i < n; i++) {
        w[i] *= decay;
        m[i] = b1 * m[i] + (1 - b1) * g[i];
        v[i] = b2 * v[i] + (1 - b2) * g[i] * g[i];
        w[i] -= step * m[i] / (sqrtf(v[i]) * root_inv + eps);
    }
}
