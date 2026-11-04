/* Per-image affine augmentation. Same transform as augment() in the notebook. */
#include <dispatch/dispatch.h>
#include <math.h>
#include <stdlib.h>
#include "nn.h"

void augment_theta(const Cfg *c, Rng *r, float t[6]) {
    float u[4];
    for (int i = 0; i < 4; i++) u[i] = 2 * rng_uniform(r) - 1;
    float angle = u[0] * c->rotation * (float)M_PI / 180;
    float scale = 1 + u[1] * c->scale_jitter;
    float cs = cosf(angle), sn = sinf(angle);
    t[0] = cs / scale, t[1] = -sn / scale, t[2] = u[2] * c->translate * 2;
    t[3] = sn / scale, t[4] = cs / scale, t[5] = u[3] * c->translate * 2;
}

/* Equivalent to affine_grid + grid_sample with bilinear sampling, align_corners=False.
   Pixels outside the image take the digit background value. */
void augment_apply(const float *src, float *dst, const float t[6], float bg) {
    for (int y = 0; y < 28; y++) {
        float gy = (2 * y + 1) / 28.0f - 1;
        for (int x = 0; x < 28; x++) {
            float gx = (2 * x + 1) / 28.0f - 1;
            float sx = ((t[0] * gx + t[1] * gy + t[2] + 1) * 28 - 1) / 2;
            float sy = ((t[3] * gx + t[4] * gy + t[5] + 1) * 28 - 1) / 2;
            int x0 = (int)floorf(sx), y0 = (int)floorf(sy);
            float fx = sx - x0, fy = sy - y0, out = 0;
            for (int j = 0; j < 2; j++)
                for (int i = 0; i < 2; i++) {
                    int px = x0 + i, py = y0 + j;
                    if (px < 0 || px > 27 || py < 0 || py > 27) continue;
                    float weight = (i ? fx : 1 - fx) * (j ? fy : 1 - fy);
                    out += weight * (src[py * 28 + px] - bg);
                }
            dst[y * 28 + x] = out + bg;
        }
    }
}

void augment_batch(const Cfg *c, Rng *r, const float *src, float *dst, int N, float bg) {
    float *theta = malloc(sizeof(float) * 6 * N);
    for (int n = 0; n < N; n++) augment_theta(c, r, theta + n * 6);
    dispatch_apply(N, dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0), ^(size_t n) {
        augment_apply(src + n * 784, dst + n * 784, theta + n * 6, bg);
    });
    free(theta);
}
