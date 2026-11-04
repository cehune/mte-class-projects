/* AdamW and OneCycle, following the formulas in torch.optim */
#include <math.h>
#include "nn.h"

static float anneal(float start, float end, float pct) {
    return end + (start - end) / 2 * (1 + cosf((float)M_PI * pct));
}

void onecycle_at(float max_lr, long s, long total, float *lr, float *b1) {
    float initial = max_lr / 25, floor_lr = initial / 1e4f;
    float peak = 0.3f * total - 1;
    if (s <= peak) {
        float pct = s / peak;
        *lr = anneal(initial, max_lr, pct);
        *b1 = anneal(0.95f, 0.85f, pct);
    } else {
        float pct = (s - peak) / (total - 1 - peak);
        *lr = anneal(max_lr, floor_lr, pct);
        *b1 = anneal(0.85f, 0.95f, pct);
    }
}

void opt_step(Net *net, float lr, float b1) {
    const float b2 = 0.999f;
    net->step++;
    float step = lr / (1 - powf(b1, net->step));
    float root_inv = 1 / sqrtf(1 - powf(b2, net->step));
    float decay = 1 - lr * net->cfg.weight_decay;
    for (int i = 0; i < net->n_params; i++) {
        Param *p = &net->p[i];
        k_adamw(p->n, p->w, p->g, p->m, p->v, b1, b2, 1e-8f, decay, step, root_inv);
    }
}
