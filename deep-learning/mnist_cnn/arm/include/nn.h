#ifndef NN_H
#define NN_H
#include <stddef.h>
#include <stdint.h>

/* Assembly kernels (src/kernels.S). r_* in kernels_ref.c share the signatures. */
#define KERNELS(P)                                                                   \
    void P##saxpy(long n, float a, const float *x, float *y);                        \
    float P##sdot(long n, const float *x, const float *y);                           \
    float P##ssum(long n, const float *x);                                           \
    void P##scale_shift(long n, float a, float b, const float *src, float *dst);     \
    void P##relu(long n, const float *src, float *dst);                              \
    void P##relu_back(long n, const float *gy, const float *x, float *gx);           \
    void P##conv3_row(long n, float *dst, const float *r0, const float *r1,          \
                      const float *r2, const float *w);                              \
    void P##conv3_grad_row(long n, float *wg, const float *g, const float *r0,       \
                           const float *r1, const float *r2);                        \
    void P##maxpool2_row(long n_out, const float *r0, const float *r1, float *dst,   \
                         int32_t *idx);                                              \
    void P##adamw(long n, float *w, const float *g, float *m, float *v, float b1,    \
                  float b2, float eps, float decay, float step, float root_inv);
KERNELS(k_)
KERNELS(r_)

/* Config, same keys as the notebook */
typedef struct {
    int n_blocks, base_filters, use_batchnorm, gap, batch_size, onecycle;
    float dropout, weight_decay, label_smoothing, lr;
    float rotation, translate, scale_jitter;
} Cfg;
void cfg_default(Cfg *c);
int cfg_load(Cfg *c, const char *path);
void cfg_save(const Cfg *c, const char *path);

/* RNG */
typedef struct { uint64_t s; } Rng;
uint32_t rng_u32(Rng *r);
float rng_uniform(Rng *r);               /* [0, 1) */
void rng_shuffle(Rng *r, int *a, int n);

/* Data */
typedef struct { int n; float *x; int *y; } Dataset; /* x is n * 784 */
Dataset idx_load(const char *images, const char *labels);
Dataset dataset_subset(const Dataset *d, const int *idx, int n);
void dataset_norm(const Dataset *train, float *mean, float *std);
void dataset_apply_norm(Dataset *d, float mean, float std);
int *split_load(const char *path, int *n);

/* Model */
enum { L_CONV, L_BN, L_RELU, L_POOL, L_GAP, L_FLAT, L_DROP, L_FC };
typedef struct { float *w, *g, *m, *v; long n; } Param;
typedef struct {
    int type, ci, hi, wi, co, ho, wo, p0, p1;
    float *x, *y, *gy, *pad, *gpad, *mask, *rmean, *rvar, *bmean, *bstd, *wg_local;
    int32_t *idx;
    float p;
} Layer;
typedef struct {
    Cfg cfg;
    int n_layers, n_params, max_batch, N, train, step;
    Layer l[64];
    Param p[32];
    float *input, *ginput;
    Rng rng;
    float mean, std;
    float loss;
    int correct;
} Net;

Net *net_build(const Cfg *c, int max_batch, uint64_t seed);
float *net_forward(Net *net, const float *x, int N, int train); /* returns logits */
void net_loss_backward(Net *net, const int *y);
long net_param_count(const Net *net);
void net_zero_grad(Net *net);
int net_save(const Net *net, const char *path);
int net_load(Net *net, const char *path);
double net_accuracy(Net *net, const Dataset *d, int *confusion);

/* Optimizer */
void onecycle_at(float max_lr, long s, long total, float *lr, float *b1);
void opt_step(Net *net, float lr, float b1);

/* Augmentation */
void augment_theta(const Cfg *c, Rng *r, float theta[6]);
void augment_apply(const float *src, float *dst, const float theta[6], float background);
void augment_batch(const Cfg *c, Rng *r, const float *src, float *dst, int N, float background);

#endif
