/* MNIST IDX files, the validation split, configs and the RNG */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nn.h"

uint32_t rng_u32(Rng *r) { /* splitmix64 */
    uint64_t z = (r->s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return (uint32_t)((z ^ (z >> 31)) >> 32);
}

float rng_uniform(Rng *r) { return (rng_u32(r) >> 8) * (1.0f / 16777216.0f); }

void rng_shuffle(Rng *r, int *a, int n) {
    for (int i = n - 1; i > 0; i--) {
        int j = rng_u32(r) % (i + 1), t = a[i];
        a[i] = a[j], a[j] = t;
    }
}

static unsigned char *read_file(const char *path, long *size) {
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(1); }
    fseek(f, 0, SEEK_END);
    *size = ftell(f);
    rewind(f);
    unsigned char *b = malloc(*size);
    if (fread(b, 1, *size, f) != (size_t)*size) { perror(path); exit(1); }
    fclose(f);
    return b;
}

static int be32(const unsigned char *p) { return p[0] << 24 | p[1] << 16 | p[2] << 8 | p[3]; }

Dataset idx_load(const char *images, const char *labels) {
    long si, sl;
    unsigned char *im = read_file(images, &si), *lb = read_file(labels, &sl);
    Dataset d = {be32(im + 4), NULL, NULL};
    if (be32(im) != 2051 || be32(lb) != 2049 || be32(lb + 4) != d.n) {
        fprintf(stderr, "bad IDX files\n");
        exit(1);
    }
    d.x = malloc((long)d.n * 784 * sizeof(float));
    d.y = malloc(d.n * sizeof(int));
    for (long i = 0; i < (long)d.n * 784; i++) d.x[i] = im[16 + i] / 255.0f;
    for (int i = 0; i < d.n; i++) d.y[i] = lb[8 + i];
    free(im), free(lb);
    return d;
}

Dataset dataset_subset(const Dataset *d, const int *idx, int n) {
    Dataset s = {n, malloc((long)n * 784 * sizeof(float)), malloc(n * sizeof(int))};
    for (int i = 0; i < n; i++) {
        memcpy(s.x + (long)i * 784, d->x + (long)idx[i] * 784, 784 * sizeof(float));
        s.y[i] = d->y[idx[i]];
    }
    return s;
}

void dataset_norm(const Dataset *d, float *mean, float *std) {
    double s1 = 0, s2 = 0;
    long M = (long)d->n * 784;
    for (long i = 0; i < M; i++) s1 += d->x[i], s2 += (double)d->x[i] * d->x[i];
    *mean = (float)(s1 / M);
    *std = (float)sqrt((s2 - s1 * s1 / M) / (M - 1));
}

void dataset_apply_norm(Dataset *d, float mean, float std) {
    k_scale_shift((long)d->n * 784, 1 / std, -mean / std, d->x, d->x);
}

int *split_load(const char *path, int *n) {
    long size;
    int *idx = (int *)read_file(path, &size);
    *n = (int)(size / sizeof(int));
    return idx;
}

void cfg_default(Cfg *c) {
    memset(c, 0, sizeof(*c));
    c->n_blocks = 3, c->base_filters = 32, c->use_batchnorm = 1, c->batch_size = 128;
    c->onecycle = 1, c->lr = 1e-3f;
}

int cfg_load(Cfg *c, const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    char key[64];
    double v;
    while (fscanf(f, "%63s %lf", key, &v) == 2) {
        if (!strcmp(key, "n_blocks")) c->n_blocks = (int)v;
        else if (!strcmp(key, "base_filters")) c->base_filters = (int)v;
        else if (!strcmp(key, "use_batchnorm")) c->use_batchnorm = (int)v;
        else if (!strcmp(key, "gap")) c->gap = (int)v;
        else if (!strcmp(key, "batch_size")) c->batch_size = (int)v;
        else if (!strcmp(key, "onecycle")) c->onecycle = (int)v;
        else if (!strcmp(key, "dropout")) c->dropout = (float)v;
        else if (!strcmp(key, "weight_decay")) c->weight_decay = (float)v;
        else if (!strcmp(key, "label_smoothing")) c->label_smoothing = (float)v;
        else if (!strcmp(key, "lr")) c->lr = (float)v;
        else if (!strcmp(key, "rotation")) c->rotation = (float)v;
        else if (!strcmp(key, "translate")) c->translate = (float)v;
        else if (!strcmp(key, "scale_jitter")) c->scale_jitter = (float)v;
        else { fprintf(stderr, "unknown config key %s\n", key); return -1; }
    }
    fclose(f);
    return 0;
}

void cfg_save(const Cfg *c, const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); exit(1); }
    fprintf(f, "n_blocks %d\nbase_filters %d\nuse_batchnorm %d\ngap %d\nbatch_size %d\nonecycle %d\n",
            c->n_blocks, c->base_filters, c->use_batchnorm, c->gap, c->batch_size, c->onecycle);
    fprintf(f, "dropout %.9g\nweight_decay %.9g\nlabel_smoothing %.9g\nlr %.9g\n", c->dropout,
            c->weight_decay, c->label_smoothing, c->lr);
    fprintf(f, "rotation %.9g\ntranslate %.9g\nscale_jitter %.9g\n", c->rotation, c->translate,
            c->scale_jitter);
    fclose(f);
}
