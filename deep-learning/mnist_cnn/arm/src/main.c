/* CLI: train, eval, predict, bench. Train writes JSON files that the dashboard reads. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include "nn.h"

#define RAW "../data/MNIST/raw/"
#define SPLIT "data_split/val_idx.bin"
#define EVAL_BATCH 512

static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

static const char *arg(int argc, char **argv, const char *key, const char *fallback) {
    for (int i = 2; i + 1 < argc; i++)
        if (!strcmp(argv[i], key)) return argv[i + 1];
    return fallback;
}

static void load_split(Dataset *train, Dataset *val) {
    Dataset all = idx_load(RAW "train-images-idx3-ubyte", RAW "train-labels-idx1-ubyte");
    int n_val;
    int *val_idx = split_load(SPLIT, &n_val);
    char *is_val = calloc(all.n, 1);
    int *train_idx = malloc(all.n * sizeof(int)), n_train = 0;
    for (int i = 0; i < n_val; i++) is_val[val_idx[i]] = 1;
    for (int i = 0; i < all.n; i++)
        if (!is_val[i]) train_idx[n_train++] = i;
    *train = dataset_subset(&all, train_idx, n_train);
    *val = dataset_subset(&all, val_idx, n_val);
}

static void json_floats(FILE *f, const float *x, long n, float scale, float shift) {
    fputc('[', f);
    for (long i = 0; i < n; i++) fprintf(f, "%s%.3f", i ? "," : "", x[i] * scale + shift);
    fputc(']', f);
}

static void json_ints(FILE *f, const int *x, long n) {
    fputc('[', f);
    for (long i = 0; i < n; i++) fprintf(f, "%s%d", i ? "," : "", x[i]);
    fputc(']', f);
}

/* Writes via a temp file, so the dashboard never reads a partial snapshot */
static void write_snapshot(Net *net, const Dataset *val, const char *dir, int epoch, double acc,
                           const int *confusion, int *errors, int n_errors, long total_epochs) {
    char tmp[512], path[512];
    snprintf(tmp, sizeof tmp, "%s/snapshot.tmp", dir);
    snprintf(path, sizeof path, "%s/snapshot.json", dir);
    FILE *f = fopen(tmp, "w");
    if (!f) return;
    fprintf(f, "{\"epoch\":%d,\"epochs\":%ld,\"val_acc\":%.5f,\"confusion\":", epoch, total_epochs, acc);
    json_ints(f, confusion, 100);
    /* First-layer filters */
    Layer *first = &net->l[0];
    int show = first->co < 32 ? first->co : 32;
    fputs(",\"filters\":[", f);
    for (int o = 0; o < show; o++) {
        fputs(o ? "," : "", f);
        json_floats(f, net->p[first->p0].w + (long)o * 9, 9, 1, 0);
    }
    /* Activations per pooled block, for one validation image */
    const float *z = net_forward(net, val->x, 1, 0);
    (void)z;
    fputs("],\"input\":", f);
    json_floats(f, val->x, 784, net->std, net->mean);
    fprintf(f, ",\"input_label\":%d,\"maps\":[", val->y[0]);
    int first_map = 1;
    for (int i = 0; i < net->n_layers; i++) {
        Layer *l = &net->l[i];
        if (l->type != L_POOL) continue;
        fprintf(f, "%s{\"size\":%d,\"channels\":[", first_map ? "" : ",", l->ho);
        int chans = l->co < 16 ? l->co : 16;
        for (int c = 0; c < chans; c++) {
            fputs(c ? "," : "", f);
            json_floats(f, l->y + (long)c * l->ho * l->wo, (long)l->ho * l->wo, 1, 0);
        }
        fputs("]}", f);
        first_map = 0;
    }
    /* Misclassified validation digits */
    fputs("],\"errors\":[", f);
    for (int i = 0; i < n_errors; i++) {
        int idx = errors[2 * i];
        fprintf(f, "%s{\"true\":%d,\"pred\":%d,\"image\":", i ? "," : "", val->y[idx], errors[2 * i + 1]);
        json_floats(f, val->x + (long)idx * 784, 784, net->std, net->mean);
        fputc('}', f);
    }
    fputs("]}\n", f);
    fclose(f);
    rename(tmp, path);
}

static int collect_errors(Net *net, const Dataset *d, int *errors, int max) {
    int found = 0;
    for (int i = 0; i < d->n && found < max; i += EVAL_BATCH) {
        int N = d->n - i < EVAL_BATCH ? d->n - i : EVAL_BATCH;
        const float *z = net_forward(net, d->x + (long)i * 784, N, 0);
        for (int n = 0; n < N && found < max; n++) {
            int arg = 0;
            for (int k = 1; k < 10; k++)
                if (z[n * 10 + k] > z[n * 10 + arg]) arg = k;
            if (arg != d->y[i + n]) errors[2 * found] = i + n, errors[2 * found + 1] = arg, found++;
        }
    }
    return found;
}

static int cmd_train(int argc, char **argv) {
    Cfg c;
    cfg_default(&c);
    const char *cfg_path = arg(argc, argv, "--config", NULL);
    if (cfg_path && cfg_load(&c, cfg_path)) { fprintf(stderr, "cannot read config\n"); return 1; }
    int epochs = atoi(arg(argc, argv, "--epochs", "30"));
    int seed = atoi(arg(argc, argv, "--seed", "577"));
    const char *name = arg(argc, argv, "--run", "default");
    int limit = atoi(arg(argc, argv, "--limit", "0"));

    char dir[256], path[512];
    mkdir("runs", 0755);
    snprintf(dir, sizeof dir, "runs/%s", name);
    mkdir(dir, 0755);
    snprintf(path, sizeof path, "%s/config.txt", dir);
    cfg_save(&c, path);
    snprintf(path, sizeof path, "%s/metrics.jsonl", dir);
    FILE *log = fopen(path, "w");

    Dataset train, val;
    load_split(&train, &val);
    if (limit && limit < train.n) train.n = limit;
    float mean, std;
    dataset_norm(&train, &mean, &std);
    dataset_apply_norm(&train, mean, std);
    dataset_apply_norm(&val, mean, std);
    float bg = -mean / std;

    int max_batch = c.batch_size > EVAL_BATCH ? c.batch_size : EVAL_BATCH;
    Net *net = net_build(&c, max_batch, seed);
    net->mean = mean, net->std = std;
    fprintf(stderr, "train %d, val %d, parameters %ld\n", train.n, val.n, net_param_count(net));

    int *order = malloc(train.n * sizeof(int));
    for (int i = 0; i < train.n; i++) order[i] = i;
    float *batch = malloc((long)c.batch_size * 784 * sizeof(float));
    float *warped = malloc((long)c.batch_size * 784 * sizeof(float));
    int *labels = malloc(c.batch_size * sizeof(int));
    long per_epoch = (train.n + c.batch_size - 1) / c.batch_size, total = per_epoch * epochs, s = 0;
    double best = 0, start = now();
    Rng aug = {(uint64_t)seed * 31 + 7};
    int errors[48];

    for (int e = 0; e < epochs; e++) {
        rng_shuffle(&net->rng, order, train.n);
        double epoch_loss = 0, t0 = now();
        int epoch_correct = 0;
        float lr = c.lr, b1 = 0.9f;
        for (int i = 0; i < train.n; i += c.batch_size) {
            int N = train.n - i < c.batch_size ? train.n - i : c.batch_size;
            for (int n = 0; n < N; n++) {
                memcpy(batch + (long)n * 784, train.x + (long)order[i + n] * 784, 784 * sizeof(float));
                labels[n] = train.y[order[i + n]];
            }
            augment_batch(&c, &aug, batch, warped, N, bg);
            if (c.onecycle) onecycle_at(c.lr, s, total, &lr, &b1);
            net_forward(net, warped, N, 1);
            net_zero_grad(net);
            net_loss_backward(net, labels);
            opt_step(net, lr, b1);
            epoch_loss += net->loss * N;
            epoch_correct += net->correct;
            if (++s % 50 == 0) {
                fprintf(log, "{\"kind\":\"step\",\"step\":%ld,\"total\":%ld,\"loss\":%.5f,\"lr\":%.6f,\"img_s\":%.0f}\n",
                        s, total, net->loss, lr, (i + N) / (now() - t0));
                fflush(log);
            }
        }
        int confusion[100] = {0};
        double acc = net_accuracy(net, &val, confusion);
        double secs = now() - t0;
        fprintf(log, "{\"kind\":\"epoch\",\"epoch\":%d,\"train_loss\":%.5f,\"train_acc\":%.5f,\"val_acc\":%.5f,\"seconds\":%.2f}\n",
                e + 1, epoch_loss / train.n, (double)epoch_correct / train.n, acc, secs);
        fflush(log);
        fprintf(stderr, "epoch %d/%d loss %.4f train %.4f val %.4f (%.1fs)\n", e + 1, epochs,
                epoch_loss / train.n, (double)epoch_correct / train.n, acc, secs);
        if (acc >= best) best = acc;
        snprintf(path, sizeof path, "%s/model.bin", dir);
        net_save(net, path);
        write_snapshot(net, &val, dir, e + 1, acc, confusion, errors, collect_errors(net, &val, errors, 24), epochs);
    }
    fprintf(stderr, "done in %.0fs, final val %.4f, best val %.4f\n", now() - start,
            net_accuracy(net, &val, NULL), best);
    fclose(log);
    return 0;
}

static Net *load_run(const char *name) {
    char path[512];
    Cfg c;
    cfg_default(&c);
    snprintf(path, sizeof path, "runs/%s/config.txt", name);
    if (cfg_load(&c, path)) { fprintf(stderr, "cannot read %s\n", path); exit(1); }
    Net *net = net_build(&c, EVAL_BATCH, 1);
    snprintf(path, sizeof path, "runs/%s/model.bin", name);
    if (net_load(net, path)) { fprintf(stderr, "cannot read %s\n", path); exit(1); }
    return net;
}

/* The test set is only ever read here, and only with an explicit --test */
static int cmd_eval(int argc, char **argv) {
    const char *name = arg(argc, argv, "--run", "default");
    Net *net = load_run(name);
    int use_test = 0;
    for (int i = 2; i < argc; i++) use_test |= !strcmp(argv[i], "--test");
    Dataset d, unused;
    if (use_test) d = idx_load(RAW "t10k-images-idx3-ubyte", RAW "t10k-labels-idx1-ubyte");
    else load_split(&unused, &d);
    dataset_apply_norm(&d, net->mean, net->std);
    int confusion[100] = {0};
    double acc = net_accuracy(net, &d, confusion);
    printf("%s accuracy: %.4f (%d images)\n", use_test ? "test" : "validation", acc, d.n);
    if (use_test) {
        char path[512];
        snprintf(path, sizeof path, "runs/%s/test.json", name);
        FILE *f = fopen(path, "w");
        fprintf(f, "{\"test_acc\":%.5f,\"confusion\":", acc);
        json_ints(f, confusion, 100);
        fputs("}\n", f);
        fclose(f);
    }
    return 0;
}

/* Reads 784 floats in [0, 1] from stdin and prints class probabilities as JSON */
static int cmd_predict(int argc, char **argv) {
    Net *net = load_run(arg(argc, argv, "--run", "default"));
    float x[784];
    for (int i = 0; i < 784; i++) {
        if (scanf("%f", &x[i]) != 1) { fprintf(stderr, "need 784 values\n"); return 1; }
        x[i] = (x[i] - net->mean) / net->std;
    }
    const float *z = net_forward(net, x, 1, 0);
    double top = z[0], sum = 0, p[10];
    for (int k = 1; k < 10; k++) if (z[k] > top) top = z[k];
    for (int k = 0; k < 10; k++) sum += p[k] = exp(z[k] - top);
    printf("{\"probs\":[");
    for (int k = 0; k < 10; k++) printf("%s%.5f", k ? "," : "", p[k] / sum);
    printf("]}\n");
    return 0;
}

static double bench_conv(void (*fn)(long, float *, const float *, const float *, const float *, const float *)) {
    float rows[3][32], dst[28], w[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    for (int r = 0; r < 3; r++)
        for (int i = 0; i < 32; i++) rows[r][i] = 0.01f * i;
    double t = now();
    for (int i = 0; i < 2000000; i++) {
        memset(dst, 0, sizeof dst);
        fn(28, dst, rows[0] + 1, rows[1] + 1, rows[2] + 1, w);
        __asm__ volatile("" : : : "memory");
    }
    return 2000000 * 28 * 9 / (now() - t) / 1e9;
}

static int cmd_bench(void) {
    double a = bench_conv(k_conv3_row), c = bench_conv(r_conv3_row);
    printf("{\"asm_gmac_s\":%.2f,\"c_gmac_s\":%.2f}\n", a, c);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: mnistarm train|eval|predict|bench [--config f] [--run name] [--epochs n] [--seed n] [--test]\n");
        return 2;
    }
    if (!strcmp(argv[1], "train")) return cmd_train(argc, argv);
    if (!strcmp(argv[1], "eval")) return cmd_eval(argc, argv);
    if (!strcmp(argv[1], "predict")) return cmd_predict(argc, argv);
    if (!strcmp(argv[1], "bench")) return cmd_bench();
    fprintf(stderr, "unknown command %s\n", argv[1]);
    return 2;
}
