/* A short run on real MNIST must learn */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nn.h"

int main(void) {
    Dataset all = idx_load("../data/MNIST/raw/train-images-idx3-ubyte", "../data/MNIST/raw/train-labels-idx1-ubyte");
    int n_val;
    int *val_idx = split_load("data_split/val_idx.bin", &n_val);
    Dataset val = dataset_subset(&all, val_idx, n_val);
    char *is_val = calloc(all.n, 1);
    for (int i = 0; i < n_val; i++) is_val[val_idx[i]] = 1;
    int *idx = malloc(4000 * sizeof(int)), n = 0;
    for (int i = 0; n < 4000; i++)
        if (!is_val[i]) idx[n++] = i;
    Dataset train = dataset_subset(&all, idx, n);
    if (is_val[idx[0]]) { printf("split overlap\n"); return 1; }

    float mean, std;
    dataset_norm(&train, &mean, &std);
    dataset_apply_norm(&train, mean, std);
    dataset_apply_norm(&val, mean, std);

    Cfg c;
    cfg_default(&c);
    c.n_blocks = 2, c.base_filters = 8, c.batch_size = 64, c.lr = 2e-3f;
    Net *net = net_build(&c, 512, 577);
    int epochs = 3;
    long total = epochs * ((n + 63) / 64), s = 0;
    for (int e = 0; e < epochs; e++)
        for (int i = 0; i < n; i += 64) {
            int N = n - i < 64 ? n - i : 64;
            float lr, b1;
            onecycle_at(c.lr, s++, total, &lr, &b1);
            net_forward(net, train.x + (long)i * 784, N, 1);
            net_zero_grad(net);
            net_loss_backward(net, train.y + i);
            opt_step(net, lr, b1);
        }
    double acc = net_accuracy(net, &val, NULL);
    printf("test_learn: validation accuracy %.4f: %s\n", acc, acc > 0.9 ? "ok" : "FAILED");
    return acc <= 0.9;
}
