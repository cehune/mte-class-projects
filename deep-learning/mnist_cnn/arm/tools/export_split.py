"""Writes the notebook's validation indices, so that the native trainer uses the same split."""
from pathlib import Path

import numpy as np
from torchvision.datasets import MNIST

SEED = 577
y = MNIST("../data", train=True, download=True).targets.numpy()

rng = np.random.default_rng(SEED)
val_idx = []
for digit in range(10):
    digit_idx = np.flatnonzero(y == digit)
    rng.shuffle(digit_idx)
    val_idx.append(digit_idx[: round(len(digit_idx) / 6)])
val_idx = np.sort(np.concatenate(val_idx)).astype(np.int32)

Path("data_split").mkdir(exist_ok=True)
val_idx.tofile("data_split/val_idx.bin")
print(len(val_idx), "validation indices")
