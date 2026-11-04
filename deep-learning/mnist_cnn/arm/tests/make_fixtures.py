"""Generates reference outputs from PyTorch. The C parity test must reproduce them.

File format, all little-endian: 13 float32 config values, then N, then x, y, the initial
parameters in order, logits, loss, gradients, parameters and BN stats after STEPS updates.
"""
import math
import struct
from pathlib import Path

import torch
import torch.nn as nn
import torch.nn.functional as F

OUT = Path("tests/fixtures")
OUT.mkdir(exist_ok=True)
STEPS = 5
TOTAL_STEPS = 20


class CNN(nn.Module):
    def __init__(self, cfg):
        super().__init__()
        layers = []
        in_ch = 1
        for block in range(cfg["n_blocks"]):
            out_ch = cfg["base_filters"] * 2**block
            for _ in range(2):
                layers.append(nn.Conv2d(in_ch, out_ch, 3, padding=1, bias=not cfg["use_batchnorm"]))
                if cfg["use_batchnorm"]:
                    layers.append(nn.BatchNorm2d(out_ch))
                layers.append(nn.ReLU())
                in_ch = out_ch
            layers.append(nn.MaxPool2d(2))
        self.features = nn.Sequential(*layers)
        size = 28 // 2 ** cfg["n_blocks"]
        self.gap = cfg["gap"]
        self.fc = nn.Linear(in_ch if self.gap else in_ch * size * size, 10)

    def forward(self, x):
        x = self.features(x)
        x = x.mean(dim=(2, 3)) if self.gap else x.flatten(1)
        return self.fc(x)


def write(f, tensor):
    f.write(tensor.detach().cpu().float().numpy().tobytes())


def params(model):
    return [p for p in model.parameters()]


def bn_stats(model):
    return [b for m in model.modules() if isinstance(m, nn.BatchNorm2d) for b in (m.running_mean, m.running_var)]


def make(name, cfg):
    torch.manual_seed(577)
    model = CNN(cfg)
    N = 16
    x = torch.randn(N, 1, 28, 28)
    y = torch.randint(0, 10, (N,))
    with open(OUT / f"{name}.bin", "wb") as f:
        f.write(struct.pack("<13f", cfg["n_blocks"], cfg["base_filters"], cfg["use_batchnorm"], cfg["gap"],
                            N, 1, 0, cfg["weight_decay"], cfg["label_smoothing"], cfg["lr"], 0, 0, 0))
        f.write(struct.pack("<i", N))
        write(f, x)
        f.write(y.numpy().astype("int32").tobytes())
        for p in params(model):
            write(f, p)
        opt = torch.optim.AdamW(model.parameters(), lr=cfg["lr"], weight_decay=cfg["weight_decay"])
        sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=cfg["lr"], total_steps=TOTAL_STEPS)
        for step in range(STEPS):
            logits = model(x)
            loss = F.cross_entropy(logits, y, label_smoothing=cfg["label_smoothing"])
            opt.zero_grad()
            loss.backward()
            if step == 0:
                write(f, logits)
                f.write(struct.pack("<f", loss.item()))
                for p in params(model):
                    write(f, p.grad)
            opt.step()
            sched.step()
        for p in params(model):
            write(f, p)
        for b in bn_stats(model):
            write(f, b)


make("bn_flatten", {"n_blocks": 2, "base_filters": 8, "use_batchnorm": 1, "gap": 0,
                    "weight_decay": 1e-3, "label_smoothing": 0.1, "lr": 1e-3})
make("plain_gap", {"n_blocks": 3, "base_filters": 4, "use_batchnorm": 0, "gap": 1,
                   "weight_decay": 1e-2, "label_smoothing": 0.0, "lr": 3e-3})

# grid_sample reference: normalized digit-like image with a background of -0.4241
torch.manual_seed(1)
bg = -0.4241
image = torch.full((1, 1, 28, 28), bg)
image[0, 0, 6:22, 10:18] = torch.rand(16, 8) * 3
angle = math.radians(13)
theta = torch.tensor([[[math.cos(angle) / 1.1, -math.sin(angle) / 1.1, 0.14],
                       [math.sin(angle) / 1.1, math.cos(angle) / 1.1, -0.2]]], dtype=torch.float32)
grid = F.affine_grid(theta, image.shape, align_corners=False)
warped = F.grid_sample(image - bg, grid, align_corners=False) + bg
with open(OUT / "augment.bin", "wb") as f:
    write(f, theta)
    write(f, image)
    write(f, warped)
print("fixtures written")
