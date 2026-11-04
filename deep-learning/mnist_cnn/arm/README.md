# MNIST CNN in ARM64 assembly

Native trainer for Apple Silicon. NEON kernels are in `src/kernels.S`, and C handles the layers, data and I/O. No PyTorch at runtime.

```bash
make test
make
./build/mnistarm train --config configs/best.txt --epochs 15 --run best
./build/mnistarm eval --run best
make viz
```

`make viz` serves the dashboard at http://localhost:8577. It updates live during training.

`./build/mnistarm eval --run best --test` reads the test set, so run it only after the model is final. `make test` needs the MNIST files, which the notebook downloads to `../data`.
