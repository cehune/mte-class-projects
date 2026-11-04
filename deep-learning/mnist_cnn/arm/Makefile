CC = /usr/bin/clang
# -fno-slp-vectorize and -fno-vectorize keep kernels_ref.c as a fair scalar comparison for bench
CFLAGS = -arch arm64 -O3 -std=gnu11 -Wall -Wextra -Wno-unused-parameter -Iinclude
PY = ../.venv/bin/python
CORE = build/kernels.o build/kernels_ref.o build/model.o build/optim.o build/data.o build/augment.o
TESTS = build/test_kernels build/test_model build/test_parity build/test_learn

all: build/mnistarm

build/%.o: src/%.S | build
	$(CC) -arch arm64 -c $< -o $@
build/kernels_ref.o: src/kernels_ref.c include/nn.h | build
	$(CC) $(CFLAGS) -fno-vectorize -fno-slp-vectorize -c $< -o $@
build/%.o: src/%.c include/nn.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/mnistarm: $(CORE) build/main.o
	$(CC) $(CFLAGS) $^ -o $@
build/test_%: tests/test_%.c $(CORE)
	$(CC) $(CFLAGS) $^ -o $@
build:
	mkdir -p build

data_split/val_idx.bin: tools/export_split.py
	$(PY) tools/export_split.py
tests/fixtures/.stamp: tests/make_fixtures.py
	$(PY) tests/make_fixtures.py
	touch $@

test: $(TESTS) data_split/val_idx.bin tests/fixtures/.stamp
	./build/test_kernels
	./build/test_model
	./build/test_parity
	./build/test_learn

bench: build/mnistarm
	./build/mnistarm bench

viz:
	$(PY) viz/server.py

clean:
	rm -rf build

.PHONY: all test bench viz clean
