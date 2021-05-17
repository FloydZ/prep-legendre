Introduction
----
This repository contains the benchmark code for the paper [Legendre PRF (Multiple) Key Attacks and the Power of Preprocessin] (TODO). 

Requirements
----
- cmake 3.20
- c++ 20 compiler
- gmp, gmpxx
 
Build Instructions:
-----
To build the project simply run:
```bash
mkdir cmake-build-release
cd cmake-build-release/
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j2
```

Run Benchmark
----
To run a benchmark, simply execute
```bash
python gen.py -r 33 --iterations 10 --bench --executable code --include main.h
```
or 
```bash
python gen.py -r 33 --iterations 10 --bench --executable code_noprep --include main.h
```
to generate data for all prime up the bit length `32`, where each data point is averaged over 10 measurements. 

Not that you have the set the `m` value in the python script `gen.py` in line 94 depending on which setting (single-key or multiple-key) you want to benchmark.

The resulting `*.log` file can be converted into an `tikz` readable format with python script `plot.py`.
