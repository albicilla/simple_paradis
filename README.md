[![Build Status](https://travis-ci.org/deftio/travis-ci-cpp-example.svg?branch=master)](https://travis-ci.org/albicilla/simple_paradis)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# PARADIS
Implementation of PARADIS - fast parallel radix sort algorithm. Paper URL http://www.vldb.org/pvldb/vol8/p1518-cho.pdf

Using this code, you can easily evaluate the execution time of PARADIS.

PARADIS: An Efficient Parallel Algorithm for In-place Radix sort (Cho et al. VLDB 2015)

MSD radix sort

## Requirements
* CMake `>= 3.50`
* C++ Compiler `>= C++17`
* OpenMP
## How to build
```sh
mkdir build
cd build
cmake ..
make
```

## How to run
You can give the number of threads and number of data by command line arguments.
```sh
./paradis <number of threads> <number of keys>
```

Example
```sh
./paradis 64 100000000
```

## Reference

Cho, M., Brand, D., Bordawekar, R., Finkler, U., Kulandaisamy, V., & Puri, R. (2015). PARADIS: an efficient parallel algorithm for in-place radix sort. Proceedings of the VLDB Endowment, 8(12), 1518–1529. https://doi.org/10.14778/2824032.2824050


