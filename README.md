[![Build Status](https://travis-ci.org/deftio/travis-ci-cpp-example.svg?branch=master)](https://travis-ci.org/albicilla/simple_paradis)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# simple_paradis
Implementation of PARADIS - fast parallel radix sort algorithm. Paper URL http://www.vldb.org/pvldb/vol8/p1518-cho.pdf

Using this code, you can easily compare the execution time of std::sort and PARADIS.

PARADIS:An Efficient Parallel Algorithm for In-place Radix sort (cho et al. VLDB 2015)の実装です。

MSD radix sortとして実装してあります。

number of data <= 2^32
number of threads には thread::hardware_concurrency() 以下の値を推奨します。

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
./paradis_ompf_repairgen <number of threads> <number of data>
```

Example
```sh
./paradis_ompf_repairgen 4 100000000
```

Example result (In my environment)
```sh
creating dataset... finish!

std::sort() is running... finish!
std::sort time 9577.141000[ms]

PARADIS is running... finish!
paradis time 2015.222000[ms]
```

## 参考文献
Cho, M., Brand, D., Bordawekar, R., Finkler, U., Kulandaisamy, V., & Puri, R. (2015). PARADIS: an efficient parallel algorithm for in-place radix sort. Proceedings of the VLDB Endowment, 8(12), 1518–1529. https://doi.org/10.14778/2824032.2824050


