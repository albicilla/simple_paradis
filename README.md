[![Build Status](https://travis-ci.org/deftio/travis-ci-cpp-example.svg?branch=master)](https://travis-ci.org/albicilla/simple-paradis)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# simple_paradis
implementation of PARADIS - fast parallel radix sort algorithm. http://www.vldb.org/pvldb/vol8/p1518-cho.pdf

Using this code, you can easily compare the execution time of std::sort and PARADIS.
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



