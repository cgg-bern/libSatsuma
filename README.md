# libSatsuma - C++ library for bidirected flow problems

## Introduction

Satsuma extends the LEMON graph library by data structures and algorithms
for handling network flow problems on general bidirected graphs.

You can find the necessary mathematical background and explanations in

**Min-Deviation-Flow in Bi-directed Graphs for T-Mesh Quantization**, Martin Heistermann, Jethro Warnett, David Bommes, ACM Trans. Graph. 2023

Please cite our paper if you use this library in academic context :-)

Project page: https://www.algohex.eu/publications/bimdf-quantization/

## Project structure

```
├── src/
│   └── libsatsuma/:     library source code
│       ├── Config/:     templates for cmake-generated headers
│       ├── Extra/:      High-level tools and serialization utilities
│       ├── Problems/:   Data structures to represent problem instances and solutions
│       ├── Reductions/: Algorithms to convert between different kinds of problems
│       └── Solvers/:    Solver Algorithms
├── tests/:              unit tests (TODO)
└── tools/
    └── example:         solves a simple fixed problem
    └── gen_figure:      solves a fixed problem and outputs intermediate reduction steps for visualization.
```

## Building and linking

### Dependencies

* [CMake](https://cmake.org) >= 3.18
* A recent C++ compiler
* LEMON
* (Optional) [Gurobi](https://gurobi.com) to solve Bi-MDF as ILP/IQP problem
* (Optional) [nlohmann-json](https://github.com/nlohmann/json/) for (de)serialization
* (Optional) Blossom-V as an alternative solver for Weighted Matching (not under free license)
* libTimekeeper (automatically downloaded)

### Commandline build

`mkdir build && cd build && cmake .. && cmake --build .`

This will automatically download all required dependencies.

## License

libSatsuma is available under the terms of the [MIT License](LICENSE).

***Note:*** By default, libSatsuma downloads, compiles and links to the non-free Blossom-V implementation. Set cmake option `SATSUMA_ENABLE_BLOSSOM5=OFF` to prevent this.


## Contact

Author: Martin Heistermann <martin.heistermann@unibe.ch>
