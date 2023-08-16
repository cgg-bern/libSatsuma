# libSatsuma - C++ library for bidirected flow problems

## Introduction

Project page: https://www.algohex.eu/publications/bimdf-quantization/

**Satsuma extends the LEMON graph library by data structures and algorithms
for handling network flow problems on general bidirected graphs.**

Most importantly, it solves the *Minimum Deviation Flow* problem in bidirected
graphs (*Bi-MDF*) -- a generalization of the *Minimum Cost Flow* problem.
We provide an approximate and an exact solver that utilized iterative refinement.

You can find the necessary mathematical background and explanations in

[**Min-Deviation-Flow in Bi-directed Graphs for T-Mesh Quantization**](https://www.algohex.eu/publications/bimdf-quantization/), **Martin Heistermann, Jethro Warnett, David Bommes**, ACM Transactions on Graphics, Volume 42, Issue 4, 2023

Please cite our paper if you use this library in academic context :-)

## Applications

We demonstrate libSatsuma as an alternative to the commercial Gurobi solver in a recent quad remeshing (or *retopology*, in 3D artist parlance) tool, [QuadWild](https://www.quadmesh.cloud/), in our [fork of the code](https://github.com/cgg-bern/quadwild-bimdf/).


## BibTeX
```
@article{Heistermann:2023:BiMDF,
  author = {Heistermann, Martin and Warnett, Jethro and Bommes, David}
  title = {Min-Deviation-Flow in Bi-directed Graphs for T-Mesh Quantization},
  journal = {ACM Trans. Graph.},
  volume = {42},
  number = {4},
  year = {2023},
  publisher = {ACM},
  address = {New York, NY, USA},
  doi = {10.1145/3592437}
}
```

## Minimum-Deviation-Flow Problem in Bidirected Graphs (Bi-MDF)

Find $x\in \mathbb{Z}^n$ that minimizes $\displaystyle\sum_{1\leq i\leq n} c_i(x_i)$, such that

* $C\cdot x = b$, where
    * $b$ is the *demand* vector $b\in\mathbb{Z}^m$, and
    * $C \in \\\{-2,-1,0,1,2\\\}^{m\times n}$ is the node-edge incidence matrix of a bidirected graph,
      i.e., the sum of absolute values of each column of $C$ is not larger than 2,
* $\forall i: l_i \leq x_i \leq u_i$, where
    $l, u \in \mathbb{Z}^n$ are lower and upper bounds (upper bounds may be infinite)

$c_i: \mathbb{Z}\to\mathbb{R}$ are convex cost functions.
Cost functions $c_i$ with $u_i=\infty$ must attain a minimal value in the interval $[l_i, \infty)$, e.g., $1/x$ would only be allowed with a finite upper bound.




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

### Prerequisites and dependencies

* [CMake](https://cmake.org) >= 3.18
* A recent C++ compiler
* LEMON (automatically downloaded)
* libTimekeeper (automatically downloaded)
* (Optional) Blossom-V as an alternative solver for Weighted Matching (not under free license)
* (Optional) [nlohmann-json](https://github.com/nlohmann/json/) for (de)serialization
* (Optional) [Gurobi](https://gurobi.com) to solve Bi-MDF as ILP/IQP problem

#### Linux (Debian/Ubuntu)

If you don't have a development environment with CMake, run the following command:
```
sudo apt install cmake build-essential
```

#### MacOS

We recommend using [Homebrew](https://brew.sh/), which will automatically install the *XCode Command Line Tools* and allow you to install CMake using `brew install cmake`

#### Windows

You can use Visual Studio for a native build. With recent VS versions, CMake support is build-in, so you can open the `CMakeLists.txt` directly.

Alternatively, you can obtain a convenient Linux based development environment using the
[Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install)

### Commandline build

To generate a build folder and compile libSatsuma, run
```
cd path/to/libSatsuma
cmake . -B build      # generate build files in build/
cmake --build build   # compile
```

This will automatically download all required dependencies and compile libSatsuma as static library, as well as the example programs.

You can now run the example program
```
./build/tools/example/example
```
(The path may differ depending on the CMake generator used)

***Note:*** To download and use the Blossom-V library, set the cmake option `SATSUMA_ENABLE_BLOSSOM5=ON`: `cmake build -DSATSUMA_ENABLE_BLOSSOM5=ON`.


### Using libSatsuma in CMake-based projects

Our recommended approach is using `add_subdirectory` on a sub-folder (possibly git submodule) containing libSatsuma,
or automatically downloading libSatsuma:

```
include(FetchContent)
if (NOT TARGET satsuma::satsuma)
    FetchContent_Declare(satsuma
        GIT_REPOSITORY https://github.com/cgg-bern/libsatsuma
        GIT_TAG main # Recommended specifing a commit or tag instead of the `main` branch
        SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/satsuma"
        )
    FetchContent_MakeAvailable(satsuma)
endif()
```

libSatsuma will now be available as a target `satsuma::satsuma`, which you can link to your program using

```
target_link_libraries(your_target PRIVATE satsuma::satsuma)
```

By default, libSatsuma is built as static library. You can set the global
[BUILD_SHARED_LIBS](https://cmake.org/cmake/help/latest/variable/BUILD_SHARED_LIBS.html)
variable to true if you prefer a shared library.

## License

libSatsuma is available under the terms of the [MIT License](LICENSE).


## Contact

I'd be happy to receive any feedback from you via github issues or email.
Please let me know if you have any issues compiling or using libSatsuma.

Author: Martin Heistermann <martin.heistermann@unibe.ch>


## Acknowledgements

This project has been developed as part of the [AlgoHex](https://www.algohex.eu/) project, which has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme ([Grant agreement No. 853343](https://cordis.europa.eu/project/id/853343)).
