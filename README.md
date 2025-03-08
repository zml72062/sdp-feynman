# Computing Master Feynman Integrals with Semi-Definite Programming

## Introduction

This repository serves as a simple implementation for the paper [Feynman Integrals from Positivity Constraints](https://arxiv.org/pdf/2303.15624). See `intro.pdf` for a brief formal introduction.

## Getting Started

If all prerequisites are already satisfied, simply run
```sh
make
```
followed by
```sh
./master <config_file.yaml>
```
to launch master integral evaluation. Please check `examples/` subdirectory for configure file format.

The dependencies `gsl` and `sdpa` are optional. Without `gsl`, the Euclidean region check is disabled. To compile without `gsl`, run
```sh
make NO_GSL=true
```
Without `sdpa`, the program cannot actually evaluate master integrals. Instead, it outputs a semi-definite program in plain text. The output format is compatible with `sdpa`, so that one may process it later with an external `sdpa` solver to get the values of master integrals. To compile without `sdpa`, run
```sh
make NO_SDPA_LIB=true
```

## Prerequisites

The "Instructions" sections are tested on Ubuntu 22.04 LTS.

- **yaml-cpp** (https://github.com/jbeder/yaml-cpp)

- **Fermat 6.5** (http://home.bway.net/lewis/)
    - **Instructions:** Download binary version of Fermat 6.5 [here](http://home.bway.net/lewis/zip.html) and extract. Then set the environment variable `FERMATPATH` to the path to Fermat executable `fer64`, namely
        ```sh
        export FERMATPATH="/path/to/ferl6/fer64"
        ```

- **GiNaC** (https://www.ginac.de/), which itself depends on **CLN** (https://www.ginac.de/CLN/)
    - **Instructions:** To build CLN and GiNaC from source, clone repositories
        ```sh
        git clone git://www.ginac.de/cln.git
        git clone git://www.ginac.de/ginac.git
        ```
        and check the respective `INSTALL` files for installation instructions.


- **Kira 2.3** (https://gitlab.com/kira-pyred/kira)
    - **Instructions:** Clone Kira repository
        ```sh
        git clone https://gitlab.com/kira-pyred/kira.git
        ```
        and check `README.rst` for installation instructions.

- **GSL (GNU Scientific Library) 2.8 (optional)** (https://www.gnu.org/software/gsl/)

- **SDPA 7.3.15 (optional)** (https://sdpa.sourceforge.net/)
    - **Instructions:** The binary version of SDPA can be installed by `apt-get`,
        ```sh
        # To install SDPA,
        sudo apt-get install sdpa
        # To install SDPA callable libraries and examples,
        sudo apt-get install libsdpa-dev
        ```