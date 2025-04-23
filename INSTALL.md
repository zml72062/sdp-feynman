# Notes on installing Kira for a non-root Linux user

1. Install GNU `m4`.
    * Fetch source code of `m4` [here](https://ftp.gnu.org/gnu/m4/), namely
        ```sh
        wget https://ftp.gnu.org/gnu/m4/m4-latest.tar.gz
        tar xvf m4-latest.tar.gz
        ```
    * At `m4` source directory, execute
        ```sh
        ./configure --prefix=/install/path
        make && make install
        ```
    **Notice:** Do not forget to add `/install/path/bin` and `/install/path/lib` to `$PATH`.
2. Install GNU `autoconf`.
    * Fetch source code of `autoconf` [here](https://ftp.gnu.org/gnu/autoconf/), namely
        ```sh
        wget https://ftp.gnu.org/gnu/autoconf/autoconf-latest.tar.gz
        tar xvf autoconf-latest.tar.gz
        ```
    * At `autoconf` source directory, execute
        ```sh
        ./configure --prefix=/install/path
        make && make install
        ```
3. Install GNU `automake`.
    * Fetch source code of `automake` [here](https://ftp.gnu.org/gnu/automake/), namely
        ```sh
        wget https://ftp.gnu.org/gnu/automake/automake-1.17.tar.gz
        tar xvf automake-1.17.tar.gz
        ```
    * At `automake` source directory, execute
        ```sh
        ./configure --prefix=/install/path
        make && make install
        ```
4. Install GNU `libtool`.
    * Fetch source code of `libtool` [here](https://www.gnu.org/software/libtool/), namely
        ```sh
        wget https://ftpmirror.gnu.org/libtool/libtool-2.5.4.tar.gz
        tar xvf libtool-2.5.4.tar.gz
        ```
    * At `libtool` source directory, execute
        ```sh
        ./configure --prefix=/install/path
        make && make install
        ```
5. Install GNU `texinfo`.
    * Fetch source code of `texinfo` [here](https://ftp.gnu.org/gnu/texinfo/), namely
        ```sh
        wget https://ftp.gnu.org/gnu/texinfo/texinfo-7.2.tar.gz
        tar xvf texinfo-7.2.tar.gz
        ```
    * At `texinfo` source directory, execute
        ```sh
        ./configure --prefix=/install/path
        make && make install
        ```
6. Install GNU `gmp`.
    * Fetch source code of `gmp` [here](https://gmplib.org/), namely
        ```sh
        wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
        tar xvf gmp-6.3.0.tar.xz
        ```
    * At `gmp` source directory, execute
        ```sh
        ./configure --prefix=/install/path
        make && make install
        ```
7. Install GNU `mpfr`.
    * Fetch source code of `mpfr` by
        ```sh
        git clone https://gitlab.inria.fr/mpfr/mpfr.git
        ```
    * At `mpfr` source directory, execute
        ```sh
        ./autogen.sh
        ./configure --prefix=/install/path CPPFLAGS=-I/install/path/include LDFLAGS=-L/install/path/lib
        make && make install
        ```
8. Install `flint`.
    * Fetch source code of `flint` [here](https://flintlib.org/downloads.html), namely
        ```sh
        wget https://flintlib.org/download/flint-3.2.1.tar.gz
        tar xvf flint-3.2.1.tar.gz
        ```
    * At `flint` source directory, execute
        ```sh
        ./configure --prefix=/install/path
        make && make install
        ```
9. Install `cln`.
    * Fetch source code of `cln` by
        ```sh
        git clone git://www.ginac.de/cln.git
        ```
    * At `cln` source directory, execute
        ```sh
        mkdir cln_build
        cd cln_build
        cmake -DCMAKE_INSTALL_PREFIX:PATH=/install/path -GNinja ..
        cmake --build .
        cmake --build . -t install
        ```
    **Notice:** One should ensure that CMake and Ninja are installed on the machine.
10. Install `ginac`.
    * Fetch source code of `ginac` by
        ```sh
        git clone git://www.ginac.de/ginac.git
        ```
    * At `ginac` source directory, execute
        ```sh
        mkdir ginac_build
        cd ginac_build
        cmake -DCMAKE_INSTALL_PREFIX:PATH=/install/path ..
        make && make install
        ```
11. Install `firefly`.
    * Fetch source code of `firefly` by
        ```sh
        git clone https://gitlab.com/firefly-library/firefly.git
        ```
    * At `firefly` source directory, execute
        ```sh
        git checkout kira-2
        mkdir build
        cd build
        cmake -DCMAKE_INSTALL_PREFIX:PATH=/install/path -DWITH_FLINT=true ..
        make && make install
        ```
12. Install Fermat.
    * Fetch binary version of Fermat [here](http://home.bway.net/lewis/zip.html), namely
        ```sh
        wget http://home.bway.net/lewis/fermat64/ferl6.tar.gz
        tar xvf ferl6.tar.gz
        ```
    * Set environment variable `$FERMATPATH` to the path to `fer64` executable.
13. Install `kira`.
    * Fetch source code of `kira` by
        ```sh
        git clone https://gitlab.com/kira-pyred/kira
        ```
    * Install `meson` if it is not already installed on the machine,
        ```sh
        pip3 install --user meson
        ```
        Export path to `meson` executable to `$PATH` if necessary.
    * At `kira` source directory, execute
        ```sh
        meson setup --prefix=/install/path --pkg-config-path=/install/path/lib/pkgconfig builddir
        cd builddir
        ninja && ninja install
        ```

## Additional note on installing SDPA-DD, SDPA-QD and SDPA-GMP for high-precision semi-definite programming computation

We thank [Dr. Mao Zeng](https://gitlab.com/zengmao) (who is author of the paper [Feynman Integrals from Positivity Constraints](https://arxiv.org/pdf/2303.15624)) for providing this note.

SDPA-DD, SDPA-QD and SDPA-GMP are high-precision variants of the [SDPA](https://sdpa.sourceforge.net/) family of semi-definite programming solvers. The suffixes DD and QD stand for double-double precision and quadruple-double precision respectively.

SDPA-DD, SDPA-QD and SDPA-GMP are not so easy to install, especially when one has no root access to the computer. Here is the trick to install them.

1. Install Julia 1.8 from the tarball:
```sh
wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.5-linux-x86_64.tar.gz
tar xvf julia-1.8.5-linux-x86_64.tar.gz
```
One may find `julia` executable at directory `julia-1.8.5/bin`.

2. Install the `SDPAFamily.jl` package using the following commands in Julia:
```julia
using Pkg
Pkg.add("SDPAFamily")
```
Julia's package manager will automatically download and compile SDPA-DD, SDPA-QD and SDPA-GMP. 

3. Run the following commands in Julia:
```julia
using SDPAFamily
SDPAFamily.sdpa_dd
SDPAFamily.sdpa_qd
SDPAFamily.sdpa_gmp
```
This will print out the path of SDPA-DD, SDPA-QD and SDPA-GMP executables. Now one can quit Julia and use the binaries as standalone programs.

    