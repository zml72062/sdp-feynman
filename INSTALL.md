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
6. Install GNU `mpfr`.
    * Fetch source code of `mpfr` by
        ```sh
        git clone https://gitlab.inria.fr/mpfr/mpfr.git
        ```
    * At `mpfr` source directory, execute
        ```sh
        ./autogen.sh
        ./configure --prefix=/install/path
        make && make install
        ```
7. Install `flint`.
    * Fetch source code of `flint` [here](https://flintlib.org/downloads.html), namely
        ```sh
        wget https://flintlib.org/download/flint-2.8.4.tar.gz
        tar xvf flint-2.8.4.tar.gz
        ```
    * At `flint` source directory, execute
        ```sh
        ./configure --prefix=/install/path
        make && make install
        ```
8. Install `cln`.
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
9. Install `ginac`.
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
10. Install `firefly`.
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
11. Install Fermat.
    * Fetch binary version of Fermat [here](http://home.bway.net/lewis/zip.html), namely
        ```sh
        wget http://home.bway.net/lewis/fermat64/ferl6.tar.gz
        tar xvf ferl6.tar.gz
        ```
    * Set environment variable `$FERMATPATH` to the path to `fer64` executable.
12. Install `kira`.
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

    