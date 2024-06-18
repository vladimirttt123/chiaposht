# Chia Proof of Space
![Build](https://github.com/Chia-Network/chiapos/actions/workflows/build-test-cplusplus.yml/badge.svg)
![Wheels](https://github.com/Chia-Network/chiapos/actions/workflows/build-wheels.yml/badge.svg)
![PyPI](https://img.shields.io/pypi/v/chiapos?logo=pypi)
![PyPI - Format](https://img.shields.io/pypi/format/chiapos?logo=pypi)
![GitHub](https://img.shields.io/github/license/Chia-Network/chiapos?logo=Github)
[![Coverage Status](https://coveralls.io/repos/github/Chia-Network/chiapos/badge.svg?branch=main)](https://coveralls.io/github/Chia-Network/chiapos?branch=main)
## Compression
### Abstract
This branch was build to provide some level of compression to exists plots of any size. 
This achived by realligning data inside plot and/or removing some data from table 1.

### Usage
Get the branch from git
```bash
git clone https://github.com/vladimirttt123/chiaposht.git -b compression
```
Build as described bellow for original project
```bash
# Requires cmake 3.14+

mkdir -p build && cd build
cmake ../
cmake --build . -- -j 6
```

Now you have executable ProofOfSpace that could be used to compress and library usually with name
chiapos.cpython-312-x86_64-linux-gnu.so (the axact name depends on python version and architecture)
In order to compress use
```bash
./ProofOfSpace compress <level> <original.plot> -f <output.plot>
```
Level is a positive number translated to number of bits to remove from table 1 of the plot. 
Level 0 does not remove any data just reallign data inside to save around 1.2-1.3% of sapce.

Now you can check timing of proof generation for plot with current compression on your system
```bash
./ProofOfSpace check 10 -f <plot.file>
```
To use compressed plots for farming you need replace the mentioned before 
library (chiapos.cpython-312-x86_64-linux-gnu.so) in original chia instalation.
The usual place for it 
```bash
<chia_install_folder>/venv/lib64/python3.12/site-packages
```
Before replacing chia services should be stopped like this from CLI
```bash
chia stop -d all
```
Than replace and run chia back as usual. The new library should support all exists plots as usual.
### Trableshooting
If compressed plots do not work, enable support of compressed plots in chia config file.

I do not have compressed plots of other type than I havn't check of other compressed plots 
will work after replacing the library but since the branch is fork of original chiapos with
very minimal changes in original code I hope it will work.

The version of chia client I checked on is 2.3.1. I think it should work also for 2.3 
and sure will not work for less than 2.1.

### Donate
If you thinks this compression is helpfull for you please consider to donate. My chia wallet is
```bash
xch1ch6s3q0enuj9wtemn473gkkvj0u8vlggypr375mk547e7aa48hmsql74e8
```

-------------------------------------------------------------------------------------------
#

Chia's proof of space is written in C++. Includes a plotter, prover, and
verifier. It exclusively runs on 64 bit architectures. Read the
[Proof of Space document](https://www.chia.net/wp-content/uploads/2022/09/Chia_Proof_of_Space_Construction_v1.1.pdf) to
learn about what proof of space is and how it works.

## C++ Usage Instructions

### Compile

```bash
# Requires cmake 3.14+

mkdir -p build && cd build
cmake ../
cmake --build . -- -j 6
```

## Static Compilation With glibc
### Statically compile ProofOfSpace
```bash
mkdir -p build && cd build
cmake -DBUILD_PROOF_OF_SPACE_STATICALLY=ON ../
cmake --build . -- -j 6
```

### Run tests

```bash
./RunTests
```

### CLI usage

```bash
./ProofOfSpace -k 25 -f "plot.dat" -m "0x1234" create
./ProofOfSpace -k 25 -f "final-plot.dat" -m "0x4567" -t TMPDIR -2 SECOND_TMPDIR create
./ProofOfSpace -f "plot.dat" prove <32 byte hex challenge>
./ProofOfSpace -k 25 verify <hex proof> <32 byte hex challenge>
./ProofOfSpace -f "plot.dat" check <iterations>
```

### Benchmark

```bash
time ./ProofOfSpace -k 25 create
```


### Hellman Attacks usage

There is an experimental implementation which implements some of the Hellman
Attacks that can provide significant space savings for the final file.


```bash
./HellmanAttacks -k 18 -f "plot.dat" -m "0x1234" create
./HellmanAttacks -f "plot.dat" check <iterations>
```

## Python

Finally, python bindings are provided in the python-bindings directory.

### Install

```bash
python3 -m venv .venv
. .venv/bin/activate
pip3 install .
```

### Run python tests

Testings uses pytest. Linting uses flake8 and mypy.

```bash
py.test ./tests -s -v
```

## ci Building
The primary build process for this repository is to use GitHub Actions to
build binary wheels for MacOS, Linux (x64 and aarch64), and Windows and publish
them with a source wheel on PyPi. See `.github/workflows/build.yml`. CMake uses
[FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html)
to download [pybind11](https://github.com/pybind/pybind11). Building is then
managed by [cibuildwheel](https://github.com/joerick/cibuildwheel). Further
installation is then available via `pip install chiapos` e.g.

## Contributing and workflow
Contributions are welcome and more details are available in chia-blockchain's
[CONTRIBUTING.md](https://github.com/Chia-Network/chia-blockchain/blob/main/CONTRIBUTING.md).

The main branch is usually the currently released latest version on PyPI.
Note that at times chiapos will be ahead of the release version that
chia-blockchain requires in it's main/release version in preparation for a
new chia-blockchain release. Please branch or fork main and then create a
pull request to the main branch. Linear merging is enforced on main and
merging requires a completed review. PRs will kick off a GitHub actions ci build
and analysis of chiapos at
[lgtm.com](https://lgtm.com/projects/g/Chia-Network/chiapos/?mode=list). Please
make sure your build is passing and that it does not increase alerts at lgtm.
