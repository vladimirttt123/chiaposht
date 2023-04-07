# Chia Proof of Space
![Build](https://github.com/Chia-Network/chiapos/workflows/Build/badge.svg)
![PyPI](https://img.shields.io/pypi/v/chiapos?logo=pypi)
![PyPI - Format](https://img.shields.io/pypi/format/chiapos?logo=pypi)
![GitHub](https://img.shields.io/github/license/Chia-Network/chiapos?logo=Github)

[![Total alerts](https://img.shields.io/lgtm/alerts/g/Chia-Network/chiapos.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/Chia-Network/chiapos/alerts/)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/Chia-Network/chiapos.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/Chia-Network/chiapos/context:python)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/Chia-Network/chiapos.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/Chia-Network/chiapos/context:cpp)

With this fork of chiapos I try to add some threads to improve plot creation time.
Important: the main changes are done for enabled bitmap plotting.
For now sortings done with threads here it is example how it improved with thread when no IO limits

```
size: 0.11GiB, threads: 2, times:  read:1.112s, total: 1.998s 
size: 0.141GiB, threads: 3, times:  read:0.494s, total: 1.033s 
size: 0.141GiB, threads: 4, times:  read:0.346s, total: 0.758s 
size: 0.141GiB, threads: 5, times:  read:0.330s, total: 0.690s 
size: 0.141GiB, threads: 6, times:  read:0.261s, total: 0.575s 
size: 0.141GiB, threads: 7, times:  read:0.239s, total: 0.536s 
size: 0.141GiB, threads: 8, times:  read:0.269s, total: 0.512s 
size: 0.141GiB, threads: 9, times:  read:0.240s, total: 0.504s 
size: 0.141GiB, threads: 10, times: read:0.203s, total: 0.463s
```
In case of enough ram next bucket presorted in parallel.

Threads added to phase2 and computational stage 2 of phase 3.

Bitmap of table 7 calculate on write of the table and rewriting of table 7 is done on read of it in phase 3
than no need of table 7 stage in phase 2.

In addition added IO compaction that can lead to faster 
creation with HDD plotting or prolonge SSD life with less writtings.
But compaction do some more work than on slow CPUs or in systems 
where IO is not a botleneck creation will be slower.
To use all implemented compactions buckets number should be at least 256.
The IO compaction allows less ram to plot on ramdisk. 
By my measures to create k32 it needed around 170GiB space in temp directory 
when using compaction with 256 (or more) buckets.

I tested ram plotting on amazon R6i.8xlarge instance with ramdisk of 180Gb,
the creation of k32 plot was taken 52 minutes.


Chia's proof of space is written in C++. Includes a plotter, prover, and
verifier. It exclusively runs on 64 bit architectures. Read the
[Proof of Space document](https://www.chia.net/assets/Chia_Proof_of_Space_Construction_v1.1.pdf) to
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
