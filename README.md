# iPRIns

![Ubuntu](https://github.com/cucpbioinfo/iprins/workflows/Ubuntu/badge.svg?branch=master)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/cucpbioinfo/iprins)](https://github.com/cucpbioinfo/iprins/releases)

Improve Precision and Recall Insertions


## Installation

### 1. Installing iPRIns from source

Requirements
1. GCC >=v5.4
2. Intel tbb-devel >=2019.9
3. HTSlib >=v1.9
4. CMake >=v3.9

```sh
git clone https://github.com/cucpbioinfo/iPRIns.git
mkdir build
cd build
cmake .. -DINSTALL_BIN_PREFIX=(directory path for the binary executable file of iPRIns) -DINCLUDE_LIBRARY_PREFIX=(directory path for include of libraries) -DLIBRARY_LINK_PREFIX=(directory path for lib of libraries)
make
make install
```

## Usage
```sh
iprins call -b (aligment) -r (reference) -o (output)
```


## License
iPRIns is distributed under the [GPLv3 License][License]. However, iPRIns includes several third-party open-source libraries, please see [Third Party Software Notices][LICENSETHIRDPARTY] for details.


[License]:LICENSE
[LICENSETHIRDPARTY]:THIRD-PARTY-LICENSE