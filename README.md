# Depth-resolved density of states calculator for VASP

## Synopsis

Computes the depth- and k-resolved density of states from the WAVECAR file.
The computed data is stored in a binary file and can then be plotted using
a simple MATLAB script.

## How to build

If FFTW is found, it will be used. Otherwise, `MKLROOT` environment variable
should be set to point to the MKL installation directory. Then:

```sh
git clone --recursive https://github.com/eugnsp/vasp_ldos.git
cd vasp_ldos
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE .. && make
```

to build the tool.

C++17 compiler is required.

## How to run

```none
vasp_ldos [options]

Options:
    -h               print help
    -o <name>        output LDOS filename (no default)
    -w <name>        input WAVECAR filename (default: "WAVECAR")
    -f <value>       Fermi level value (default: 0)
    -c <comment>     arbitrary text comment (default: none)
```

If no output filename is given, `WAVECAR` file basic information is displayed
and the program terminates.

## Output file format

Header:

| Data type    | Size      |  Description                                           |
|:-------------|:---------:|:-------------------------------------------------------|
| `char[500]`  | `500`     | Text header (tail-padded with spaces)                  |
| `uint32`     | `4`       | File format version (= `103`)                          |
| `double[3]`  | `24`      | Real space basis <code>a<sub>i</sub></code>            |
| `double[3]`  | `24`      | Reciprocal space basis <code>b<sub>i</sub></code>      |
| `uint32`     | `4`       | Number of spin projection (`1` or `2`)                 |
| `uint32`     | `4`       | Number of k-points (`nkpt`)                            |
| `uint32`     | `4`       | Number of energy bands (`nb`)                          |
| `uint32`     | `4`       | Number of layers (`nl`)                                |
| `double`	   | `8`       | Supercell height	                                    |
| `double`     | `8`       | Fermi level (specified by the `-f` option)	            |
| `double`     | `8`       | Minimum value of `E(k)`                                |
| `double`     | `8`       | Maximum value of `E(k)`                                |

Then `nkpt` blocks follow:

| Data type        | Size          |  Description                                        |
|:-----------------|:-------------:|:----------------------------------------------------|
| `double[3]`      | `24`          | `k` point                                           |
| `double[nb]`     | `8 * nb`      | Energies <code>E<sub>n</sub></code>                 |
| `double[nb]`     | `8 * nb`      | Occupations <code>nocc<sub>n</sub></code>           |
| `float[nb * nl]` | `4 * nb * nl` | LDOS <code>&rho;<sub>n</sub>(z<sub>l</sub>)</code>  |

## External dependencies

* [Intel MKL](https://software.intel.com/en-us/mkl) or [FFTW](http://www.fftw.org/)

## License

This code is distributed under GNU General Public License v3.0.
