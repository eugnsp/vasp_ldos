# Depth-resolved density of states calculator for VASP

## Synopsis

Computes the depth- and k-resolved density of states from the WAVECAR file.
The computed data is stored in a binary file and can then be plotted using
a simple MATLAB script.

## Results

(to be written)

## How to run

```none
vasp_ldos [options]

Options:
    -h               print help
    -w <name>        input WAVECAR filename
    -l <name>        output LDOS filename
    -f <value>       Fermi level value
    -c <comment>     arbitrary text comment
```

If no output filename is given, WAVECAR file basic information is displayed
and the program terminates.

## Output file format

| Data type   |  Description                                       |
|:------------|:---------------------------------------------------|
| `char[500]` | Text header tail-padded with spaces                |
| `uint32`    | File format version (= 103)                        |
| `uint32`    | Number of spin projection (1 or 2)                 |
| `uint32`    | Number of k-points                                 |
| `uint32`    | Number of energy bands                             |
| `uint32`    | Number of layers                                   |
| `double`	  | Supercell height	                               |
| `double`    | Fermi level (`-f` option)			               |
| `double`    | Minimum value of `E(k)`                            |
| `double`    | Maximum value of `E(k)`                            |

(to be written)

## External dependencies

* [Intel MKL](https://software.intel.com/en-us/mkl) or [FFTW](http://www.fftw.org/)

## License

This code is distributed under GNU General Public License v3.0.
