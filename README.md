# wallforce
An analysis tool for gromacs trajectories that computes the average force a
harmonic wall applies on a selection of molecules/atoms.

## Installation

Source Gromacs installation first. Tested with Gromacs 2019.2

You will need to make sure that you use the same C++ compiler
and C++ Standard Library as the one that was used for compiling
GROMACS.

### cmake
```bash
mkdir build
cd build
cmake ..
make
```

### Makefile
```bash
make -f Makefile.pkg
```
