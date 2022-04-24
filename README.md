# wallforce
An analysis tool for gromacs trajectories that computes the average force a
harmonic wall applies on a selection of molecules/atoms.

I did not test for the effect of walls close to periodic boundaries. In my 
simulations the walls are not directly at the edge of the box.

There might be bugs!

## Installation

Source Gromacs installation first. Tested with Gromacs 2019.2. There is also an Gromacs 2021 branch!

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

## Usage

Example:
```
wallforce -quiet -f traj_comp.xtc -s topol.tpr -axis z -wallr 6.313665 -wallk 1000 <<< "name OW" -o wallforce.xvg
```

You can get help options with `wallforce -h`.

Remember to specify the following arguments:
```
 -wallr  <real>             (0)
           Wall position as distance from the origin
 -axis   <enum>             (z)
           Axis normal to the wall: x, y, z, xn, yn, zn
 -wallk  <real>             (0)
           Wall force constant
```
