# MMFT Fuel Cell

* What is this repository
* It is developed by Chair for Design Automation at the Technical University of Munich
* Munich MicroFluidic Toolkit (MMFT)
* It allows to simulate the diffusion and advection of hydrogen in a unidirectional flow over an anode with nanostructures.
* Additionally, the source code for the geometry generation is included, so
* Citations ?
* We used palabos, which is an open source LBM solver, which can be installed here. [cite]
* For more information about our work on Microfluidics, please visit cda.cit.tum.de /research/Microfluidics
* If you have any questions, feel free to contact us via microfluidics... or creating an issue on GitHub

## System Requirements
* Same as for palabos
* Python 3.10 to generate the geometries
* tested only for ubuntu 20.6

## Usage


#### Palabos or the fluid solver
Install Palabos from here, for more info see.

In file XXX specify the location of your palabos library

then do

`cmake`

`make`

for a test run please do

`./runall.sh 1 50`

where 1 is nproc, and 50 is N, the amount of cells in the z direction

#### Geometry generation

To generate a shape with dimensions, please see the help output for each shape, e.g.

`python3 geometry -cube`

The available shapes are

* cube
* Step
* Groove
* Offset (for the offset cubes)

## Example


## Reference
