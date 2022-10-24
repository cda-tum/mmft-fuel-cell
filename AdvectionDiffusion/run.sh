#!/bin/bash
nproc=$1
N=$2

if [ ! -d "./tmpCube/" ]; then
    mkdir ./tmpCube/
fi
mpirun -np ${nproc} ./AdvectionDiffusion ${N} ChannelGeometries/Cube.stl ./tmpCube/ &> Cube.out
