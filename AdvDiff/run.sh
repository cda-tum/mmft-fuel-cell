#!/bin/bash
nproc=$1
N=$2

cW = 30
cL = 30
h = 40
l = 30
w = 30

if [ ! -d "./tmpCube_cW${cW}cL${cL}H${h}L${l}W${w}/" ]; then
    mkdir ./tmpCube_cW${cW}cL${cL}H${h}L${l}W${w}/
fi
mpirun -np ${nproc} ./customDiffusion3D ${N} channels3D/Cube_cW${cW}_cL${cL}_H${h}_L${l}_W${w}.stl ./tmpCube_cW${cW}cL${cL}H${h}L${l}W${w}/ &> Cube_cW${cW}cL${cL}H${h}L${l}W${w}.out
