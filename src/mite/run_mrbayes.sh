#!/bin/sh

for dir in ../../output/nexus/*; do
    mpirun -np $1 mb "$dir/mite.nex"
done
