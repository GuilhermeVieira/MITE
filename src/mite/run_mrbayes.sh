#!/bin/sh

if [ -z "$1" ]
then
    echo "You have to specify the number of processors that will be used to run MrBayes!"
else
    for dir in ../../output/nexus/*; do
        printf '\n'
        echo $dir
        mpirun -np $1 mb "$dir/mite.nex"
    done
fi
