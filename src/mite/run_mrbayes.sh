#!/bin/sh

for dir in ../../output/nexus/*; do
    mb "$dir/mite.nex"
done
