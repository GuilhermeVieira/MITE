#!/bin/bash

for dir in ../output/nexus/*/; do
    mb "$dir/mite.nex"
done
