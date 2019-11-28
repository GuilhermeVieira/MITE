#!/bin/sh

if [ -z "$1" ]
then
  echo "Specify a path"
else
  for file in $1/*; do
    sbatch $file
    sleep 1
  done
fi
