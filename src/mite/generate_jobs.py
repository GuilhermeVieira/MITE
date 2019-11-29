#!/usr/bin/python

import math
import os
import sys

path = sys.argv[1] + '/'
start = int(sys.argv[2])
end = int(sys.argv[3])
ncores = int(sys.argv[4])

script = """#!/bin/bash

#SBATCH -o {0}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=9000

python3.6 run_optimization.py {1} {2} ../../output/nexus/ --window_width=2 --window_height=2 -b --f=0.25

exit 0
"""

interval = math.ceil((end - start) / ncores)

if not os.path.exists(os.path.dirname(path)):
    os.makedirs(os.path.dirname(path))

if not os.path.exists(os.path.dirname(path + "out/")):
    os.makedirs(os.path.dirname(path + "/out/"))

for i in range(start, end, interval):
    j = i + interval - 1

    if j > end:
        j = end

    filename = 'job_' + str(i) + '-' + str(j)
    job = open(path + filename + '.job', "w")
    output = path + "out/" + filename + '.out'
    job.write(script.format(output, i, j))
