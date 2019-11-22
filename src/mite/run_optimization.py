#!/usr/bin/env python3

##################################################################################
##                                                                              ##
##   This module implements optimization routines based on the                  ##
##   multiresolution method.                                                    ##
##                                                                              ##
##   This file is part of the featsel program                                   ##
##   Copyright (C) 2019 Gustavo Mendes Maciel                                   ##
##                                                                              ##
##   This program is free software: you can redistribute it and/or modify       ##
##   it under the terms of the GNU General Public License as published by       ##
##   the Free Software Foundation, either version 3 of the License, or          ##
##   (at your option) any later version.                                        ##
##                                                                              ##
##   This program is distributed in the hope that it will be useful,            ##
##   but WITHOUT ANY WARRANTY; without even the implied warranty of             ##
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              ##
##   GNU General Public License for more details.                               ##
##                                                                              ##
##   You should have received a copy of the GNU General Public License          ##
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.      ##
##                                                                              ##
##################################################################################

import argparse
import math
import os
import subprocess

from MiteToNexusWriter import MiteToNexusWriter
import create_nexus_files as nexus
#import run_cadm as cadm

dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
log_path = dir_path + '../../log/'

max_pcount_row = 20 
max_pcount_col = 20

# Returns a uniform partition that will be used in all ion maps.
def generate_uniform_partition(shape, pcount_row, pcount_col):
    partition = []
    step_row = int(math.ceil(shape[0] / pcount_row))
    step_col = int(math.ceil(shape[1] / pcount_col))

    for i in range(0, shape[0], step_row):
        for j in range(0, shape[1], step_col):
            slicex = slice(i, i + step_row)
            slicey = slice(j, j + step_col)
            partition.append((slicex, slicey))

    return partition

# Constructs the NEXUS file from the partitioned ion maps.
def construct_nexus(mtnw, args, partition):
    nexus.run(
        mtnw, args.window_width, args.window_height, args.binary,
        args.f, partition
    )

    return

# Runs MrBayes with the NEXUS file as input.
def run_mrbayes(nexus_path, run_num):
    subprocess.call(
        ['mb', dir_path + nexus_path + '/mite.nex', '>',
         log_path + 'mb/mb_log' + run_num + '.txt']
    )
    return

# Runs the CADM test, comparing the constructed tree with the mtDNA tree.
def run_cadm():
    return

# Assesses the results from the CADM test.
def assess_results():
    return

# Run optimization procedure
def run_basic_optimization(args):
    mtnw = MiteToNexusWriter(args.xml_path, args.nexus_path, args.binary)
    shape = mtnw.mites[0].matrix.shape

    for i, j in zip(range(1, max_pcount_row + 1), range(1, max_pcount_col + 1)):
        partition = generate_uniform_partition(shape, i, j)
        construct_nexus(mtnw, args, partition)
        run_mrbayes(args.nexus_path, i)
        run_cadm()
        assess_results()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "nexus_path", type=str, help="the path to store generated nexus "
        "files"
    )
    parser.add_argument(
        "--xml_path", type=str, help="the path to the XML files that will "
        "be used to construct the ion maps",
        default=dir_path+"/../../input/ion_map/xml/"
    )
    parser.add_argument(
        "--window_width", type=int, help="the width of the window used to "
        "reduce the MITEs", required=True
    )
    parser.add_argument(
        "--window_height", type=int, help="the height of the window used to "
        "reduce the MITEs", required=True
    )
    parser.add_argument(
        "--binary", action="store_true", help="use binary MITEs to create the "
        "NEXUS files"
    )
    parser.add_argument(
        "--f", type=float, help="value for the frequency of 1's necessary to "
        "get 1 in the new window"
    )

    args = parser.parse_args()

    if (args.binary and args.f is None):
        parser.error("The --binary argument requires --f to be set")

    if (not args.binary and args.f is not None):
        parser.error("The argument --f requires the --binary flag to be set")

    run_basic_optimization(args)
