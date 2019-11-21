#!/usr/bin/env python3

##################################################################################
##                                                                              ##
##   This module implements the following pipeline:                             ##
##   1. Generate NEXUS files from ion maps                                      ##
##   2. Run MrBayes for each NEXUS file created                                 ##
##   3. Run CADM test for each tree generated from MrBayes                      ##
##                                                                              ##
##   This file is part of the featsel program                                   ##
##   Copyright (C) 2018 Victor Wichmann Raposo                                  ##
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
import numpy as np
import os
import subprocess

from MiteToNexusWriter import MiteToNexusWriter
import create_nexus_files as nexus
import run_cadm as cadm

def __print(output, string):
    print(string)
    output.write(string)

def create_nexus_files(mtnw, w, h, binary, f):
    print("Creating the NEXUS files...")
    nexus.run(mtnw, w, h, binary, f)

def run_mb(nproc):
    process = subprocess.call(["./run_mrbayes.sh", str(nproc)])

def run_cadm():
    cadm.run()

def run(mtnw, w, h, binary, f):
    create_nexus_files(
        mtnw, args.window_width, args.window_height, args.binary, args.f
    )
    run_mb(args.nproc)
    run_cadm()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--nproc", type=int, help="number of processors that will be used in "
        "MrBayes analysis", required=True
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

    dir_path = os.path.dirname(os.path.realpath(__file__))
    input_path = dir_path + '/../../input/ion_map/xml/'
    output_path = dir_path + '/test/'
    mtnw = MiteToNexusWriter(input_path, output_path, args.binary)
    run(mtnw, args.window_width, args.window_height, args.binary, args.f)
