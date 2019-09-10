#!/usr/bin/env python3

##################################################################################
##                                                                              ##
##   Module that creates a number of different nexus files based on a           ##
##   combination of different arguments.                                        ##
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

import os
dir_path = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.insert(0, dir_path + '/../src/')

import argparse
import numpy as np
from MiteToNexusWriter import MiteToNexusWriter

def __write_nexus(niter, f, mode, all):
    print('Creating nexus file... (niter=' + str(niter) +
          ', f=' + str(f) + ')')

    if mode == 'both':
        print('    mode=tr')
        mtnw.write_nexus(niter, f, mode='tr', all=all)
        print('    mode=mz')
        mtnw.write_nexus(niter, f, mode='mz', all=all)
    else:
        print('    mode=' + mode)
        mtnw.write_nexus(niter, f, mode=mode, all=all)

parser = argparse.ArgumentParser()
parser.add_argument("niter_start", type=int,
                    help="number of iterations for the first nexus creation")
parser.add_argument("niter_end", type=int,
                    help="number of iterations for the last nexus creation")
parser.add_argument("f_start", type=float,
                    help="f value for the first nexus creation")
parser.add_argument("f_end", type=float,
                    help="f value for the last nexus creation")
parser.add_argument("df", type=float,
                    help="f value step between two consecutive nexus creation")
parser.add_argument("-a", "--all", action="store_true",
                    help="use all LC-MS runs")
parser.add_argument("-m", "--mode", choices=['tr', 'mz', 'both'],
                    default='both',
                    help="define the way that the reduced matrix will be"
                         " further reduced into a string")
args = parser.parse_args()

input_path = dir_path + '/../input/ion_map/xml/'
output_path = dir_path + '/../output/nexus/'
mtnw = MiteToNexusWriter(input_path, output_path)

for niter in range(args.niter_start, args.niter_end + 1):
    if args.df > 0:
        for f in np.arange(args.f_start, args.f_end + args.df, args.df):
            __write_nexus(niter, f, args.mode, args.all)
    else:
        __write_nexus(niter, args.f_start, args.mode, args.all)
