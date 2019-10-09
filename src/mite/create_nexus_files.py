#!/usr/bin/env python3

##################################################################################
##                                                                              ##
##   Module that creates a number of different NEXUS files based on a           ##
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


import argparse
import numpy as np
import os

from MiteToNexusWriter import MiteToNexusWriter


def __write_nexus(w, h, min_niter, append_by_row, binary, f):
    print('Creating NEXUS file... (w=' + str(w) + ', h=' + str(h) +
          ', min_niter=' + str(min_niter) +
          ', append_by_row=' + str(append_by_row) +
          ', binary=' + str(binary) + ', f=' + str(f), end = '')
    print(')')

    mtnw.write_nexus(w, h, min_niter, append_by_row=append_by_row, binary=binary, f=f)

parser = argparse.ArgumentParser()

parser.add_argument("window_width", type=int,
                    help="the width of the window used to reduce the MITEs")
parser.add_argument("window_height", type=int,
                    help="the height of the window used to reduce the MITEs")
parser.add_argument("--min_niter_start", type=int, default=1,
                    help="start value for the minimum number of iterations")
parser.add_argument("--min_niter_end", type=int, default=1,
                    help="end value for the minimum number of iterations")
parser.add_argument("--append_by_row", action="store_true",
                    help="append all the MITE's rows to form the final string")
parser.add_argument("--binary", action="store_true",
                    help="use binary MITEs to create the NEXUS files")
parser.add_argument("--f_start", type=float,
                    help="start value for the frequency of 1's necessary to "
                         "get 1 in the new window")
parser.add_argument("--f_end", type=float,
                    help="start value for the frequency of 1's necessary to "
                         "get 1 in the new window")
parser.add_argument("--df", type=float,
                    help="step between the f value of two consecutive NEXUS "
                         "file creation")

args = parser.parse_args()

if (args.binary and
    (args.f_start is None or args.f_end is None or args.df is None)):
    parser.error("The --binary argument requires --f_start, --f_end and --df "
                 "arguments")

if (not args.binary and
    (args.f_start is not None or args.f_end is not None or args.df is not None)):
    parser.error("The arguments --f_start, --f_end and --df requires the "
                 "--binary argument")

dir_path = os.path.dirname(os.path.realpath(__file__))
input_path = dir_path + '/../../input/ion_map/xml/aligned/'
output_path = dir_path + '/../../output/nexus/'
mtnw = MiteToNexusWriter(input_path, output_path)

for min_niter in range(args.min_niter_start, args.min_niter_end + 1):
    if args.binary and args.df > 0:
        for f in np.arange(args.f_start, args.f_end + args.df, args.df):
            __write_nexus(args.window_width, args.window_height, min_niter,
                          args.append_by_row, args.binary, f)
    else:
        __write_nexus(args.window_width, args.window_height, min_niter,
                      args.append_by_row, args.binary, args.f_start)
