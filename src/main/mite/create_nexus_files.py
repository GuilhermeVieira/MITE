#!/usr/bin/env python3

##################################################################################
##                                                                              ##
##   Module that creates a NEXUS file based on arguments passed on command      ##
##   line.                                                                      ##
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
import logging
import numpy as np
import os

from mite.MiteToNexusWriter import MiteToNexusWriter

def run(mtnw, w, h, binary, f, partition, name="tres-especies"):
    logging.info(
        'Creating NEXUS file... (w=' + str(w) + ', h=' + str(h) +
        ', binary=' + str(binary) + ', f=' + str(f) + ')'
    )

    return mtnw.write_nexus(w, h, partition, name, f=f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("window_width", type=int,
                        help="the width of the window used to reduce the MITEs")
    parser.add_argument("window_height", type=int,
                        help="the height of the window used to reduce the MITEs")
    parser.add_argument("--binary", action="store_true",
                        help="use binary MITEs to create the NEXUS files")
    parser.add_argument("--f", type=float,
                        help="value for the frequency of 1's necessary to get 1 "
                             "in the new window")

    args = parser.parse_args()

    if (args.binary and args.f is None):
        parser.error("The --binary argument requires --f to be set")

    if (not args.binary and args.f is not None):
        parser.error("The argument --f requires the --binary flag to be set")

    dir_path = os.path.dirname(os.path.realpath(__file__))
    input_path = dir_path + '/../../input/ion_map/xml/'
    output_path = dir_path + '/../../output/nexus/'
    mtnw = MiteToNexusWriter(input_path, output_path, args.binary)
    m = mtnw.mites[0].matrix.shape[0]
    n = mtnw.mites[0].matrix.shape[1]
    partition = [(slice(0, m), slice(0, n))]
    run(
        mtnw, args.window_width, args.window_height, args.binary, args.f,
        partition
    )
