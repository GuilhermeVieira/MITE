#!/usr/bin/env python3

##################################################################################
##                                                                              ##
##   This module performs the CADM test.                                        ##
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
import os
import sys
sys.path.append("..")

from CADM import CADM

def __print(output, string):
    #print(string)
    output.write(string)

def run(nexus_path='../../output/nexus/bob/mite.nex.con.tre', reports_path='../../reports/bob/bob_reports.txt'):
    gen = '../../input/genomic/genomic.nex.con.tre'

    with open(reports_path, 'w') as report:
        __print(report,
            "\nReport comparing topology using the CADM test of the resulting"
            "trees with the tree obtained with the genomic data.\n")
        __print(report, 30*("-") + "\n")
        __print(report, "Target: MITE\n")
        tree = '{0}'.format(nexus_path)
        res, stats = CADM(gen, tree, '-t')
        __print(report, str(res) + "\n")
        __print(report, 30*("-") + "\n")

    return stats

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "nexus_filepath", type=str, metavar="NEXUS_FILEPATH", help="the path "
        "of the tree NEXUS file which will be tested against the genomic tree"
    )
    parser.add_argument(
        "output", type=str, metavar="OUTPUT_PATH", help="the output file path"
    )
    args = parser.parse_args()

    run(args.nexus_filepath, args.output)
