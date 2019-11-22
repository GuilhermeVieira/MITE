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

import ntpath
import os
import sys
sys.path.append("..")

from CADM import CADM

def __print(output, string):
    #print(string)
    output.write(string)

def run(nexus_path, reports_path):
    gen = '../../input/genomic/genomic.nex.con.tre'
    basename = ntpath.basename(nexus_path)

    with open(reports_path, 'w') as report:
        __print(report,
            "\nReport comparing topology using the CADM test of the resulting"
            "trees with the tree obtained with the genomic data.\n")
        __print(report, 30*("-") + "\n")
        __print(report, "Target: MITE\n")
        tree = '{0}/mite.nex.con.tre'.format(nexus_path)
        __print(report, str(CADM(gen, tree, '-t')) + "\n")
        __print(report, 30*("-") + "\n")

if __name__ == '__main__':
    run()
