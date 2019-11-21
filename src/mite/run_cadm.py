#!/usr/bin/env python3

##################################################################################
##                                                                              ##
##   This module performs the CADM test.                                        ##
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
import sys
sys.path.append("..")

from CADM import CADM

nexus_dirname = '../../output/nexus/'
nexus_directory = os.fsencode(nexus_dirname)
gen = '../../input/genomic/genomic.nex.con.tre'

def __print(output, string):
    print(string)
    output.write(string)

with open('../../reports/report_mite.txt', 'w') as report:
    __print(report,
        "Report comparing topology using the CADM test of the resulting"
        "trees with the tree obtained with the genomic data.\n")
    __print(report, 30*("-") + "\n")

    for d in os.listdir(nexus_directory):
        dname = os.fsdecode(d)

        if (not os.path.isdir(nexus_dirname + dname) or dname.startswith('.')):
            continue

        __print(report, "Target: " + dname + "\n")
        tree = '{0}{1}/mite.nex.con.tre'.format(nexus_dirname, dname)
        __print(report, str(CADM(gen, tree, '-t')) + "\n")
        __print(report, 30*("-") + "\n")
