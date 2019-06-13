#!/usr/bin/env python3

##################################################################################
##                                                                              ##
##   Script that plots the ion maps given XML files generated from SuperHirn.   ##
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

import matplotlib.pylab      as plt
import numpy                 as np
import os
import pandas                as pd
import sys
import xml.etree.ElementTree as ET
from scipy.sparse import coo_matrix

mz_round = 2
tr_round = 2

if(len(sys.argv) != 2):
    print('Usage:', os.path.basename(sys.executable), sys.argv[0],
            'files_path')
    print('')
    print('\t%-16s %s' % ('ARGUMENT', 'DESCRIPTION'))
    print('\t%-16s %s' % ('files_path', 'the path to XML files'))
    sys.exit(0)

files_path = sys.argv[1]
if not files_path.endswith('/'):
    files_path += '/'

for filename in sorted(os.listdir(files_path)):
    basename, ext = os.path.splitext(filename)
    if (ext == '.xml'):
        tree = ET.parse(files_path + filename)
        root = tree.getroot()
#        features_count = int(root.find('LC_MS_RUN').get('number_of_features'))
        features_count = 0
        for feature in root.iter('MS1_FEATURE'):
            features_count += 1
        row = np.empty(features_count)
        col = np.empty(features_count)
        data = np.ones(features_count, dtype=bool)
        index = 0
        for feature in root.iter('MS1_FEATURE'):
            row[index] = feature.get('m_z')
            col[index] = feature.get('Tr')
            index += 1
        row *= 10 ** mz_round
        col *= 10 ** tr_round
        row = np.trunc(row)
        col = np.trunc(col)
        mite = coo_matrix((data, (row, col)),
                shape=(int(np.amax(row) + 1), int(np.amax(col) + 1)))
        plt.spy(mite, markersize=0.1, aspect='auto', color='black')
        plt.title(filename)
        plt.xlabel('m/z')
        plt.ylabel('retention time')
        plt.show()
