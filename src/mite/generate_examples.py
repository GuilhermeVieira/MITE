#!/usr/bin/env python3

##################################################################################
##                                                                              ##
##   Script that generates examples of ion map CSV files given an original      ##
##   file, the number of examples to be generated and the probability of        ##
##   changing a value from 1 to 0.                                              ##
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

import os
import pandas as pd
import random
import sys

if (len(sys.argv) != 4):
    print('Usage:', os.path.basename(sys.executable), sys.argv[0],
            'original_file file_count probability')
    print('')
    print('\t%-16s %s' % ('ARGUMENT', 'DESCRIPTION'))
    print('\t%-16s %s' % ('orginal_file', 'original CSV file'))
    print('\t%-16s %s' % ('file_count', 'number of files to be generated'))
    print('\t%-16s %s' % ('probability', 'probability of changing a value from'
        ' 1 to 0'))
    sys.exit(0)

filename, file_extension = os.path.splitext(sys.argv[1])

try:
    file_count = int(sys.argv[2])
    p = float(sys.argv[3])
except ValueError:
    print('Error: file_count must be an integer and p must be a float!')
    sys.exit(0)

print('Original CSV file:', filename + file_extension)
print('Number of files to be generated:', file_count)
print('Probability of change:', p)

original_mite = pd.read_csv(filename + file_extension,
                            index_col='#',
                            header=0,
                            names=['#', 'm/z', 'Retention Time (min)'])

print('\nFiles generated:')

for i in range(file_count):
    new_mite = original_mite.copy()
    for index, row in new_mite.iterrows():
        if random.random() < p:
            new_mite.drop(index, inplace=True)
    full_name = filename + '_' + str(i) + file_extension
    new_mite.to_csv(full_name)
    print(full_name)
