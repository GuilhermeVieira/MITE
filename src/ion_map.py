##################################################################################
##                                                                              ##
##   Module that manipulates the ion intensity maps (or matrix of intensity by  ##
##   time of elution) extracted from mass spectrometry raw files.               ##
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

import decimal
import numpy as np
import os
import xml.etree.ElementTree as ET
from scipy.sparse import coo_matrix, csr_matrix

# Resize the matrices if their shapes are different.
def __resizeMatrices(mite1, mite2):
    if (mite1.shape != mite2.shape):
        s0 = max(mite1.shape[0], mite2.shape[0])
        s1 = max(mite1.shape[1], mite2.shape[1])
        mite1.resize((s0, s1))
        mite2.resize((s0, s1))

# Returns an ion intensity map in COO format given a XML file
def constructIonMap(filepath):
    basename, ext = os.path.splitext(filepath)

    if (ext == '.xml'):
        tree = ET.parse(filepath)
        root = tree.getroot()
        features_count = int(root.find('LC_MS_RUN').get('number_of_features'))
        row = np.empty(features_count)
        col = np.empty(features_count)
        data = np.ones(features_count, dtype=bool)
        tr_round = np.empty(features_count)
        mz_round = np.empty(features_count)
        index = 0

        for feature in root.iter('MS1_FEATURE'):
            tr = feature.get('Tr')
            mz = feature.get('m_z')
            row[index] = tr
            col[index] = mz
            tr_round[index] = abs(decimal.Decimal(tr).as_tuple().exponent)
            mz_round[index] = abs(decimal.Decimal(mz).as_tuple().exponent)
            index += 1

        row *= 10 ** np.amax(mz_round)
        col *= 10 ** np.amax(tr_round)
        row = np.trunc(row)
        col = np.trunc(col)
        shape = (int(np.amax(row) + 1), int(np.amax(col) + 1))
        mite = coo_matrix((data, (row, col)), shape)

        return mite

    return None

# Returns the intersection of two ion intensity maps
def ionMapIntersection(mite1, mite2):
    m1, m2 = mite1.copy(), mite2.copy()
    __resizeMatrices(m1, m2)
    intersec = m1.multiply(m2)
    return intersec

# Returns the symmetric difference of two ion intensity maps
def ionMapSymmetricDiff(mite1, mite2):
    m1, m2 = mite1.copy(), mite2.copy()
    __resizeMatrices(m1, m2)
    sym_diff = (m1 - m2) + (m2 - m1)
    return sym_diff

# Returns a distance matrix given a list of mites
def calculateDistMatrix(mites):
    mites_inter = []
    for m in mites:
        inter = ionMapIntersection(m[0], m[1])
        mites_inter.append(inter)
    dist_matrix = np.zeros([len(mites_inter), len(mites_inter)])
    np.fill_diagonal(dist_matrix, 0)
    for i in range(len(mites_inter)):
        for j in range(i + 1, len(mites_inter)):
            sym_diff = ionMapSymmetricDiff(mites_inter[i], mites_inter[j])
            dist_matrix[i][j] = sym_diff.count_nonzero()
    return dist_matrix
