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

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

# Resize and returns the matrices if their shapes are different.
def __resizeMatrices(mite1, mite2):
    if (mite1.shape != mite2.shape):
        s0 = max(mite1.shape[0], mite2.shape[0])
        s1 = max(mite1.shape[1], mite2.shape[1])
        mite1.resize((s0, s1))
        mite2.resize((s0, s1))
    return mite1, mite2

# Returns an ion intensity map in CSR format given a CSV file
def constructIonMap(filepath):
    df = pd.read_csv(filepath,
                     index_col=0,
                     na_filter=False,
                     compression='bz2')
    mite = csr_matrix(df.to_numpy(dtype=bool))
    return mite

# Returns the intersection of two ion intensity maps
def ionMapIntersection(mite1, mite2):
    mite1, mite2 = __resizeMatrices(mite1, mite2)
    intersec = mite1.multiply(mite2)
    return intersec

# Returns the symmetric difference of two ion intensity maps
def ionMapSymmetricDiff(mite1, mite2):
    mite1, mite2 = __resizeMatrices(mite1, mite2)
    sym_diff = (mite1 - mite2) + (mite2 - mite1)
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
