##################################################################################
##                                                                              ##
##   Module that constructs the ion intensity map (or matrix of intensity by    ##
##   time of elution) from a CSV file containing the m/z and the elution time   ##
##   of each ion detected.                                                      ##
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
from scipy.sparse import coo_matrix, csr_matrix, isspmatrix_csr

# Returns a ion intensity map given a CSV file
def constructIonMap(infile, mz_round=2, etime_round=2):
    df = pd.read_csv(infile,
                     index_col='#',
                     header=0,
                     names=['#', 'mz', 'etime'])
    row = np.array(df.mz.values)
    col = np.array(df.etime.values)
    data = np.ones(row.shape)
    row *= 10 ** mz_round
    col *= 10 ** etime_round
    row = np.trunc(row)
    col = np.trunc(col)
    mite = coo_matrix((data, (row, col)),
                      shape=(int(np.amax(row) + 1), int(np.amax(col) + 1)))
    return mite

# Returns the intersection of two ion intensity maps
def ionMapIntersection(mite1, mite2):
    if (not isspmatrix_csr(mite1)):
        mite1 = mite1.tocsr()
    if (not isspmatrix_csr(mite2)):
        mite2 = mite2.tocsr()
    intersection = mite1.astype(bool).multiply(mite2.astype(bool))
    return intersection 
