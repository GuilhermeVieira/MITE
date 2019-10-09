##################################################################################
##                                                                              ##
##   Class that implements the ion intensity maps (or matrix of intensity by    ##
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
import math
import numpy as np
import os
from scipy.sparse import coo_matrix, csr_matrix, lil_matrix

from MiteXMLReader import MiteXMLReader


class Mite:

    # Class constructor
    def __init__(self, filepath, binary=False):
        xml_reader = MiteXMLReader(filepath)

        self.binary = binary
        self.tr_min = xml_reader.get_tr_min()
        self.tr_max = xml_reader.get_tr_max()
        self.mz_min = xml_reader.get_mz_min()
        self.mz_max = xml_reader.get_mz_max()

        data, row, col, tr_round, mz_round = xml_reader.construct_ionmap(binary=binary)
        shape = (int((self.tr_max - self.tr_min) * 10 ** tr_round),
                int((self.mz_max - self.mz_min) * 10 ** mz_round))
        self.matrix = coo_matrix((data, (row, col)), shape)

    # Returns the window (and its position) to which a matrix element belongs
    # and the position of this element in the window
    def __get_window(self, matrix, pos, w, h):
        new_i, win_i = divmod(pos[0], h)
        new_j, win_j = divmod(pos[1], w)
        i = pos[0] - pos[0] % h
        j = pos[1] - pos[1] % w
        return matrix[i:i+h, j:j+w], (new_i, new_j), (win_i, win_j)

    # Returns a new matrix with reduced dimensionality
    def __reduce_dim(self, old, f, w, h):
        s0 = math.ceil(old.shape[0] / w)
        s1 = math.ceil(old.shape[1] / h)
        new = lil_matrix((s0, s1))
        iarray, jarray = old.nonzero()

        if len(iarray) != 0:
            for i, j in np.nditer([iarray, jarray]):
                win, new_pos, win_pos = self.__get_window(old, (i, j), w, h)
                if not self.binary:
                    if win.nnz != 0:
                        new[new_pos] = int(round(np.sum(win) / win.nnz))
                    else:
                        new[new_pos] = 0
                else:
                    if win.nnz / (win.shape[0] * win.shape[0]) >= f:
                        new[new_pos] = True

        return new.tocsr(copy=True)


    # Returns a new matrix with reduced dimensionality
    def reduce_dim(self, w, h, min_niter=1, max_size=math.inf, f=0.0):
        if w >= self.matrix.shape[0]:
            raise ValueError('w is too big!')
        if h >= self.matrix.shape[1]:
            raise ValueError('h is too big!')

        reduced = self.matrix.tocsr(copy=True)
        x = 0

        while x < min_niter or reduced.shape[0] * reduced.shape[1] > max_size:
            reduced = self.__reduce_dim(reduced, f, w, h)
            x += 1

        return reduced 

    # Returns the intersection with another mite
    def intersect(self, mite):
        m1, m2 = self.matrix.tocsr(copy=True), mite.matrix.tocsr(copy=True)
        intersection = m1.multiply(m2)
        return intersection

    # Returns the symmetric difference with another mite
    def calculate_symmdiff(self, mite):
        m1, m2 = self.matrix.tocsr(copy=True), mite.matrix.tocsr(copy=True)
        symmdiff = (m1 - m2) + (m2 - m1)
        return symmdiff
