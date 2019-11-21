#!/usr/bin/env python3

##################################################################################
##                                                                              ##
##   This module implements functions needed by the multiresolution pipeline.   ##
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

import math
import numpy as np
from scipy.sparse import csr_matrix, lil_matrix

# Returns the window (and its position) to which a matrix element belongs
# and the position of this element in the window
def __get_window(matrix, pos, w, h):
    new_i, win_i = divmod(pos[0], h)
    new_j, win_j = divmod(pos[1], w)
    i = pos[0] - pos[0] % h
    j = pos[1] - pos[1] % w

    return matrix[i:i+h, j:j+w], (new_i, new_j), (win_i, win_j)

# Returns a new matrix with reduced dimensionality
def __reduce_dim(old, f, w, h, binary):
    s0 = math.ceil(old.shape[0] / w)
    s1 = math.ceil(old.shape[1] / h)
    new = lil_matrix((s0, s1))
    iarray, jarray = old.nonzero()

    if len(iarray) != 0:
        for i, j in np.nditer([iarray, jarray]):
            win, new_pos, win_pos = __get_window(old, (i, j), w, h)
            if not binary:
                if win.nnz != 0:
                    new[new_pos] = int(round(np.sum(win) / win.nnz))
                else:
                    new[new_pos] = 0
            else:
                if win.nnz / (win.shape[0] * win.shape[0]) >= f:
                    new[new_pos] = True

    return new.tocsr(copy=True)

# Returns a new matrix with reduced dimensionality
def reduce_dim(matrix, w, h, binary, max_size=math.inf, f=0.0):
    if w >= matrix.shape[0]:
        raise ValueError('w is too big!')
    if h >= matrix.shape[1]:
        raise ValueError('h is too big!')

    reduced = matrix.tocsr(copy=True)
    x = 0

    while reduced.shape[0] * reduced.shape[1] > max_size:
        reduced = __reduce_dim(reduced, f, w, h, binary)
        x += 1

    return reduced 
