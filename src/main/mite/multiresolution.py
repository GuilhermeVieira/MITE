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

import logging
import math
from typing import List

import numpy as np
from scipy.sparse import lil_matrix

from Mite import Mite


class Partitioner:

    def __init__(self, w, h):
        self.__w = w
        self.__h = h
        self.__max_token_length = 99990

    def partition_mites(self, mites: List[Mite]):
        parts_list = []

        for mite in mites:
            parts_list.append(self.__get_parts(mite.matrix))

        flattened_mites = self.__flatten_parts(parts_list, self.__w, self.__h)

        return flattened_mites

    # Flattens the matrices parts
    def __flatten_parts(self, parts_list, w, h):
        flattened_mites = []
        for parts in parts_list:
            flattened_mite = np.empty(0)
            for p in parts:
                r = self.reduce_dim(p, w, h, max_size=(self.__max_token_length / len(parts)))
                flattened_mite = np.concatenate((flattened_mite, self.__matrix2array(r)))
            flattened_mites.append(flattened_mite)

        flattened_mites = np.array(flattened_mites, dtype=int)
        #flattened_mites = self.__remove_equal_columns(flattened_mites);

        return flattened_mites

    # Returns the window (and its position) to which a matrix element belongs
    # and the position of this element in the window
    def __get_window(self, matrix, pos, w, h):
        new_i, win_i = divmod(pos[0], h)
        new_j, win_j = divmod(pos[1], w)
        i = pos[0] - pos[0] % h
        j = pos[1] - pos[1] % w

        return matrix[i:i+h, j:j+w], (new_i, new_j), (win_i, win_j)

    # Returns a new matrix with reduced dimensionality
    def __reduce_dim(self, old, w, h):
        s0 = math.ceil(old.shape[0] / w)
        s1 = math.ceil(old.shape[1] / h)
        new = lil_matrix((s0, s1))
        iarray, jarray = old.nonzero()

        if len(iarray) != 0:
            for i, j in np.nditer([iarray, jarray]):
                win, new_pos, win_pos = self.__get_window(old, (i, j), w, h)
                # Previous partition algorithm discretization
                '''
                if not binary:
                    if win.nnz != 0:
                        new[new_pos] = int(round(np.sum(win) / win.nnz))
                    else:
                        new[new_pos] = 0
                else:
                    if win.nnz / (win.shape[0] * win.shape[0]) >= f:
                        new[new_pos] = True
                '''
                # SHOULD IT BE DIVIDED THE TOTAL NUMBER OF ELEMENTS? OR NO NNZ ELEMENTS?
                if win.nnz != 0:
                    new[new_pos] = int(round(np.sum(win) / win.nnz))
                else:
                    new[new_pos] = 0.0

        return new.tocsr(copy=True)

    # Returns a new matrix with reduced dimensionality (wrapper)
    def reduce_dim(self, matrix, w, h, max_size=math.inf):
        if w >= matrix.shape[0]:
            w = matrix.shape[0]
        if h >= matrix.shape[1]:
            h = matrix.shape[1]

        reduced = matrix.copy()
        x = 0

        while reduced.shape[0] * reduced.shape[1] > max_size:
            reduced = self.__reduce_dim(reduced, w, h)
            x += 1

        return reduced

    # Partitioning code that was not being used
    '''
    # Returns a partition of a matrix, such that the number of parts in row and
    # column are equal
    def generate_uniform_partition(self, shape, pcount_row, pcount_col):
        logging.info(
            'Generating an partition dividing rows and cols in ' + str(pcount_row)
            + 'x' + str(pcount_col) + ' parts'
        )
        partition = []
        step_row = int(math.ceil(shape[0] / pcount_row))
        step_col = int(math.ceil(shape[1] / pcount_col))

        for i in range(0, shape[0], step_row):
            for j in range(0, shape[1], step_col):
                slicex = slice(i, i + step_row)
                slicey = slice(j, j + step_col)
                partition.append((slicex, slicey))

        return partition

    # Returns a partition with square types
    def generate_square_partition(self, shape, part_size_row, part_size_col):
        logging.info(
            'Generating a partition with parts of size ' + str(part_size_row) +
            'x' + str(part_size_col)
        )
        partition = []

        for i in range(0, shape[0], part_size_row):
            for j in range(0, shape[1], part_size_col):
                slicex = slice(i, i + part_size_row)
                slicey = slice(j, j + part_size_col)
                partition.append((slicex, slicey))

        return partition
    
    '''

    # Returns the parts of the partition applied to the matrix
    def __get_parts(self, matrix):

        m = matrix.shape[0]
        n = matrix.shape[1]
        partition = [(slice(0, m), slice(0, n))]

        parts = []

        #print(partition)

        for p in partition:
            parts.append(matrix[p[0], p[1]])

        return parts

# Reduces a 2d matrix to a string
    def __matrix2array(self, matrix):
        array = np.empty(0)

        for i in range(0, matrix.shape[0]):
            a = matrix.getrow(i).toarray().astype(int).flatten()
            array = np.concatenate((array, a))

        return array

    # Deletes all matrix columns that all values are the same
    def __remove_equal_columns(self, matrix):
        index = np.argwhere(np.all(matrix == matrix[0, :], axis=0))
        matrix = np.delete(matrix, index, axis=1)
        return matrix

