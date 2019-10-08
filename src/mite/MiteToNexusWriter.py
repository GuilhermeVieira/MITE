##################################################################################
##                                                                              ##
##   Module that writes a nexus file from dimensionally reduced mites.          ##
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
dir_path = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.insert(0, dir_path + '/../src/')

import numpy as np
from Mite import Mite
from nexus import NexusWriter


class MiteToNexusWriter:

    def __init__(self, input_path, output_path):
        self.input_path = input_path
        self.output_path = output_path
        self.files_aux = sorted(os.listdir(self.input_path))
        self.max_token_length = 99990

    # Reduces an 1d array to a string
    def __array2string(self, array):
        s = ''

        for i in range(0, len(array)):
            s += str(array[i])

        return s

    # Reduces a 2d matrix to a string
    def __matrix2array(self, matrix, append_by_row):
        array = np.empty(0, dtype=bool)
        dim = 1

        if not append_by_row:
            matrix = matrix.tocsc()
        else:
            matrix = matrix.tocsr()
            dim = 0

        funcdict = {
            True: matrix.getrow,
            False: matrix.getcol
        }

        for i in range(0, matrix.shape[dim]):
            a = funcdict[append_by_row](i).toarray().astype(int).flatten()
            array = np.concatenate((array, a))

        return array

    # Deletes all matrix columns that all values are the same
    def __remove_equal_columns(self, matrix):
        index = np.argwhere(np.all(matrix == matrix[0, :], axis=0))
        matrix = np.delete(matrix, index, axis=1)
        return matrix

    # Writes the nexus file
    def write_nexus(self, f, w, h, min_niter, append_by_row=False, binary=False):
        dirname = 'f=' + str(f) + '_w=' + str(w) + '_h=' + str(h)
        dirname += '_min-niter=' + str(min_niter)
        nw = NexusWriter()
        run_name = []
        flattened_mites = []

        if append_by_row:
            dirname += '_append-by-row'

        if binary:
            dirname += '_binary'

        for i in range(0, len(self.files_aux)):
            basename, extension = os.path.splitext(self.files_aux[i])
            run_name.append(basename);
            m = Mite(self.input_path + self.files_aux[i], binary=binary)
            r = m.reduce_dim(f, w, h, min_niter=min_niter,
                             max_size=self.max_token_length)
            flattened_mites.append(self.__matrix2array(r, append_by_row))

        flattened_mites = np.array(flattened_mites, dtype=int)
        flattened_mites = self.__remove_equal_columns(flattened_mites);

        for i in range(0, flattened_mites.shape[0]):
            s = self.__array2string(flattened_mites[i])
            nw.add(run_name[i], 'ion_maps', 'Standard', s)

        if not os.path.exists(self.output_path + dirname):
            os.makedirs(self.output_path + dirname)

        dirname += '/'
        nw.writeFile(self.output_path + dirname + 'mite.nex')
