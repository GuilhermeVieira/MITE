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

import numpy as np
import os
import sys

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, dir_path + '/../../src/')

import multiresolution as mr
from Mite import Mite
from nexus import NexusWriter

class MiteToNexusWriter:

    def __init__(self, input_path, output_path, binary):
        self.input_path = input_path
        self.output_path = output_path
        self.binary = binary
        self.files_aux = sorted(os.listdir(self.input_path))
        self.max_token_length = 99990
        self.mites, self.run_name = self.__load_mites()

    # Loads all the mites
    def __load_mites(self):
        print("Loading ion maps...")

        mites = []
        run_name = []

        for i in range(0, len(self.files_aux)):
            basename, extension = os.path.splitext(self.files_aux[i])

            if not basename.startswith('.'):
                run_name.append(basename);
                mites.append(
                    Mite(self.input_path + self.files_aux[i], binary=self.binary)
                )

        return mites, run_name

    # Transforms an 1d array in a string
    def __array2string(self, array):
        s = ''

        for i in range(0, len(array)):
            s += str(array[i])

        return s

    # Reduces a 2d matrix to a string
    def __matrix2array(self, matrix):
        array = np.empty(0)
        dim = 1
        matrix = matrix.tocsc()

        for i in range(0, matrix.shape[dim]):
            a = matrix.getcol(i).toarray().astype(int).flatten()
            array = np.concatenate((array, a))

        return array

    # Deletes all matrix columns that all values are the same
    def __remove_equal_columns(self, matrix):
        index = np.argwhere(np.all(matrix == matrix[0, :], axis=0))
        matrix = np.delete(matrix, index, axis=1)
        return matrix

    # Sets the dir name to write the nexus file
    def __set_dirname(self, w, h, f):
        dirname = 'w=' + str(w) + '__h=' + str(h)

        if self.binary:
            dirname += '__binary' + '__f=' + str(f)

        return dirname

    # Flattens the mites
    def __flatten_mites(self, w, h, f):
        flattened_mites = []

        for m in self.mites:
            r = mr.reduce_dim(
                m.matrix, w, h, m.binary, max_size=self.max_token_length, f=f
            )
            flattened_mites.append(self.__matrix2array(r))

        flattened_mites = np.array(flattened_mites, dtype=int)
        flattened_mites = self.__remove_equal_columns(flattened_mites);

        return flattened_mites

    # Writes the nexus file
    def write_nexus(self, w, h, f=0.0):
        nw = NexusWriter()
        dirname = self.__set_dirname(w, h, f)
        flattened_mites = self.__flatten_mites(w, h, f)

        for i in range(0, flattened_mites.shape[0]):
            s = self.__array2string(flattened_mites[i])
            nw.add(self.run_name[i], 'ion_maps', 'Standard', s)

        if not os.path.exists(self.output_path + dirname):
            os.makedirs(self.output_path + dirname)

        dirname += '/'
        nw.writeFile(self.output_path + dirname + 'mite.nex')
