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

    # Reduces a 2d matrix into a string
    def __matrix2string(self, matrix, mode, f=0.25):
        if mode == 'mz':
            matrix = matrix.tocsc()

        s = ''
        dimdict = {
            'tr': 0,
            'mz': 1
        }
        funcdict = {
            'tr': matrix.getrow,
            'mz': matrix.getcol
        }

        for i in range(0, matrix.shape[dimdict[mode]]):
            ratio = funcdict[mode](i).count_nonzero()
            ratio /= matrix.shape[dimdict[mode]]

            if ratio >= f:
                s += '1'
            else:
                s += '0'

        return s

    # Writes the nexus file
    def write_nexus(self, niter, f, mode='tr', all=False):
        step = 2
        dirname = 'mode-' + mode
        dirname += '_f-' + str(f)
        dirname += '_niter-' + str(niter)
        nw = NexusWriter()

        if all:
            step = 1
            dirname += '_all'

        for i in range(0, len(self.files_aux), step):
            basename, extension = os.path.splitext(self.files_aux[i])
            m = Mite(self.input_path + self.files_aux[i])
            r = m.reduce_dim(niter=niter, f=f)
            nw.add(basename, 'ion_maps', 'Standard',
                   self.__matrix2string(r, mode, f=f))

        if not os.path.exists(self.output_path + dirname):
            os.makedirs(self.output_path + dirname)

        dirname += '/'
        nw.writeFile(self.output_path + dirname + 'mite.nex')
