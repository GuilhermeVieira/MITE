##################################################################################
##                                                                              ##
##   Class that implements the ion intensity maps extracted from mass           ##
##   spectrometry raw files.                                                    ##
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

import matplotlib.pylab as plt
from scipy.sparse import coo_matrix, csr_matrix

from mite.MiteXMLReader import MiteXMLReader

class Mite:

    # Class constructor
    def __init__(self, filepath, binary=False):
        xml_reader = MiteXMLReader(filepath)

        self.filepath = filepath
        self.binary = binary
        self.tr_min = xml_reader.get_tr_min()
        self.tr_max = xml_reader.get_tr_max()
        self.mz_min = xml_reader.get_mz_min()
        self.mz_max = xml_reader.get_mz_max()

        data, row, col, tr_round, mz_round = xml_reader.construct_ionmap(binary=binary)
        shape = (int((self.tr_max - self.tr_min) * 10 ** tr_round),
                int((self.mz_max - self.mz_min) * 10 ** mz_round))
        self.matrix = coo_matrix((data, (row, col)), shape)
        self.matrix = csr_matrix(self.matrix)

    # Plots the mite
    def plot(self):
        plt.spy(self.matrix, markersize=0.1, aspect='auto', color='black')
        plt.title(self.filepath)
        plt.xlabel('m/z')
        plt.ylabel('retention time')
        plt.savefig(self.filepath + '.png')
        plt.clf()
        plt.cla()
        plt.close()
