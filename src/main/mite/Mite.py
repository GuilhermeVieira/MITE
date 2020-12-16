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
from pathlib import Path
from typing import Dict, List
from scipy.sparse import csr_matrix
import numpy as np
import matplotlib.pylab as plt


class Mite:

    def __init__(self, name: str, matrix: csr_matrix, ms_min_max_values: Dict[str, float], global_quartiles: np.quantile):
        self.__name: str = name
        self.__matrix: csr_matrix = matrix
        self.__ms_min_max_values: Dict[str, float] = ms_min_max_values
        self.__global_quartiles: np.quantile = global_quartiles
        self.__partitioned: List[float] = []

    @property
    def name(self):
        return self.__name

    @property
    def matrix(self):
        return self.__matrix

    @property
    def width(self):
        return self.__matrix.shape[1]

    @property
    def height(self):
        return self.__matrix.shape[0]

    @property
    def global_quartiles(self):
        return self.__global_quartiles

    @property
    def partitioned(self):
        if len(self.__partitioned) == 0:
            raise ValueError("Mite was not partitioned yet")
        return self.__partitioned

    @partitioned.setter
    def partitioned(self, flattened_mite: List[float]):
        self.__partitioned = flattened_mite

    def plot(self, output_directory: Path):
        plt.spy(self.__matrix, markersize=0.1, aspect='auto', color='black')
        plt.title(self.__name)
        plt.xlabel('m/z')
        plt.ylabel('retention time')
        plt.savefig(f'{str(output_directory)}/{self.__name}.png')
        plt.clf()
        plt.cla()
        plt.close()
