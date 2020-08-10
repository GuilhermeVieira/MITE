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
from typing import Dict

import numpy as np
import matplotlib.pylab as plt
from scipy.sparse import csr_matrix


class Mite:

    def __init__(self, name: str,
                 matrix: csr_matrix, ms_boundary_values: Dict[str, float], global_quartiles: np.quantile):
        self.__name = name
        self.__matrix: csr_matrix = matrix
        self.__ms_boundary_values: Dict[str, float] = ms_boundary_values
        self.__global_quartiles = global_quartiles

    def plot(self, output_directory: Path):
        plt.spy(self.__matrix, markersize=0.1, aspect='auto', color='black')
        plt.title(self.__name)
        plt.xlabel('m/z')
        plt.ylabel('retention time')
        plt.savefig(f'{output_directory}/{self.__name}.png')
        plt.clf()
        plt.cla()
        plt.close()
