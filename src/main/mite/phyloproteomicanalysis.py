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

from abc import abstractmethod, ABC
from pathlib import Path
from typing import List, Dict
from discretizer import QuartileDiscretizer, BinaryDiscretizer, Discretizer
from Mite import Mite
from nexus import NexusWriter


class PhyloproteomicAnalyser(ABC):

    @abstractmethod
    def execute(self, mites: List[Mite]):
        pass


class HierarchicalClusteringPhyloproteomicAnalyser(PhyloproteomicAnalyser):

    def execute(self, mites: List[Mite]):
        pass


class MrBayesPhyloproteomicAnalyser(PhyloproteomicAnalyser):

    def __init__(self, output_path: Path, binary: bool = False):
        self.__binary: bool = binary
        self.__output_path: Path = output_path

    def execute(self, mites: List[Mite]):
        discretizer: Discretizer = BinaryDiscretizer() if self.__binary else QuartileDiscretizer(mites[0].global_quartiles)
        discretized_mites: Dict[str, List[int]] = {}

        for mite in mites:
            discretized_mites[mite.name] = discretizer.discretize(mite.partitioned)

        self.__write_nexus(discretized_mites)

    # Transforms an 1d array in a string
    def __array2string(self, array: List) -> str:
        s: str = ''

        for i in range(0, len(array)):
            s += str(array[i])

        return s

    # Writes the nexus file
    def __write_nexus(self, mites: Dict[str, List[int]]):
        nw: NexusWriter = NexusWriter()

        for name, array in mites.items():
            nw.add(name, 'ion_maps', 'Standard', self.__array2string(array))

        if not self.__output_path.is_dir():
            self.__output_path.mkdir()

        nw.writeFile(str(self.__output_path) + '/mite.nex')
