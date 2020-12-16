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

import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

from discretizer import QuartileDiscretizer, BinaryDiscretizer, Discretizer
from matplotlib import pyplot as plt
from Mite import Mite
from nexus import NexusWriter

import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
from rpy2.robjects.packages import importr

class PhyloproteomicAnalyser(ABC):

    @abstractmethod
    def execute(self, mites: List[Mite]):
        pass


class HierarchicalClusteringPhyloproteomicAnalyser(PhyloproteomicAnalyser):

    def execute(self, mites: List[Mite]):
        heat_map = []
        labels = []

        for mite in mites:
            labels.append(mite.name)
            heat_map.append(mite.partitioned)

        highest_intensity_sum = 0
        for i in range(len(heat_map)):
            for j in range(len(heat_map[i])):
                if heat_map[i][j] > highest_intensity_sum:
                    highest_intensity_sum = heat_map[i][j]

        heat_map = np.array(heat_map)

        # for i in range(len(heat_map)):
        #     for j in range(len(heat_map[i])):
        #         heat_map[i][j] /= highest_intensity_sum

        # for i in range(len(heat_map1)):
        #     for j in range(len(heat_map[i])):
        #         heat_map[i][j] = 1.0 if heat_map[i][j] != 0.0 else 0.0

        #heat_map = self.__remove_equal_columns(heat_map)

        self.__build_hierarchical_clustering(heat_map, labels)

        base = importr("base")
        pvclust = importr("pvclust")
        stats = importr("stats")
        graphics = importr("graphics")
        parallel = importr("parallel")

        numpy2ri.activate()

        result = stats.hclust(stats.dist(heat_map, method="euclidean"), method="average")
        graphics.plot(result, labels=labels)

        #result = pvclust.pvclust(heat_map.transpose(), nboot=1000, method_dist="canberra", method_hclust="average", parallel=True)
        #graphics.plot(result)
        numpy2ri.activate()

    # def __remove_equal_columns(self, matrix):
    #     index = np.argwhere(np.all(matrix < 0.5, axis=0))
    #     index = np.argwhere(np.any(matrix > 0.2, axis=0))
    #     print("Percentage of removed itens: " + str(len(index)/len(matrix[0])))
    #     matrix = np.delete(matrix, index, axis=1)
    #     return matrix

    def __build_hierarchical_clustering(self, occurrence_matrix, labels):
        metric = 'euclidean'
        method = 'average'

        Z = linkage(occurrence_matrix, method=method, metric=metric, optimal_ordering=False)

        print(Z)
        plt.figure(figsize=(10, 10))
        plt.title("Snakes Hierarchical Clustering Dendrogram - " + metric + ", " + " " + method)
        plt.xlabel("Snakes")
        plt.ylabel("Distance")

        SE = dendrogram(Z, labels=labels, leaf_rotation=0.0, leaf_font_size=9.0)
        save_name = "dendogram-snakes-" + metric + "-" + method
        plt.savefig(save_name + ".svg", format="svg")


class MrBayesPhyloproteomicAnalyser(PhyloproteomicAnalyser):

    def __init__(self, output_path: Path, binary: bool = False):
        self.__binary: bool = binary
        self.__output_path: Path = output_path

    def execute(self, mites: List[Mite]):
        discretizer: Discretizer = BinaryDiscretizer() if self.__binary else QuartileDiscretizer(mites[0].global_quartiles)
        discretized_mites: Dict[str, List[int]] = {}

        for mite in mites:
            discretized_mites[mite.name] = discretizer.discretize(mite.partitioned)

        discretized_mites = self.__remove_equal_columns(discretized_mites)

        self.__write_nexus(discretized_mites)

    # Transforms an 1d array in a string
    def __array2string(self, array: List) -> str:
        s: str = ''

        for i in range(0, len(array)):
            s += str(array[i])

        return s

    # Deletes all matrix columns that all values are the same
    def __remove_equal_columns(self, mites: Dict[str, List[int]]):
        matrix = np.array([value for value in mites.values()])
        index = np.argwhere(np.all(matrix == matrix[0, :], axis=0))
        matrix = np.delete(matrix, index, axis=1)

        a = {}
        index = 0
        for key in mites.keys():
            a[key] = matrix[index].tolist()
            index += 1

        return a

    # Writes the nexus file
    def __write_nexus(self, mites: Dict[str, List[int]]):
        nw: NexusWriter = NexusWriter()

        for name, array in mites.items():
            nw.add(name, 'ion_maps', 'Standard', self.__array2string(array))

        if not self.__output_path.is_dir():
            self.__output_path.mkdir()

        nw.writeFile(str(self.__output_path) + '/mite.nex')


    # def mink(self, X):
    #     return pdist(X, 'minkowski', p=2)
