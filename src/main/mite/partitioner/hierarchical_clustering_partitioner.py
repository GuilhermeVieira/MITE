import numpy as np
from typing import List
from Mite import Mite
from partitioner.partitioner import Partitioner
from scipy.sparse import find as find_non_zero_elements


class HierarchicalClusteringPartitioner(Partitioner):

    def __init__(self, window_side_size_in_pixels: int):
        self.__window_side_size_in_pixels = window_side_size_in_pixels

    def partition_mites(self, mites: List[Mite]) -> None:
        for mite in mites:
            mite.partitioned = self.__partition_mite(mite)

    def __partition_mite(self, mite: Mite) -> np.ndarray:
        if not self.__can_partition_mite(mite):
            raise ValueError("MITEs can not be partitioned: Height or width is not divisible by given window side size")

        partitioned_mite = self.__create_new_empty_mite_matrix(mite)

        rows, columns, values = find_non_zero_elements(mite.matrix)

        for k in range(len(values)):
            index = self.__get_new_flattened_index(mite, rows[k], columns[k]) #int(((rows[k] // self.__window_side_size_in_pixels) * mite.height / self.__window_side_size_in_pixels) + columns[k] // self.__window_side_size_in_pixels)
            partitioned_mite[index] += values[k]

        return partitioned_mite

    def __can_partition_mite(self, mite: Mite) -> bool:
        return self.__is_divisible_by(mite.width, self.__window_side_size_in_pixels) and \
               self.__is_divisible_by(mite.height, self.__window_side_size_in_pixels)

    @staticmethod
    def __is_divisible_by(dividend: int, divisor: int) -> bool:
        return dividend % divisor == 0

    def __create_new_empty_mite_matrix(self, mite: Mite) -> np.ndarray:
        return np.zeros(int((mite.width / self.__window_side_size_in_pixels) * (mite.height / self.__window_side_size_in_pixels)))

    def __get_new_flattened_index(self, mite: Mite, row: int, column: int) -> int:
        new_row = row // self.__window_side_size_in_pixels
        new_column = column // self.__window_side_size_in_pixels
        new_row_flattened_index = (new_row * mite.width) / self.__window_side_size_in_pixels
        new_index = new_row_flattened_index + new_column
        return int(new_index)
