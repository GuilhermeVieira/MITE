from typing import List

import numpy as np
from abc import ABC, abstractmethod


class Discretizer(ABC):

    @abstractmethod
    def discretize_number(self, number: float) -> int:
        pass

    def discretize(self, array: List[float]) -> List[int]:
        discretized_array: List[int] = []

        for number in array:
            discretized_array.append(self.discretize_number(number))

        return discretized_array


class BinaryDiscretizer(Discretizer, ABC):

    def discretize_number(self, number: float) -> int:
        return 1 if number else 0


class QuartileDiscretizer(Discretizer):

    def __init__(self, quartiles: np.ndarray):
        self.__quartiles: np.ndarray = quartiles

    def discretize_number(self, number: float) -> int:
        super(QuartileDiscretizer, self).discretize_number(number)
        if number == 0.0:
            return 0
        elif number < self.__quartiles[0]:
            return 1
        elif number < self.__quartiles[1]:
            return 2
        elif number < self.__quartiles[2]:
            return 3
        else:
            return 4
