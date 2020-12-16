from abc import abstractmethod, ABC
from typing import List
from Mite import Mite


class Partitioner(ABC):

    @abstractmethod
    def partition_mites(self, mites: List[Mite]):
        pass
