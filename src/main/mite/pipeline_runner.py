from abc import abstractmethod, ABC
from pathlib import Path
from typing import List

from Mite import Mite
from XMLReader import MiteLoader


class PhyloproteomicAnalyser(ABC):

    @abstractmethod
    def execute(self, mites: List[Mite]):
        pass


class MrBayesPhyloproteomicAnalyser(PhyloproteomicAnalyser):

    def execute(self, mites: List[Mite]):
        pass


class HierarchicalClusteringPhyloproteomicAnalyser(PhyloproteomicAnalyser):

    def execute(self, mites: List[Mite]):
        pass


class PipelineRunner:

    def __init__(self, mite_loader: MiteLoader,
                 phyloproteomic_analyser: PhyloproteomicAnalyser):
        self.__mites: List[Mite] = mite_loader.load_mites()
        phyloproteomic_analyser.execute(self.__mites)

    def plot_mites(self, output_directory: Path):
        for mite in self.__mites:
            mite.plot(output_directory)


if __name__ == '__main__':
    mite_loader = MiteLoader(Path('/user/src/mite/input/xml/'))
    mrbayes_phyloproteomic_analyser = MrBayesPhyloproteomicAnalyser()
    pipeline_runner = PipelineRunner(mite_loader, mrbayes_phyloproteomic_analyser)
    output_directory = Path('/user/src/mite/output/plots')
    output_directory.mkdir()
    pipeline_runner.plot_mites(output_directory)
