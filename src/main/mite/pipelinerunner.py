from abc import abstractmethod, ABC
from pathlib import Path
from typing import List

from Mite import Mite
from XMLReader import MiteLoader
from phyloproteomicanalysis import PhyloproteomicAnalyser, MrBayesPhyloproteomicAnalyser
from multiresolution import Partitioner


class PipelineRunner:

    def __init__(self, mite_loader: MiteLoader,
                 partitioner: Partitioner,
                 phyloproteomic_analyser: PhyloproteomicAnalyser):

        self.__mites: List[Mite] = mite_loader.load_mites()
        self.__partitioned_mites: List[float] = partitioner.partition_mites(self.__mites)

        for i in range(len(self.__mites)):
            self.__mites[i].partitioned = self.__partitioned_mites[i]

        phyloproteomic_analyser.execute(self.__mites)

    def plot_mites(self, output_directory: Path):
        for mite in self.__mites:
            mite.plot(output_directory)


if __name__ == '__main__':
    mite_loader = MiteLoader(Path('/user/src/mite/input/xml/'))
    mrbayes_phyloproteomic_analyser = MrBayesPhyloproteomicAnalyser(Path('/user/src/mite/output/nexus/tres-especies'), binary=True)

    w = 16
    h = 16

    partitioner = Partitioner(w, h)

    pipeline_runner = PipelineRunner(mite_loader, partitioner, mrbayes_phyloproteomic_analyser)
    output_directory = Path('/user/src/mite/output/plots')
    # output_directory.mkdir()
    # pipeline_runner.plot_mites(output_directory)
