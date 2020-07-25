import pathlib
from enum import Enum

class Format(Enum):
    RAW = '.raw'
    MZXML = '.mzXML'

class Taxon:
    
    def __init__(self, path):
        self._path = path
        self._name = pathlib.Path(path).stem
    
    @property
    def path(self):
        return self._path

    @property
    def name(self):
        return self._name

import filesystem as fs

class SuperHirnManager:

    def __init__(self):
        self.NUM_MIN_TAXA = 2
    
    def process(self, dirpath, raw=False):
        format = Format.RAW if raw else Format.MZXML
        self.detected_taxa = self.__detect_taxa(dirpath, format)

        if len(self.detected_taxa) < self.NUM_MIN_TAXA:
            print("Invalid numeber of taxa")

        if raw:
            self.__convert_to_mzxml()
        
    def __detect_taxa(self, dirpath, format):
        detected_taxa = []

        for dir in fs.find_subdirs(dirpath):
            if fs.contains_format(dir, format):
                detected_taxa.append(Taxon(dir)) 

        return detected_taxa

    def __convert_to_mzxml(self):
        print("CALL MSCONVERT")

if __name__ == '__main__':
    dir = '/home/guilherme/MEGA/workspace/TCC_data/mzXML_files'

    manager = SuperHirnManager()
    manager.process(dir, True)

    for taxon in manager.detected_taxa:
        print(taxon.name)