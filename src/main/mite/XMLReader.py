##################################################################################
##                                                                              ##
##   Class that reads data from XML files based on the APML format.             ##
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


import decimal
from pathlib import Path
from typing import List, Dict

import numpy as np
import os
import xml.etree.ElementTree as ET

from scipy.sparse import coo_matrix, csr_matrix

from Mite import Mite


class XMLReader:
    def __init__(self, filepath: Path):
        if filepath.suffix == '.xml':
            self.filepath = filepath
            tree = ET.parse(self.filepath)
            root = tree.getroot()
            self.header = root.find("LC_MS_RUN")
            self.features_obj = self.header.find("LC_MS_FEATURES")
            self.features_count = int(self.header.get("number_of_features"))
        else:
            raise ValueError(str(filepath) + ' is not a XML file!')

    # Returns the minimum retention time value
    def get_tr_min(self):
        return float(self.header.get("tr_min"))

    # Returns the maximum retention time value
    def get_tr_max(self):
        return float(self.header.get("tr_max"))

    # Returns the minimum mass to charge value
    def get_mz_min(self):
        return float(self.header.get("m_z_min"))

    # Returns the maximum mass to charge value
    def get_mz_max(self):
        return float(self.header.get("m_z_max"))

    # Returns the ion map base
    def build_ion_map(self):
        array_size = self.features_count
        row = np.empty(array_size)
        col = np.empty(array_size)
        tr_round = np.empty(array_size)
        mz_round = np.empty(array_size)
        data = np.empty(array_size, dtype=float)

        index = 0
        for feature in self.features_obj.findall("MS1_FEATURE"):
            tr = feature.get('Tr')
            mz = feature.get('m_z')
            tr_round[index] = abs(decimal.Decimal(tr).as_tuple().exponent)
            mz_round[index] = abs(decimal.Decimal(mz).as_tuple().exponent)
            tr = float(tr) - self.get_tr_min()
            mz = float(mz) - self.get_mz_min()
            row[index] = tr
            col[index] = mz

            intensity = float(feature.find("LC_INFO").get("AREA"))
            data[index] = intensity

            index += 1

        row *= 10 ** np.amax(tr_round)
        col *= 10 ** np.amax(mz_round)
        row = np.trunc(row)
        col = np.trunc(col)

        return data, row, col, np.amax(tr_round), np.amax(mz_round)


class MiteLoader:

    def __init__(self, xml_directory: Path):
        self.__xml_directory: Path = xml_directory
        self.__global_quartiles: np.ndarray = self.__get_global_quartiles()

    def __get_global_quartiles(self):
        dirpath = str(self.__xml_directory)
        files = sorted(os.listdir(dirpath))
        intensities = []

        for file in files:
            if not file.startswith('.'):
                tree = ET.parse(dirpath + "/" + file)
                root = tree.getroot()
                header = root.find("LC_MS_RUN")
                features_obj = header.find("LC_MS_FEATURES")

                for feature in features_obj.findall("MS1_FEATURE"):
                    intensities.append(float(feature.find("LC_INFO").get("AREA")))

        return np.quantile(intensities, [0.25, 0.50, 0.75])

    def __build_mite(self, mite_path: Path) -> Mite:
        xml_reader = XMLReader(mite_path)
        ms_min_max_values: Dict[str, float] = {
            "tr_min": xml_reader.get_tr_min(),
            "tr_max": xml_reader.get_tr_max(),
            "mz_min": xml_reader.get_mz_min(),
            "mz_max": xml_reader.get_mz_max()
        }
        data, row, col, tr_round, mz_round = xml_reader.build_ion_map()
        shape = (int((xml_reader.get_tr_max() - xml_reader.get_tr_min()) * 10 ** tr_round),
                 int((xml_reader.get_mz_max() - xml_reader.get_mz_min()) * 10 ** mz_round))
        matrix = coo_matrix((data, (row, col)), shape)
        matrix = csr_matrix(matrix)

        return Mite(mite_path.stem, matrix, ms_min_max_values, self.__global_quartiles)

    def load_mites(self) -> List[Mite]:
        mites: List[Mite] = []

        for file in self.__xml_directory.iterdir():
            if file.is_file() and file.suffix == '.xml':
                mites.append(self.__build_mite(file))

        return mites
