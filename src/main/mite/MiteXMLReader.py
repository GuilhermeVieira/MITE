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
import numpy as np
import os
import xml.etree.ElementTree as ET


class MiteXMLReader:
    def __init__(self, path):
        self.filepath = path
        basename, extension = os.path.splitext(self.filepath)

        if extension == '.xml':
            tree = ET.parse(self.filepath)
            root = tree.getroot()
            self.header = root.find("LC_MS_RUN")
            self.features_obj = self.header.find("LC_MS_FEATURES")
            self.features_count = int(self.header.get("number_of_features"))
        else:
            raise ValueError(path + ' is not a XML file!')

    # Returns the three local quartiles of intensity values
    # def __get_local_quartiles(self):
    #     intensities = []
    #
    #     for feature in self.features_obj.findall("MS1_FEATURE"):
    #         intensities.append(float(feature.find("LC_INFO").get("AREA")))
    #
    #     return np.quantile(intensities, [0.25, 0.50, 0.75])

    # Returns the three global quartiles of intensity values
    def __get_global_quartiles(self):
        dirpath = os.path.dirname(self.filepath)
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

    # Constructs and returns the ion map base
    def construct_ionmap(self, binary=False):
        array_size = self.features_count
        row = np.empty(array_size)
        col = np.empty(array_size)
        tr_round = np.empty(array_size)
        mz_round = np.empty(array_size)
        data = None
        quartiles = None

        if binary:
            data = np.ones(array_size, dtype=bool)
        else:
            data = np.empty(array_size, dtype=int)
            quartiles = self.__get_global_quartiles()

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

            if not binary:
                intensity = float(feature.find("LC_INFO").get("AREA"))

                if intensity < quartiles[0]:
                    data[index] = 1
                elif intensity < quartiles[1]:
                    data[index] = 2
                elif intensity < quartiles[2]:
                    data[index] = 3
                else:
                    data[index] = 4

            index += 1

        row *= 10 ** np.amax(tr_round)
        col *= 10 ** np.amax(mz_round)
        row = np.trunc(row)
        col = np.trunc(col)

        return data, row, col, np.amax(tr_round), np.amax(mz_round)
