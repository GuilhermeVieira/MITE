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
        if os.path.isfile(path):
            self.filepath = path
            basename, extension = os.path.splitext(self.filepath)

            if extension == '.xml':
                self.tree = ET.parse(self.filepath)
                self.root = self.tree.getroot()
                self.run_header = self.root.find("LC_MS_RUN")
                self.features = self.run_header.find("LC_MS_FEATURES")
            else:
                raise ValueError(path + ' is not a XML file!')

        #elif os.path.isdir(path):

        else:
            raise ValueError(path + " is not a file or directory!")

    def get_features_count(self):
        return int(self.run_header.get("number_of_features"))

    def get_tr_min(self):
        return float(self.run_header.get("tr_min"))

    def get_tr_max(self):
        return float(self.run_header.get("tr_max"))

    def get_mz_min(self):
        return float(self.run_header.get("m_z_min"))

    def get_mz_max(self):
        return float(self.run_header.get("m_z_max"))

    def construct_ionmap(self, binary=False):
        array_size = self.get_features_count()
        row = np.empty(array_size)
        col = np.empty(array_size)
        tr_round = np.empty(array_size)
        mz_round = np.empty(array_size)

        if binary:
            data = np.ones(array_size, dtype=bool)
        else:
            data = np.empty(array_size, dtype=int)
            # calcular quartis

        index = 0
        for feature in self.features.findall("MS1_FEATURE"):
            tr = feature.get('Tr')
            mz = feature.get('m_z')
            tr_round[index] = abs(decimal.Decimal(tr).as_tuple().exponent)
            mz_round[index] = abs(decimal.Decimal(mz).as_tuple().exponent)
            tr = float(tr) - self.get_tr_min()
            mz = float(mz) - self.get_mz_min()
            row[index] = tr
            col[index] = mz

            #if not binary:
                #data[index] = 

            index += 1

        row *= 10 ** np.amax(tr_round)
        col *= 10 ** np.amax(mz_round)
        row = np.trunc(row)
        col = np.trunc(col)

        return data, row, col, np.amax(tr_round), np.amax(mz_round)
