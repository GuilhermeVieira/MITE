import os
dir_path = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.insert(0, dir_path + '/../src/mite/')

import matplotlib.pylab as plt
import pytest
from Mite import Mite

input_path = dir_path + '/../input/ion_map/xml/'
mites = []

def test_construction():
    for f in sorted(os.listdir(input_path)):
        basename, extension = os.path.splitext(f)

        if extension == '.xml':
            mites.append(Mite(input_path + f, binary=True))

def test_plot():
    for mite in mites:
        mite.plot()
