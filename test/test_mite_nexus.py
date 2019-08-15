import os
dir_path = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.insert(0, dir_path + '/../src/')

import numpy as np
import pytest
from mite import Mite
from nexus import NexusWriter

input_path = dir_path + '/../input/ion_map/xml/'

def __array2string(array):
    s = ''

    for i in range(0, len(array)):
        s += str(int(array[i]))

    return s

def test_mite2nex():
    nw = NexusWriter()
    files_aux = sorted(os.listdir(input_path))
    print('\n')

    for i in range(0, len(files_aux)):
        print('Processing file ' + files_aux[i])
        basename, extension = os.path.splitext(files_aux[i])
        m = Mite(input_path + files_aux[i])
        r = m.reduce_dim(niter=20, f=0.2)
        nw.add(basename, 'ion_maps', 'Standard',
                __array2string(r.toarray().flatten()))

    nw.writeFile('mite.nex')
