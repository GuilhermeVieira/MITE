import os
dir_path = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.insert(0, dir_path + '/../src/')

import matplotlib.pylab as plt
import pytest
from Mite import Mite

input_path = dir_path + '/../input/ion_map/xml/'
files = []
mites = []

def __plotExamples(matrices, files):
    for index, (m, f) in enumerate(zip(matrices, files)):
        plt.subplot(1, len(matrices), index + 1)
        plt.spy(m, markersize=0.1, aspect='auto', color='black')
        plt.title(f)
        plt.xlabel('m/z')
        plt.ylabel('retention time')

    plt.show()

def test_construction():
    files_aux = sorted(os.listdir(input_path))

    for i in range(0, len(files_aux), 2):
        mite1 = Mite(input_path + files_aux[i])
        mite2 = Mite(input_path + files_aux[i+1])
        mites.append((mite1, mite2))
        files.append((files_aux[i], files_aux[i+1]))

    __plotExamples([mites[0][0].matrix], [files[0][0]])

def test_reduction():
    niter, w, h, f = 8, 2, 2, 0.2
    reduced = mites[0][0].reduce_dim(niter, w=w, h=h, f=f)
    __plotExamples([mites[0][0].matrix, reduced], [files[0][0],
        'Reduced matrix (niter=' + str(niter) + ', w=' + str(w) +
        ', h=' + str(h) + ', f=' + str(f) + ')'])

def test_intersection():
    intersection = mites[0][0].intersect(mites[0][1])
    count = intersection.count_nonzero()
    __plotExamples([mites[0][0].matrix, mites[0][1].matrix, intersection],
                   [files[0][0], files[0][1],
                    'Intersection (count = ' + str(count) + ')'])

def test_symmetric_diff():
    symmdiff = mites[0][0].calculate_symmdiff(mites[0][1])
    count = symmdiff.count_nonzero()
    __plotExamples([mites[0][0].matrix, mites[0][1].matrix, symmdiff],
                   [files[0][0], files[0][1],
                    'Symmetric difference (count = ' + str(count) + ')'])
