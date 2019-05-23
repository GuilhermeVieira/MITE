import os
dir_path = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.insert(0, dir_path + '/../src/')

import ion_map as im
import matplotlib.pylab as plt
import pytest
from scipy.sparse import isspmatrix_coo, isspmatrix_csr

def __plotExamples(matrices, files):
    for index, (m, f) in enumerate(zip(matrices, files)):
        plt.subplot(1, len(matrices), index + 1)
        plt.spy(m, markersize=0.5, aspect='auto')
        plt.title(f)
        plt.xlabel('retention time')
        plt.ylabel('m/z')
    plt.show()

def test_construction():
    filename = dir_path + '/example_pgiq.csv'
    mite = im.constructIonMap(filename, 2, 2)
    assert isspmatrix_coo(mite)
    __plotExamples([mite], [filename])

def test_intersection():
    file1 = dir_path + '/example_pgiq.csv'
    file2 = dir_path + '/example_pgiq.csv'
    mite1 = im.constructIonMap(file1, 2, 2)
    mite2 = im.constructIonMap(file2, 2, 2)
    intersection = im.ionMapIntersection(mite1, mite2)
    count = intersection.count_nonzero()
    assert isspmatrix_csr(intersection)
    __plotExamples([mite1, mite2, intersection],
                   [file1, file2, 'Intersection (count = ' + str(count) + ')'])

def test_symmetricDiff():
    file1 = dir_path + '/example_pgiq.csv'
    file2 = dir_path + '/example_pgiq.csv'
    mite1 = im.constructIonMap(file1, 2, 2)
    mite2 = im.constructIonMap(file2, 2, 2)
    symmetric_diff = im.ionMapSymmetricDiff(mite1, mite2)
    count = symmetric_diff.count_nonzero()
    assert isspmatrix_csr(symmetric_diff)
    __plotExamples([mite1, mite2, symmetric_diff],
                   [file1, file2, 'Symmetric difference (count = ' + str(count)
                       + ')'])

def test_distMatrix():
    files = []
    mites = []
    filepath = dir_path + '/../input/ion_map/'
    for filename in os.listdir(filepath):
        files.append(filename)
    files.sort()
    for i in range(0, len(files), 2):
        mite1 = im.constructIonMap(filepath + files[i], 2, 2)
        mite2 = im.constructIonMap(filepath + files[i+1], 2, 2)
        intersection = im.ionMapIntersection(mite1, mite2)
        mites.append(intersection)
    dist_matrix = im.calculateDistMatrix(mites)
