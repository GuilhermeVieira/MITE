import sys
sys.path.insert(0, '../src/')
import ion_map as im
import matplotlib.pylab as plt
import pytest
from scipy.sparse import isspmatrix_coo, isspmatrix_csr

import scipy

def __plotExamples(matrices, files):
    for index, (m, f) in enumerate(zip(matrices, files)):
        plt.subplot(1, len(matrices), index + 1)
        plt.spy(m, markersize=0.5, aspect='auto')
        plt.title(f)
        plt.xlabel('retention time')
        plt.ylabel('m/z')
    plt.show()

def test_construction():
    filename = 'example_pgiq.csv'
    mite = im.constructIonMap(filename, 2, 2)
    assert isspmatrix_coo(mite)
    __plotExamples([mite], [filename])

def test_intersection():
    file1 = 'example_pgiq.csv'
    file2 = 'example_pgiq.csv'
    mite1 = im.constructIonMap(file1, 2, 2)
    mite2 = im.constructIonMap(file2, 2, 2)
    intersection, count = im.ionMapIntersection(mite1, mite2)
    assert isspmatrix_csr(intersection)
    __plotExamples([mite1, mite2, intersection],
                   [file1, file2, 'Intersection (count = ' + str(count) + ')'])

def test_symmetricDiff():
    file1 = 'example_pgiq.csv'
    file2 = 'example_pgiq.csv'
    mite1 = im.constructIonMap(file1, 2, 2)
    mite2 = im.constructIonMap(file2, 2, 2)
    symmetric_diff, count = im.ionMapSymmetricDiff(mite1, mite2)
    assert isspmatrix_csr(symmetric_diff)
    __plotExamples([mite1, mite2, symmetric_diff],
                   [file1, file2, 'Symmetric difference (count = ' + str(count)
                       + ')'])
