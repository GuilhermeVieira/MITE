import sys
sys.path.insert(0, '../src/')
import ion_map as im
import matplotlib.pylab as plt
import pytest
from scipy.sparse import isspmatrix_coo, isspmatrix_csr

def test_construction():
    file = 'example_pgiq.csv'
    mite = im.constructIonMap(file)
    assert isspmatrix_coo(mite)

    plt.spy(mite, markersize=0.5, aspect='auto')
    plt.title(file)
    plt.xlabel('retention time')
    plt.ylabel('m/z')
    plt.show()

def test_intersection():
    file1 = 'example_pgiq.csv'
    file2 = 'example_pgiq.csv'
    mite1 = im.constructIonMap(file1)
    mite2 = im.constructIonMap(file2)
    intersection = im.ionMapIntersection(mite1, mite2)
    assert isspmatrix_csr(intersection)

    plt.subplot(131)
    plt.spy(mite1, markersize=0.5, aspect='auto')
    plt.title(file1)
    plt.xlabel('retention time')
    plt.ylabel('m/z')
    plt.subplot(132)
    plt.spy(mite2, markersize=0.5, aspect='auto')
    plt.title(file2)
    plt.xlabel('retention time')
    plt.subplot(133)
    plt.spy(intersection, markersize=0.5, aspect='auto')
    plt.title('Intersection')
    plt.xlabel('retention time')
    plt.show()
