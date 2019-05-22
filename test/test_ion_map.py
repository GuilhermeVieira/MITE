import sys
sys.path.insert(0, '../src/')
import ion_map as im
import matplotlib.pylab as plt
import pytest
from scipy.sparse import isspmatrix_coo, isspmatrix_csr

def test_construction():
    mite = im.constructIonMap('example_pgiq.csv')
    assert isspmatrix_coo(mite)

    plt.spy(mite, markersize=0.5, aspect='auto')
    plt.show()

def test_intersection():
    mite1 = im.constructIonMap('example_pgiq.csv')
    mite2 = im.constructIonMap('example_pgiq.csv')
    intersection = im.ionMapIntersection(mite1, mite2)
    assert isspmatrix_csr(intersection)

    plt.subplot(131)
    plt.spy(mite1, markersize=0.5, aspect='auto')
    plt.subplot(132)
    plt.spy(mite2, markersize=0.5, aspect='auto')
    plt.subplot(133)
    plt.spy(intersection, markersize=0.5, aspect='auto')
    plt.show()
