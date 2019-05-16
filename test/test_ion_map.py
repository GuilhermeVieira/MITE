import sys
sys.path.insert(0, '../src/')
import ion_map as im
import matplotlib.pylab as plt
import pytest

def test_ion_map():
    mite = im.constructIonMap('example_pgiq.csv')
    plt.spy(mite, markersize=0.5, aspect='auto')
    plt.show()
