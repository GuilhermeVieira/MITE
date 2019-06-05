import os
dir_path = os.path.dirname(os.path.realpath(__file__))
input_path = dir_path + '/../input/ion_map/1/'

import sys
sys.path.insert(0, dir_path + '/../src/')

import ion_map          as im
import matplotlib.pylab as plt
import pytest
from scipy.sparse            import isspmatrix_coo, isspmatrix_csr
from scipy.spatial.distance  import squareform
from scipy.cluster.hierarchy import dendrogram, linkage

def __plotExamples(matrices, files):
    for index, (m, f) in enumerate(zip(matrices, files)):
        plt.subplot(1, len(matrices), index + 1)
        plt.spy(m, markersize=0.25, aspect='auto', color='black')
        plt.title(f)
        plt.xlabel('retention time')
        plt.ylabel('m/z')
    plt.show()

def test_construction():
    filename = input_path + 'OB_COTIARA_20131206_04.csv'
    mite = im.constructIonMap(filename)
    assert isspmatrix_csr(mite)
    __plotExamples([mite], [filename])

def test_intersection():
    file1 = dir_path + '/example_pgiq_0.csv'
    file2 = dir_path + '/example_pgiq_1.csv'
    mite1 = im.constructIonMap(file1)
    mite2 = im.constructIonMap(file2)
    intersection = im.ionMapIntersection(mite1, mite2)
    count = intersection.count_nonzero()
    assert isspmatrix_csr(intersection)
    __plotExamples([mite1, mite2, intersection],
                   [file1, file2, 'Intersection (count = ' + str(count) + ')'])

def test_symmetricDiff():
    file1 = dir_path + '/example_pgiq_0.csv'
    file2 = dir_path + '/example_pgiq_1.csv'
    mite1 = im.constructIonMap(file1)
    mite2 = im.constructIonMap(file2)
    symmetric_diff = im.ionMapSymmetricDiff(mite1, mite2)
    count = symmetric_diff.count_nonzero()
    assert isspmatrix_csr(symmetric_diff)
    __plotExamples([mite1, mite2, symmetric_diff],
                   [file1, file2, 'Symmetric difference (count = ' + str(count)
                       + ')'])

def test_distMatrix():
    files = []
    mites = []
    filepath = dir_path + '/'
    for filename in os.listdir(filepath):
        if filename.endswith('.csv'):
            files.append(filename)
    files.sort()
    files.pop(0)
    for f in files:
        mite = im.constructIonMap(filepath + f)
        mites.append(mite)
    dist_matrix = im.calculateDistMatrix(mites)
    Z = linkage(squareform(dist_matrix, checks=False), 'average')
    fig = plt.figure(figsize=(25, 10))
    dn = dendrogram(Z)
    plt.show()
