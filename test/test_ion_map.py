import os
dir_path = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.insert(0, dir_path + '/../src/')

import ion_map          as im
import matplotlib.pylab as plt
import pytest
from scipy.sparse            import isspmatrix_coo, isspmatrix_csr
from scipy.spatial.distance  import squareform
from scipy.cluster.hierarchy import dendrogram, linkage

input_path = dir_path + '/../input/ion_map/1/'
files = []
mites = []

def __plotExamples(matrices, files):
    for index, (m, f) in enumerate(zip(matrices, files)):
        plt.subplot(1, len(matrices), index + 1)
        plt.spy(m, markersize=0.1, aspect='auto', color='black')
        plt.title(f)
        plt.xlabel('retention time')
        plt.ylabel('m/z')
    plt.show()

def test_construction():
    files_aux = sorted(os.listdir(input_path))
    for i in range(0, len(files_aux), 2):
        mite1 = im.constructIonMap(input_path + files_aux[i])
        mite2 = im.constructIonMap(input_path + files_aux[i+1])
        mites.append((mite1, mite2))
        files.append((files_aux[i], files_aux[i+1]))
    __plotExamples([mites[0][0]], [files[0][0]])

def test_intersection():
    intersec = im.ionMapIntersection(mites[0][0], mites[0][1])
    count = intersec.count_nonzero()
    __plotExamples([mites[0][0], mites[0][1], intersec],
                   [files[0][0], files[0][1],
                    'Intersection (count = ' + str(count) + ')'])

def test_symmetricDiff():
    sym_diff = im.ionMapSymmetricDiff(mites[0][0], mites[0][1])
    count = sym_diff.count_nonzero()
    __plotExamples([mites[0][0], mites[0][1], sym_diff],
                   [files[0][0], files[0][1],
                    'Symmetric difference (count = ' + str(count) + ')'])

def test_distMatrix():
    dist_matrix = im.calculateDistMatrix(mites)
    Z = linkage(squareform(dist_matrix, checks=False),
                'average',
                optimal_ordering=True)
    fig = plt.figure(figsize=(25, 10))
    dn = dendrogram(Z)
    plt.show()
