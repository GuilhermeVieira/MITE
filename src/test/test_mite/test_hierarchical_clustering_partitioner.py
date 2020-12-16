from unittest.mock import Mock

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
import pytest
from partitioner import HierarchicalClusteringPartitioner


def test_should_not_partition():
    partitioner = HierarchicalClusteringPartitioner(window_side_size_in_pixels=33)
    mite_mock = Mock()
    mite_mock.height = 330
    mite_mock.width = 500

    partitioner = HierarchicalClusteringPartitioner(window_side_size_in_pixels=33)

    with pytest.raises(ValueError):
        partitioner.partition_mites([mite_mock])

    mite_mock.height = 500
    mite_mock.width = 330

    with pytest.raises(ValueError):
        partitioner.partition_mites([mite_mock])


def test_hierarchical_clustering_partitioner():
    partitioner = HierarchicalClusteringPartitioner(window_side_size_in_pixels=5)
    mite_mock = Mock()
    mite_mock.height = 10
    mite_mock.width = 20
    row = [0, 9, 5, 7, 7]
    col = [0, 19, 0, 4, 5]
    data = [1.0, 2.0, 2.0, 2.0, 1.0]
    mite_mock.matrix = csr_matrix(coo_matrix((data, (row, col)), (10, 20)))

    partitioner.partition_mites([mite_mock])

    assert (mite_mock.partitioned == np.array([1.0, 0.0, 0.0, 0.0, 4.0, 1.0, 0.0, 2.0])).all()
