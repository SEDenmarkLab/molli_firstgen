import numpy as np


def distance(p1, p2):
    """
    Returns Euclidean distance between two points
    """
    _p1, _p2 = np.array(p1), np.array(p2)
    return np.sqrt(np.sum((_p1 - _p2) ** 2))