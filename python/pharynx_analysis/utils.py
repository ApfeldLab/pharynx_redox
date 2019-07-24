from collections import Counter

import numpy as np


def create_occurrence_count_tuples(l: iter) -> [(object, int)]:
    # TODO: test
    """ Given a list of things, return a list of tuples (item, nth_occurrence)

    The n_th occurrence is 0-indexed

    :param l: an iterable
    :return: [(item, nth_occurrence), ...]
    """
    count_tuples = []
    c = Counter()
    for item in l:
        c.update([item])
        count_tuples.append((item, c[item] - 1))
    return count_tuples


def figure_to_np_array(fig):
    """Given a figure, return a 3D numpy array, where the dimensions are (height, width, RGB)"""
    fig.tight_layout(pad=0)
    fig.canvas.draw()
    data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    return data
