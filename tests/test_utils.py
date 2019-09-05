import os
import pickle
import numpy as np

from pharynx_analysis import profile_processing as pp, utils


class TestUtils:
    def test_create_occurrence_count_tuples(self):
        l = ["410", "470", "410", "470"]
        expected = [("410", 0), ("470", 0), ("410", 1), ("470", 1)]

        assert utils.create_occurrence_count_tuples(l) == expected
