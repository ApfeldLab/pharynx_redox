import numpy as np
from skimage import feature


def lbp_feat(img, radius=12, method='uniform'):
    n_points = 8 * radius

    lbp = feature.local_binary_pattern(img, n_points, radius, method)
    n_bins = int(lbp.max() + 1) // 4
    return np.histogram(lbp, normed=True, bins=n_bins, range=(0, n_bins))[0]


def hog_feat(img, block_norm='L2-Hys'):
    return feature.hog(img, feature_vector=True, block_norm='L1', pixels_per_cell=(8, 8))


def get_feature_matrix(img_set):
    """Returns a matrix computed by applying the given function to the img set"""
    func_list = [hog_feat, lbp_feat]
    return np.concatenate([[func(img) for img in img_set] for func in func_list], axis=1)
