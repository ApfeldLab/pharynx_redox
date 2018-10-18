import numpy as np
from skimage.external.tifffile import tifffile
from sklearn.externals import joblib
import feats
import os
import argparse

# parser = argparse.ArgumentParser(description='Classify ratio images')
# parser.add_argument('filename', type=str, nargs=1,
#                     help='the TIFF file containing the ratio images')
# args = parser.parse_args()

clf = joblib.load('classifiers/adaboost_clf_03_20_18')

# TIFF_FILE = '/Users/sean/Desktop/R410_410_HD233_High_Movement_03_16_by5.tif'
TIFF_FILE = '/Users/sean/code/wormAnalysis/proc_data/r_410_470.tiff'

img_set = tifffile.imread(TIFF_FILE)

X = feats.get_feature_matrix(img_set)
np.savetxt('/Users/sean/Desktop/feats.csv', X, delimiter=',')

# predictions = np.array(clf.predict(X[:, :584]), dtype=np.uint)

# FILENAME = os.path.basename(TIFF_FILE)
# DIR_NAME = os.path.dirname(TIFF_FILE)

# np.savetxt(DIR_NAME + '/MVMT_PRED_' + os.path.splitext(FILENAME)[0] + '.csv', predictions, delimiter=',', fmt='%1u')
