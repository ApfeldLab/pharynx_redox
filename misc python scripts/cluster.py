import numpy as np
import pandas as pd
import sys
from skimage.external.tifffile import tifffile
from sklearn import preprocessing
from sklearn.externals import joblib
import feats
import sklearn.cluster

TIFF_FILE = sys.argv[1]
N_CLUSTERS = int(sys.argv[2])

img_set = tifffile.imread(TIFF_FILE)

X = feats.get_feature_matrix(img_set)
X_scaled = preprocessing.scale(X)

# clusters = sklearn.cluster.AgglomerativeClustering(n_clusters=N_CLUSTERS).fit_predict(X)
clusters = sklearn.cluster.KMeans(n_clusters=N_CLUSTERS).fit_predict(X=X_scaled)
df = pd.DataFrame(clusters)

for n in range(N_CLUSTERS):
    print(','.join(map(str, df[df[0] == n].index.values + 1)))

