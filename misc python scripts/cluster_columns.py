from sklearn.cluster import KMeans
import pandas as pd
import numpy as np

LOG_FILE = '/Users/sean/Desktop/2018_06_27_HD233_SAY47_SEAN/2018_06_27_HD233_SAY47_SEAN.LOG'

def cluster_cols(logfile_path, strains, n_channels=3):
    df = pd.read_csv(logfile_path, quotechar='"')
    x = df[' "Plane Stage Position X"'][1::n_channels]
    y = df[' "Plane Stage Position Y"'][1::n_channels]
    x = np.asarray(x).reshape(-1, 1)
    kmeans = KMeans(n_clusters=len(strains)).fit(x)
    cols = pd.unique(kmeans.labels_)
    print(cols)


if __name__ == '__main__':
    cluster_cols(LOG_FILE, ['HD233', 'SAY47', 'HD233'])