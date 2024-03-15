
import numpy as np
from dipy.segment.clustering import QuickBundles
from dipy.io.pickles import save_pickle
import scipy.io as spio
from dipy.segment.metric import ResampleFeature
from dipy.segment.metric import VectorOfEndpointsFeature
from dipy.segment.metric import AveragePointwiseEuclideanMetric

mat = spio.loadmat('pl_cur.mat', squeeze_me=True)
streams = mat['pl_cur']
feature = ResampleFeature(nb_points=124)
metric = AveragePointwiseEuclideanMetric(feature=feature)
qb = QuickBundles(threshold=4, metric=metric)
pl = [i for i in streams]
clusters = qb.cluster(pl)
len(clusters)
print("Nb. clusters:", len(clusters))
print("Cluster sizes:", list(map(len, clusters)))
pli = [i.indices for i in clusters]
pli_array = np.array([np.array(i) for i in pli])
spio.savemat('pli_array.mat', {'pli_array': pli_array})
pl_centroid = [i.centroid for i in clusters]
pl_centroid_array = np.array([np.array(i) for i in pl_centroid])
spio.savemat('pl_centroid_array.mat', {'pl_centroid_array': pl_centroid_array})