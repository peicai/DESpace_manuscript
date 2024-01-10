# pip install python-igraph
# pip install louvain
# pip install tables
# #pip install git+https://github.com/BiomedicalMachineLearning/stLearn.git@a5c1ee6f5ef3c97d89bab46cf88054889a0f2cac
# pip install -U stlearn

import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime
#import pandas.util.testing as tm
import pandas.testing as tm
import stlearn as st
import random, torch

paths = ["./Simulation/Output/151507/split_mixture_patch_{}".format(i) for i in [2, 4, 6, 8, 10, 12, 14]]

counts_files = np.ravel([path + "/probs_0.5_0.9FALSE/py_count.csv" for path in paths])
loc = np.ravel([path + "/probs_0.5_0.9FALSE/py_meta.csv" for path in paths])
sample_output = np.ravel([path + "/probs_0.5_0.9FALSE/stLeaern_results.csv" for path in paths])
times_output = np.ravel([path + '/probs_0.5_0.9FALSE/stLearn_runtime.csv' for path in paths])
num_cluster = [2, 4, 6, 8, 10, 12, 14]
resolution = [0.4, 0.7, 1, 1.2, 1.5, 1.7, 2]
df = pd.DataFrame({'counts':counts_files, 'meta':loc, 'output':sample_output, 'time': times_output,
                   'num_cluster': num_cluster, 'res': resolution})

for index, row in df.iterrows():
    r_seed=t_seed=n_seed=123
    #Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    counts = pd.read_csv(row['counts'], index_col=0) # load counts
    #print(counts.head())
    sample_info = pd.read_csv(row['meta'], index_col=0) # load meta data with spatial coordinates
    #print(sample_info.head())
    print(sample_info.columns)
    sample_info = sample_info.iloc[:,5:7]
    #print(sample_info)
    #sample_info = sample_info.rename(columns={'col': 'imagecol', 'row': 'imagerow'})
    print(sample_info.head())
    adata = st.create_stlearn(count=counts.T,spatial=sample_info,library_id="LIBD", scale=1,background_color="white")
    st.pp.filter_genes(adata,min_cells=0)
    start= datetime.now()
    st.pp.normalize_total(adata)
    st.pp.log1p(adata)
    data = adata
    st.em.run_pca(data,n_comps=50)
    # K-means clustering on stSME normalised PCA
    data_SME_LIBD = data.copy()
    st.tl.clustering.kmeans(data_SME_LIBD,n_clusters=row['num_cluster'], use_data="X_pca", key_added="X_pca_kmeans")
    st.pl.cluster_plot(data_SME_LIBD, use_label="X_pca_kmeans")
    # louvain clustering on stSME normalised data
    st.pp.neighbors(data_SME_LIBD,n_neighbors=17,use_rep='X_pca')
    st.tl.clustering.louvain(data_SME_LIBD, resolution=row['res'])
    st.pl.cluster_plot(data_SME_LIBD,use_label="louvain")
    obsm_data=pd.DataFrame(data_SME_LIBD.obs)
    print(obsm_data.head())
    end = datetime.now()
    # total time taken
    print(f"Runtime of the program is {end - start}")
    time = pd.DataFrame({end - start})
    time.to_csv(row['time']) # Save spatial results
    obsm_data.to_csv(row['output']) # Save spatial results
    print("Finished =", row['counts'])
