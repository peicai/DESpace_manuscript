import pandas as pd
import pandas.util.testing as tm
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime
import pandas as pd
import pandas.util.testing as tm
import random, torch
import stlearn as st
def stLearn(path):
   # counts_files = glob.glob(path+ "mel" + '[0-4]' + "_rep1_count.csv")
    #counts_files = glob.glob(path+  "*[0-2]"+"_count.csv")
    counts_files = glob.glob(path+ '151[5-6][0-7][0-9]' + "_count.csv")
    loc = np.ravel([[x[:-10] + '_location.csv'] for x in counts_files])
    sample_output = np.ravel([[x[:-10] + '_stLearn_results.csv'] for x in counts_files])
    times_output = np.ravel([[x[:-10] + '_stLearn_runtime.csv'] for x in counts_files])
    df = pd.DataFrame({'counts':counts_files, 'meta':loc, 'output':sample_output, 'time':times_output})
    for index, row in df.iterrows():
        np.random.seed(123)
        random.seed(123)
        torch.manual_seed(123)
        counts = pd.read_csv(row['counts'], index_col=0) # load counts
        print(counts.head())
        sample_info = pd.read_csv(row['meta'], index_col=0) # load meta data with spatial coordinates
        sample_info = sample_info.iloc[:,0:2]
        print(sample_info.head())
        sample_info = sample_info.rename(columns={'col': 'imagecol', 'row': 'imagerow'})
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
        st.tl.clustering.kmeans(data_SME_LIBD,n_clusters=5, use_data="X_pca", random_state = int('123'), key_added="X_pca_kmeans")
        st.pl.cluster_plot(data_SME_LIBD, use_label="X_pca_kmeans")
        obsm_data=pd.DataFrame(data_SME_LIBD.obs)
        print(obsm_data.head())
        end = datetime.now()
        time = pd.DataFrame({end - start})
        time.to_csv(row['time']) # Save spatial results
        # total time taken
        print(f"Runtime of the program is {end - start}")
        obsm_data.to_csv(row['output']) # Save spatial results
        print("Finished =", row['counts'])
    return


## Change 'counts_files' path and 'n_clusters=' in the function

path = "/home/peicai/project/SV/data/real_data_analysis/LIBD/"
stLearn(path)
path = "./DESpace_data/Simulation/Output/LIBD/bottom_patch/probs_0.5_0.9FALSE"
stLearn(path)

path1 = "./DESpace_data/Simulation/Output/LIBD/circle_patch/probs_0.5_0.9FALSE"
stLearn(path1)

path3 = "./DESpace_data/Simulation/Output/LIBD/Manual_clusters_patch/probs_0.5_0.9FALSE"
stLearn(path3)

path1 = "./DESpace_data/Simulation/Output/LIBD/mixture_patch/probs_0.5_0.9FALSE"
stLearn(path1)

path3 = "./DESpace_data/Simulation/Output/LIBD/mixture_reverse_patch/probs_0.5_0.9FALSE"
stLearn(path3)

path = "./DESpace_data/Simulation/Output/LIBD/bottom_patch/probs_0.6_0.9FALSE"
stLearn(path)

path1 = "./DESpace_data/Simulation/Output/LIBD/circle_patch/probs_0.6_0.9FALSE"
stLearn(path1)

path3 = "./DESpace_data/Simulation/Output/LIBD/Manual_clusters_patch/probs_0.6_0.9FALSE"
stLearn(path3)

path1 = "./DESpace_data/Simulation/Output/LIBD/mixture_patch/probs_0.6_0.9FALSE"
stLearn(path1)

path3 = "./DESpace_data/Simulation/Output/LIBD/mixture_reverse_patch/probs_0.6_0.9FALSE"
stLearn(path3)



## n_clusters=3
path1 = "./DESpace_data/Simulation/Output/melanoma/mixture_patch/probs_0.5_0.9FALSE/"
stLearn(path1)

path2 = "./DESpace_data/Simulation/Output/melanoma/mixture_reverse_patch/probs_0.5_0.9FALSE/"
stLearn(path2)

path3 = "./DESpace_data/Simulation/Output/melanoma/right_patch/probs_0.5_0.9FALSE/"
stLearn(path3)

path4 = "./DESpace_data/Simulation/Output/melanoma/circle_patch/probs_0.5_0.9FALSE/"
stLearn(path4)

path3 = "./DESpace_data/Simulation/Output/melanoma/right_patch/probs_0.6_0.9FALSE/"
stLearn(path3)

path4 = "./DESpace_data/Simulation/Output/melanoma/circle_patch/probs_0.6_0.9FALSE/"
stLearn(path4)


path1 = "./DESpace_data/Simulation/Output/melanoma/mixture_patch/probs_0.6_0.6FALSE/"
stLearn(path1)

path1 = "./DESpace_data/Simulation/Output/melanoma/mixture_reverse_patch/probs_0.6_0.6FALSE/"
stLearn(path1)

path3 = "./DESpace_data/Simulation/Output/melanoma/BayesSpace_clusters_patch/probs_0.5_0.9FALSE/"
stLearn(path3)

path4 = "./DESpace_data/Simulation/Output/melanoma/BayesSpace_clusters_patch/probs_0.6_0.9FALSE/"
stLearn(path4)


## mouse cerebellum
# n_clusters=4
path = "./DESpace_data/Simulation/Output/mouse cerebellum/BayesSpace_clusters_patch/probs_0.6_0.9FALSE/"
stLearn(path)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/BayesSpace_clusters_patch/probs_0.5_0.9FALSE/"
stLearn(path)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/right_patch/probs_0.6_0.9FALSE/"
stLearn(path)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/right_patch/probs_0.5_0.9FALSE/"
stLearn(path)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/circle_patch/probs_0.5_0.9FALSE/"
stLearn(path)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/circle_patch/probs_0.6_0.9FALSE/"
stLearn(path)

path1 = "./DESpace_data/Simulation/Output/mouse cerebellum/mixture_patch/probs_0.6_0.6FALSE/"
stLearn(path1)

path1 = "./DESpace_data/Simulation/Output/mouse cerebellum/mixture_reverse_patch/probs_0.6_0.6FALSE/"
stLearn(path1)

path1 = "./DESpace_data/Simulation/Output/mouse cerebellum/mixture_patch/probs_0.5_0.9FALSE/"
stLearn(path1)

path1 = "./DESpace_data/Simulation/Output/mouse cerebellum/mixture_reverse_patch/probs_0.5_0.9FALSE/"
stLearn(path1)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/"
# counts_files = glob.glob(path+  "100_count.csv")
stLearn(path)
