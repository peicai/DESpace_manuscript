import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime
import glob
#pip install NaiveDE
#pip install SpatialDE
import NaiveDE
import SpatialDE
from math import sqrt
from joblib import Parallel, delayed
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import spatial

def xpercent_scale():
    gca().set_xticklabels(['{:.0f}%'.format(x*100) for x in gca().get_xticks()])
    


def spatialDE(path):
    counts_files = glob.glob(path+  "*[0-2]"+"_count.csv")
    new_counts = np.ravel([[x[:-10] + '_new_count.csv'] for x in counts_files])
    sample_info_new = np.ravel([[x[:-10] + '_new_meta.csv'] for x in counts_files])
    sample_info_files = np.ravel([[x[:-10] + '_meta.csv'] for x in counts_files])
    sample_output = np.ravel([[x[:-10] +'_spatialDE_results' '.csv'] for x in counts_files])
    times_output = np.ravel([[x[:-10] + '_spatialDE_runtime.csv'] for x in counts_files])
    df = pd.DataFrame({'counts':counts_files, 'meta':sample_info_files, 'counts_new':new_counts, 'meta_new':sample_info_new,'output':sample_output, 'time':times_output})
    for index, row in df.iterrows():
        counts = pd.read_csv(row['counts'], index_col=0) # load counts
        counts = counts.T
        counts.to_csv(row['counts_new'],index=False)
        counts = pd.read_csv(row['counts_new'])
        sample_info = pd.read_csv(row['meta'], index_col=0) # load meta data with spatial coordinates
        sample_info.to_csv(row['meta_new'],index=False)
        sample_info = pd.read_csv(row['meta_new'])
        print(counts.shape)
        print(sample_info.shape)
        norm_expr = NaiveDE.stabilize(counts).T # remove tech variation
        resid_expr = NaiveDE.regress_out(sample_info, norm_expr, 'np.log(total)').T
        #X = sample_info[['imagerow', 'imagecol']]
        X = sample_info[['row', 'col']]
        start= datetime.now()
        results = SpatialDE.run(X, resid_expr)
        end = datetime.now()
        time = pd.DataFrame({end - start})
        time.to_csv(row['time']) # Save spatial results
        # total time taken
        print(f"Runtime of the program is {end - start}")
        results.to_csv(row['output']) # Save spatial results
        print("Finished =", row['counts'])
    return


### LIBD (need to change the 'counts_files' path in the function)

path = "./DESpace_data/Simulation/Output/LIBD/Manual_clusters_patch/probs_0.6_0.9FALSE/"
spatialDE(path)

path = "./DESpace_data/Simulation/Output/LIBD/Manual_clusters_patch/probs_0.5_0.9FALSE/"
spatialDE(path)

path = "./DESpace_data/Simulation/Output/LIBD/bottom_patch/probs_0.6_0.9FALSE/"
spatialDE(path)

path = "./DESpace_data/Simulation/Output/LIBD/bottom_patch/probs_0.5_0.9FALSE/"
spatialDE(path)

path = "./DESpace_data/Simulation/Output/LIBD/circle_patch/probs_0.5_0.9FALSE/"
spatialDE(path)

path = "./DESpace_data/Simulation/Output/LIBD/circle_patch/probs_0.6_0.9FALSE/"
spatialDE(path)

path1 = "./DESpace_data/Simulation/Output/LIBD/mixture_patch/probs_0.6_0.6FALSE/"
spatialDE(path1)

path1 = "./DESpace_data/Simulation/Output/LIBD/mixture_reverse_patch/probs_0.6_0.6FALSE/"
spatialDE(path1)

path1 = "./DESpace_data/Simulation/Output/LIBD/mixture_patch/probs_0.5_0.9FALSE/"
spatialDE(path1)

path1 = "./DESpace_data/Simulation/Output/LIBD/mixture_reverse_patch/probs_0.5_0.9FALSE/"
spatialDE(path1)

path = "/home/peicai/project/SV/data/real_data_analysis/LIBD/"
spatialDE(path)
### melanoma

path = "./DESpace_data/Simulation/Output/melanoma/BayesSpace_clusters_patch/probs_0.6_0.9FALSE/"
spatialDE(path)

path = "./DESpace_data/Simulation/Output/melanoma/BayesSpace_clusters_patch/probs_0.5_0.9FALSE/"
spatialDE(path)

path = "./DESpace_data/Simulation/Output/melanoma/right_patch/probs_0.6_0.9FALSE/"
spatialDE(path)

path = "./DESpace_data/Simulation/Output/melanoma/right_patch/probs_0.5_0.9FALSE/"
spatialDE(path)

path = "./DESpace_data/Simulation/Output/melanoma/circle_patch/probs_0.5_0.9FALSE/"
spatialDE(path)

path = "./DESpace_data/Simulation/Output/melanoma/circle_patch/probs_0.6_0.9FALSE/"
spatialDE(path)

path1 = "./DESpace_data/Simulation/Output/melanoma/mixture_patch/probs_0.6_0.6FALSE/"
spatialDE(path1)

path1 = "./DESpace_data/Simulation/Output/melanoma/mixture_reverse_patch/probs_0.6_0.6FALSE/"
spatialDE(path1)

path1 = "./DESpace_data/Simulation/Output/melanoma/mixture_patch/probs_0.5_0.9FALSE/"
spatialDE(path1)

path1 = "./DESpace_data/Simulation/Output/melanoma/mixture_reverse_patch/probs_0.5_0.9FALSE/"
spatialDE(path1)

path = "/home/peicai/project/SV/data/real_data_analysis/melanoma/"
spatialDE(path)

### mouse cerebellum (need to change the 'counts_files' path in the function)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/BayesSpace_clusters_patch/probs_0.6_0.9FALSE/"
spatialDE(path)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/BayesSpace_clusters_patch/probs_0.5_0.9FALSE/"
spatialDE(path)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/right_patch/probs_0.6_0.9FALSE/"
spatialDE(path)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/right_patch/probs_0.5_0.9FALSE/"
spatialDE(path)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/circle_patch/probs_0.5_0.9FALSE/"
spatialDE(path)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/circle_patch/probs_0.6_0.9FALSE/"
spatialDE(path)

path1 = "./DESpace_data/Simulation/Output/mouse cerebellum/mixture_patch/probs_0.6_0.6FALSE/"
spatialDE(path1)

path1 = "./DESpace_data/Simulation/Output/mouse cerebellum/mixture_reverse_patch/probs_0.6_0.6FALSE/"
spatialDE(path1)

path1 = "./DESpace_data/Simulation/Output/mouse cerebellum/mixture_patch/probs_0.5_0.9FALSE/"
spatialDE(path1)

path1 = "./DESpace_data/Simulation/Output/mouse cerebellum/mixture_reverse_patch/probs_0.5_0.9FALSE/"
spatialDE(path1)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/"
# counts_files = glob.glob(path+  "100_count.csv")
spatialDE(path)
