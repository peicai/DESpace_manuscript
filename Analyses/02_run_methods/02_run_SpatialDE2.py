import numpy as np
import scipy
import pandas as pd
import scanpy as sc
import anndata as ad
import SpatialDE
from tqdm.auto import trange, tqdm
import random#, torch
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from IPython.display import set_matplotlib_formats
#%matplotlib inline
#set_matplotlib_formats('png')
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 200
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
from plotnine import *
import mizani
from io import StringIO
from itertools import chain
from bioservices import BioMart, UniProt
import goatools
from datetime import datetime
import glob
from math import sqrt
from joblib import Parallel, delayed

def spatialDE2(path):
  counts_files = glob.glob(path+ '151[5-6][0-7][0-9]' + "_count.csv")
  sample_info_files = np.ravel([[x[:-9] + 'meta.csv'] for x in counts_files])
  sample_output = np.ravel([[x[:-9] + 'spatialDE2_results.csv'] for x in counts_files])
  sample_output_pkl = np.ravel([[x[:-9] + 'spatialDE2_results.pkl'] for x in counts_files])
  loc_file = np.ravel([[x[:-9] + 'location.csv'] for x in counts_files])
  times_output = np.ravel([[x[:-9] + 'spatialDE2_runtime.csv'] for x in counts_files])
  df = pd.DataFrame({'counts':counts_files, 'meta':sample_info_files, 'loc':loc_file, 'output':sample_output, 'pkl':sample_output_pkl, 'time':times_output})
  for index, row in df.iterrows():
      random.seed(123)
      np.random.seed(123)
      counts = pd.read_csv(row['counts'], index_col=0) # load counts
      counts = counts.T[counts.sum(0) >= 0].T  # Filter practically unobserved genes
      sample_info = pd.read_csv(row['meta'], index_col=0) # load meta data with spatial coordinates
      counts = counts.T.loc[sample_info.index].T  # Align count matrix with metadata table
      loc = pd.read_csv(row['loc'], index_col=0)
      loc = loc.iloc[:,0:2]
      loc = loc.rename(columns={'col': 'imagecol', 'row': 'imagerow'})
      X = sample_info[['row', 'col']]
      #print(counts.T.index)
      #print(sample_info.index)
      adata = sc.AnnData(X = counts.T, obs = sample_info, obsm = X)
      print(adata)
      adata.obsm['spatial'] = X
      start= datetime.now()
      svg_full, _ = SpatialDE.test(adata, omnibus=True)
      end = datetime.now()
      # total time taken
      print(f"Runtime of the program is {end - start}")
      svg_full.to_csv(row['output'])
      svg_full["total_counts"] = np.asarray(adata.X.sum(axis=0)).squeeze()
      #svg_full.to_pickle(row['pkl'])
      time = pd.DataFrame({end - start})
      time.to_csv(row['time']) # Save spatial results
      #obsm_data.to_csv(row['output']) # Save spatial results
      print("Finished =", row['counts'])
  return

### LIBD (need to change the 'counts_files' path in the function)

path = "./DESpace_data/Simulation/Output/LIBD/Manual_clusters_patch/probs_0.6_0.9FALSE/"
spatialDE2(path)

path = "./DESpace_data/Simulation/Output/LIBD/Manual_clusters_patch/probs_0.5_0.9FALSE/"
spatialDE2(path)

path = "./DESpace_data/Simulation/Output/LIBD/bottom_patch/probs_0.6_0.9FALSE/"
spatialDE2(path)

path = "./DESpace_data/Simulation/Output/LIBD/bottom_patch/probs_0.5_0.9FALSE/"
spatialDE2(path)

path = "./DESpace_data/Simulation/Output/LIBD/circle_patch/probs_0.5_0.9FALSE/"
spatialDE2(path)

path = "./DESpace_data/Simulation/Output/LIBD/circle_patch/probs_0.6_0.9FALSE/"
spatialDE2(path)

path1 = "./DESpace_data/Simulation/Output/LIBD/mixture_patch/probs_0.6_0.6FALSE/"
spatialDE2(path1)

path1 = "./DESpace_data/Simulation/Output/LIBD/mixture_reverse_patch/probs_0.6_0.6FALSE/"
spatialDE2(path1)

path1 = "./DESpace_data/Simulation/Output/LIBD/mixture_patch/probs_0.5_0.9FALSE/"
spatialDE2(path1)

path1 = "./DESpace_data/Simulation/Output/LIBD/mixture_reverse_patch/probs_0.5_0.9FALSE/"
spatialDE2(path1)

path = "/home/peicai/project/SV/data/real_data_analysis/LIBD/"
spatialDE2(path)
### melanoma

path = "./DESpace_data/Simulation/Output/melanoma/BayesSpace_clusters_patch/probs_0.6_0.9FALSE/"
spatialDE2(path)

path = "./DESpace_data/Simulation/Output/melanoma/BayesSpace_clusters_patch/probs_0.5_0.9FALSE/"
spatialDE2(path)

path = "./DESpace_data/Simulation/Output/melanoma/right_patch/probs_0.6_0.9FALSE/"
spatialDE2(path)

path = "./DESpace_data/Simulation/Output/melanoma/right_patch/probs_0.5_0.9FALSE/"
spatialDE2(path)

path = "./DESpace_data/Simulation/Output/melanoma/circle_patch/probs_0.5_0.9FALSE/"
spatialDE2(path)

path = "./DESpace_data/Simulation/Output/melanoma/circle_patch/probs_0.6_0.9FALSE/"
spatialDE2(path)

path1 = "./DESpace_data/Simulation/Output/melanoma/mixture_patch/probs_0.6_0.6FALSE/"
spatialDE2(path1)

path1 = "./DESpace_data/Simulation/Output/melanoma/mixture_reverse_patch/probs_0.6_0.6FALSE/"
spatialDE2(path1)

path1 = "./DESpace_data/Simulation/Output/melanoma/mixture_patch/probs_0.5_0.9FALSE/"
spatialDE2(path1)

path1 = "./DESpace_data/Simulation/Output/melanoma/mixture_reverse_patch/probs_0.5_0.9FALSE/"
spatialDE2(path1)

path = "/home/peicai/project/SV/data/real_data_analysis/melanoma/"
spatialDE2(path)

### mouse cerebellum (need to change the 'counts_files' path in the function)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/BayesSpace_clusters_patch/probs_0.6_0.9FALSE/"
spatialDE2(path)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/BayesSpace_clusters_patch/probs_0.5_0.9FALSE/"
spatialDE2(path)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/right_patch/probs_0.6_0.9FALSE/"
spatialDE2(path)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/right_patch/probs_0.5_0.9FALSE/"
spatialDE2(path)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/circle_patch/probs_0.5_0.9FALSE/"
spatialDE2(path)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/circle_patch/probs_0.6_0.9FALSE/"
spatialDE2(path)

path1 = "./DESpace_data/Simulation/Output/mouse cerebellum/mixture_patch/probs_0.6_0.6FALSE/"
spatialDE2(path1)

path1 = "./DESpace_data/Simulation/Output/mouse cerebellum/mixture_reverse_patch/probs_0.6_0.6FALSE/"
spatialDE2(path1)

path1 = "./DESpace_data/Simulation/Output/mouse cerebellum/mixture_patch/probs_0.5_0.9FALSE/"
spatialDE2(path1)

path1 = "./DESpace_data/Simulation/Output/mouse cerebellum/mixture_reverse_patch/probs_0.5_0.9FALSE/"
spatialDE2(path1)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/"
spatialDE2(path)
