import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
#In order to read in image data, we need to install some package. Here we recommend package "opencv"
#inatll opencv in python
#!pip3 install opencv-python
import cv2
import anndata
import h5py
import datetime
import glob
import time


## Define function
def fitModel2(adata, target, row, ratio):
  # identify SVGs
    x_array=adata.obs['imagerow'].tolist()
    y_array=adata.obs['imagecol'].tolist()
    adata.obs["x_array"] = adata.obs["imagerow"]
    adata.obs["y_array"] = adata.obs["imagecol"]
    raw = adata
    from datetime import datetime
    startTime= datetime.now()
    raw.X=(raw.X.A if issparse(raw.X) else raw.X)
    sc.pp.log1p(raw)
  #Search radius such that each spot in the target domain has approximately 10 neighbors on average
    adj_2d=spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)
    start, end= np.quantile(adj_2d[adj_2d!=0],q=0.001), np.quantile(adj_2d[adj_2d!=0],q=0.1)
    r=spg.search_radius(target_cluster=target, cell_id=adata.obs.index.tolist(), x=x_array, y=y_array, pred=adata.obs["pred"].tolist(), start=start, end=end, num_min=10, num_max=14,  max_run=100)
    #Detect neighboring domains
    nbr_domians=spg.find_neighbor_clusters(target_cluster=target,
                                   cell_id=raw.obs.index.tolist(),
                                   x=raw.obs["x_array"].tolist(),
                                   y=raw.obs["y_array"].tolist(),
                                   pred=raw.obs["pred"].tolist(),
                                   radius=r,
                                   ratio=ratio)
    print(nbr_domians)
    #nbr_domians=nbr_domians[0:3]
    raw.var.index=raw.var['gene_ids']
    de_genes_info=spg.rank_genes_groups(input_adata=raw,
                                target_cluster=target,
                                nbr_list=nbr_domians,
                                label_col="pred",
                                adj_nbr=True,
                                log=True)
    from datetime import datetime
    endTime = datetime.now()
    time = pd.DataFrame({endTime - startTime})
    print(time)
    time.to_csv(row['time']) # Save spatial results
    de_genes_info["target_dmain"]=target
    de_genes_info["time"]=time
    de_genes_info["neighbors"]=str(nbr_domians)
    de_genes_info.to_csv(row['output1'])
    print(target)
    return de_genes_info


def preprocess1(adata, adj, l, n_clusters, x_pixel, y_pixel, row):
  #Set seed
    r_seed=t_seed=n_seed=123
    #Search for suitable resolution
    res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
    ## Run SpaGCN
    clf=spg.SpaGCN()
    clf.set_l(l)
    #Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    #Run
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob=clf.predict()
    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')
    #Do cluster refinement(optional)
    x_array=x_pixel
    y_array=y_pixel
    #shape="hexagon" for Visium data, "square" for ST data.
    adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
    refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
    adata.obs["refined_pred"]=refined_pred
    adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
    adata.obs["x_array"] = adata.obs["imagerow"]
    adata.obs["y_array"] = adata.obs["imagecol"]
    adata.obs["x_pixel"] = adata.obs["imagerow"]
    adata.obs["y_pixel"] = adata.obs["imagecol"]
    plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"]
    #Plot spatial domains
    domains="pred"
    num_celltype=len(adata.obs[domains].unique())
    adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
    ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,size=100000/adata.shape[0])
    plt.savefig(row['plot_pred'], dpi=600)
    plt.close()
    #Plot refined spatial domains
    domains="refined_pred"
    num_celltype=len(adata.obs[domains].unique())
    adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
    ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,size=100000/adata.shape[0])
    plt.savefig(row['plot_refined'], dpi=600)
    plt.close()
    return adata

def input(row):
    counts = pd.read_csv(row['counts'], index_col=0) # load counts
    counts = counts.T[counts.sum(0) >= 0].T  # Filter practically unobserved genes
    sample_info = pd.read_csv(row['meta'], index_col=0) # load meta data with spatial coordinates
    counts = counts.T.loc[sample_info.index].T  # Align count matrix with metadata table
    X = sample_info[['row', 'col']]
    adata = sc.AnnData(X = counts.T, obs = sample_info)
    print(adata)
    sc.pp.log1p(adata)
    adata.obs['imagerow'] = adata.obs['row']
    adata.obs['imagecol'] = adata.obs['col']
    x_pixel=adata.obs['imagerow'].tolist()
    y_pixel=adata.obs['imagecol'].tolist()
    adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
    adata.var = pd.DataFrame(counts.index.values, columns=['gene_ids'])
    ## Expression data preprocessing
    adata.var_names_make_unique()
    p=0.5
    #Find the l value given p
    l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
    return adata, x_pixel, y_pixel, adj, l;



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


def spaGCN(path, n_clusters):
 counts_files = glob.glob(path+ "mel" + '[0-4]' + "_rep1_count.csv")
  sample_info_files = np.ravel([[x[:-9] + 'meta.csv'] for x in counts_files])
  sample_output = np.ravel([[x[:-9] + 'SpaGCN_results.csv'] for x in counts_files])
  sample_output0 = np.ravel([[x[:-9] + 'SpaGCN_results0.csv'] for x in counts_files])
  sample_output1 = np.ravel([[x[:-9] + 'SpaGCN_results1.csv'] for x in counts_files])
  sample_output2 = np.ravel([[x[:-9] + 'SpaGCN_results2.csv'] for x in counts_files])
  sample_output3 = np.ravel([[x[:-9] + 'SpaGCN_results3.csv'] for x in counts_files])
  sample_output4 = np.ravel([[x[:-9] + 'SpaGCN_results4.csv'] for x in counts_files])
  sample_output5 = np.ravel([[x[:-9] + 'SpaGCN_results5.csv'] for x in counts_files])
  sample_output6 = np.ravel([[x[:-9] + 'SpaGCN_results6.csv'] for x in counts_files])
  loc_file = np.ravel([[x[:-9] + 'location.csv'] for x in counts_files])
  times_output = np.ravel([[x[:-9] + 'SpaGCN_runtime.csv'] for x in counts_files])
  plot1 = np.ravel([[x[:-9] + 'SpaGCN_pred.png'] for x in counts_files])
  plot2 = np.ravel([[x[:-9] + 'SpaGCN_refined_pred.png'] for x in counts_files])
  df = pd.DataFrame({'counts':counts_files, 'meta':sample_info_files, 'loc':loc_file, 'output':sample_output, 'output0':sample_output0, 'output1':sample_output1,'output2':sample_output2,'output3':sample_output3,'output4':sample_output4,'output5':sample_output5, 'output6':sample_output6, 'time':times_output, 'plot_pred':plot1, 'plot_refined':plot2})
  np.random.seed(123)

  ratio = 1/10
  for index, row in df.iterrows():
    adata, x_pixel, y_pixel, adj, l = input(row)
    adata = preprocess1(adata, adj, l, n_clusters, x_pixel, y_pixel, row)
    data = []
    for j in range(n_clusters-1):
      target = j
      de_genes_info = fitModel2(adata, target, row, ratio)
      de_genes_info['target'] = target
      data.append(de_genes_info)
      string = "output"
      string += str(j)
      de_genes_info.to_csv(row[string])
    dff2 = pd.concat(data, axis=1)
    dff2['min_pvals_adj'] = dff2[['pvals_adj']].min(axis=1)
    dff2.to_csv(row['output'])
    print("Finished =", row['counts'])
  return
print("Sart running SpaGCN ... ")

path = "./DESpace/Real/melanoma/"
spaGCN(path=path, n_clusters = 3)

## Simulation
path = "./DESpace_data/Simulation/Output/melanoma/BayesSpace_clusters_patch/probs_0.5_0.9FALSE/"
spaGCN(path=path, n_clusters = 2)

path = "./DESpace_data/Simulation/Output/melanoma/right_patch/probs_0.5_0.9FALSE/"
spaGCN(path=path, n_clusters = 2)

path = "./DESpace_data/Simulation/Output/melanoma/circle_patch/probs_0.5_0.9FALSE/"
spaGCN(path=path, n_clusters = 2)

path = "./DESpace_data/Simulation/Output/melanoma/mixture_patch/probs_0.5_0.9FALSE/"
spaGCN(path=path, n_clusters = 5)

path = "./DESpace_data/Simulation/Output/melanoma/mixture_reverse_patch/probs_0.5_0.9FALSE/"
spaGCN(path=path, n_clusters = 5)

path = "./DESpace_data/Simulation/Output/melanoma/BayesSpace_clusters_patch/probs_0.6_0.9FALSE/"
spaGCN(path=path, n_clusters = 2)

print("Finish the first one")
path = "./DESpace_data/Simulation/Output/melanoma/right_patch/probs_0.6_0.9FALSE/"
spaGCN(path=path, n_clusters = 2)

path = "./DESpace_data/Simulation/Output/melanoma/circle_patch/probs_0.6_0.9FALSE/"
spaGCN(path=path, n_clusters = 2)

path = "./DESpace_data/Simulation/Output/melanoma/mixture_patch/probs_0.6_0.6FALSE/"
spaGCN(path=path, n_clusters = 5)

path = "./DESpace_data/Simulation/Output/melanoma/mixture_reverse_patch/probs_0.6_0.6FALSE/"
spaGCN(path=path, n_clusters = 5)


def spaGCN(path, n_clusters):
 # counts_files = glob.glob(path+ "mel" + '[0-4]' + "_rep1_count.csv")
#  counts_files = glob.glob(path+  "*[1-2]"+"_count.csv")
  counts_files = glob.glob(path+ '151[5-6][0-7][0-9]' + "_count.csv")
  sample_info_files = np.ravel([[x[:-9] + 'meta.csv'] for x in counts_files])
  sample_output = np.ravel([[x[:-9] + 'SpaGCN_results.csv'] for x in counts_files])
  sample_output0 = np.ravel([[x[:-9] + 'SpaGCN_results0.csv'] for x in counts_files])
  sample_output1 = np.ravel([[x[:-9] + 'SpaGCN_results1.csv'] for x in counts_files])
  sample_output2 = np.ravel([[x[:-9] + 'SpaGCN_results2.csv'] for x in counts_files])
  sample_output3 = np.ravel([[x[:-9] + 'SpaGCN_results3.csv'] for x in counts_files])
  sample_output4 = np.ravel([[x[:-9] + 'SpaGCN_results4.csv'] for x in counts_files])
  sample_output5 = np.ravel([[x[:-9] + 'SpaGCN_results5.csv'] for x in counts_files])
  sample_output6 = np.ravel([[x[:-9] + 'SpaGCN_results6.csv'] for x in counts_files])
  loc_file = np.ravel([[x[:-9] + 'location.csv'] for x in counts_files])
  times_output = np.ravel([[x[:-9] + 'SpaGCN_runtime.csv'] for x in counts_files])
  plot1 = np.ravel([[x[:-9] + 'SpaGCN_pred.png'] for x in counts_files])
  plot2 = np.ravel([[x[:-9] + 'SpaGCN_refined_pred.png'] for x in counts_files])
  df = pd.DataFrame({'counts':counts_files, 'meta':sample_info_files, 'loc':loc_file, 'output':sample_output, 'output0':sample_output0, 'output1':sample_output1,'output2':sample_output2,'output3':sample_output3,'output4':sample_output4,'output5':sample_output5, 'output6':sample_output6, 'time':times_output, 'plot_pred':plot1, 'plot_refined':plot2})
  np.random.seed(123)

  ratio = 1/10
  for index, row in df.iterrows():
    adata, x_pixel, y_pixel, adj, l = input(row)
    adata = preprocess1(adata, adj, l, n_clusters, x_pixel, y_pixel, row)
    data = []
    for j in range(n_clusters-1):
      target = j
      de_genes_info = fitModel2(adata, target, row, ratio)
      de_genes_info['target'] = target
      data.append(de_genes_info)
      string = "output"
      string += str(j)
      de_genes_info.to_csv(row[string])
    dff2 = pd.concat(data, axis=1)
    dff2['min_pvals_adj'] = dff2[['pvals_adj']].min(axis=1)
    dff2.to_csv(row['output'])
    print("Finished =", row['counts'])
  return
print("Sart running SpaGCN ... ")

path = "./DESpace_data/Real/LIBD/"
spaGCN(path=path, n_clusters = 7)

## Simulation
path = "./DESpace_data/Simulation/Output/LIBD/Manual_clusters_patch/probs_0.5_0.9FALSE/"
spaGCN(path=path, n_clusters = 2)

path = "./DESpace_data/Simulation/Output/LIBD/bottom_patch/probs_0.5_0.9FALSE/"
spaGCN(path=path, n_clusters = 2)

path = "./DESpace_data/Simulation/Output/LIBD/circle_patch/probs_0.5_0.9FALSE/"
spaGCN(path=path, n_clusters = 2)

path = "./DESpace_data/Simulation/Output/LIBD/mixture_patch/probs_0.5_0.9FALSE/"
spaGCN(path=path, n_clusters = 5)

path = "./DESpace_data/Simulation/Output/LIBD/mixture_reverse_patch/probs_0.5_0.9FALSE/"
spaGCN(path=path, n_clusters = 5)

path = "./DESpace_data/Simulation/Output/LIBD/Manual_clusters_patch/probs_0.6_0.9FALSE/"
spaGCN(path=path, n_clusters = 2)

print("Finish the first one")
path = "./DESpace_data/Simulation/Output/LIBD/bottom_patch/probs_0.6_0.9FALSE/"
spaGCN(path=path, n_clusters = 2)

path = "./DESpace_data/Simulation/Output/LIBD/circle_patch/probs_0.6_0.9FALSE/"
spaGCN(path=path, n_clusters = 2)

path = "./DESpace_data/Simulation/Output/LIBD/mixture_patch/probs_0.6_0.6FALSE/"
spaGCN(path=path, n_clusters = 5)

path = "./DESpace_data/Simulation/Output/LIBD/mixture_reverse_patch/probs_0.6_0.6FALSE/"
spaGCN(path=path, n_clusters = 5)


def preprocess1(adata, adj, l, n_clusters, x_pixel, y_pixel, row):
  #Set seed
    r_seed=t_seed=n_seed=123
    #Search for suitable resolution
    res=spg.search_res(adata, adj, l, n_clusters, start=0.01, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
    ## Run SpaGCN
    clf=spg.SpaGCN()
    clf.set_l(l)
    #Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    #Run
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob=clf.predict()
    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')
    #Do cluster refinement(optional)
    x_array=x_pixel
    y_array=y_pixel
    #shape="hexagon" for Visium data, "square" for ST data.
    adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
    refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
    adata.obs["refined_pred"]=refined_pred
    adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
    adata.obs["x_array"] = adata.obs["imagerow"]
    adata.obs["y_array"] = adata.obs["imagecol"]
    adata.obs["x_pixel"] = adata.obs["imagerow"]
    adata.obs["y_pixel"] = adata.obs["imagecol"]
    plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"]
    #Plot spatial domains
    domains="pred"
    num_celltype=len(adata.obs[domains].unique())
    adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
    ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,size=100000/adata.shape[0])
    plt.savefig(row['plot_pred'], dpi=600)
    plt.close()
    #Plot refined spatial domains
    domains="refined_pred"
    num_celltype=len(adata.obs[domains].unique())
    adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
    ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,size=100000/adata.shape[0])
    plt.savefig(row['plot_refined'], dpi=600)
    plt.close()
    return adata

def spaGCN(path, n_clusters):
  counts_files = glob.glob(path+  "100_count.csv")
  sample_info_files = np.ravel([[x[:-9] + 'meta.csv'] for x in counts_files])
  sample_output = np.ravel([[x[:-9] + 'SpaGCN_results.csv'] for x in counts_files])
  sample_output0 = np.ravel([[x[:-9] + 'SpaGCN_results0.csv'] for x in counts_files])
  sample_output1 = np.ravel([[x[:-9] + 'SpaGCN_results1.csv'] for x in counts_files])
  sample_output2 = np.ravel([[x[:-9] + 'SpaGCN_results2.csv'] for x in counts_files])
  sample_output3 = np.ravel([[x[:-9] + 'SpaGCN_results3.csv'] for x in counts_files])
  sample_output4 = np.ravel([[x[:-9] + 'SpaGCN_results4.csv'] for x in counts_files])
  loc_file = np.ravel([[x[:-9] + 'location.csv'] for x in counts_files])
  times_output = np.ravel([[x[:-9] + 'SpaGCN_runtime.csv'] for x in counts_files])
  plot1 = np.ravel([[x[:-9] + 'SpaGCN_pred.png'] for x in counts_files])
  plot2 = np.ravel([[x[:-9] + 'SpaGCN_refined_pred.png'] for x in counts_files])
  df = pd.DataFrame({'counts':counts_files, 'meta':sample_info_files, 'loc':loc_file, 'output':sample_output, 'output0':sample_output0, 'output1':sample_output1,'output2':sample_output2,'output3':sample_output3,'output4':sample_output4,'time':times_output, 'plot_pred':plot1, 'plot_refined':plot2})
  np.random.seed(123)

  ratio = 1/4
  for index, row in df.iterrows():
    adata, x_pixel, y_pixel, adj, l = input(row)
    adata = preprocess1(adata, adj, l, n_clusters, x_pixel, y_pixel, row)
    data = []
    for j in range(n_clusters-1):
      target = j
      de_genes_info = fitModel2(adata, target, row, ratio)
      de_genes_info['target'] = target
      data.append(de_genes_info)
      string = "output"
      string += str(j)
      de_genes_info.to_csv(row[string])
      print("Finish ", row[string])
    dff2 = pd.concat(data, axis=1)
    dff2['min_pvals_adj'] = dff2[['pvals_adj']].min(axis=1)
    dff2.to_csv(row['output'])
    print("Finished =", row['counts'])
  return
print("Sart running SpaGCN ... ")

path = "./DESpace_data/Real/mouse_cerebellum/"
spaGCN(path=path, n_clusters = 4)
## Simulation
### mouse cerebellum (need to change the 'counts_files' path in the function)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/stLearn_clusters_patch/probs_0.6_0.9FALSE/"
spaGCN(path=path, n_clusters = 4)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/stLearn_clusters/probs_0.5_0.9FALSE/"
spaGCN(path=path, n_clusters = 4)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/right_patch/probs_0.6_0.9FALSE/"
spaGCN(path=path, n_clusters = 2)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/right_patch/probs_0.5_0.9FALSE/"
spaGCN(path=path, n_clusters = 2)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/circle_patch/probs_0.5_0.9FALSE/"
spaGCN(path=path, n_clusters = 2)

path = "./DESpace_data/Simulation/Output/mouse cerebellum/circle_patch/probs_0.6_0.9FALSE/"
spaGCN(path=path, n_clusters = 2)

path1 = "./DESpace_data/Simulation/Output/mouse cerebellum/mixture_patch/probs_0.6_0.6FALSE/"
spaGCN(path=path, n_clusters = 4)

path1 = "./DESpace_data/Simulation/Output/mouse cerebellum/mixture_reverse_patch/probs_0.6_0.6FALSE/"
spaGCN(path=path, n_clusters = 4)

path1 = "./DESpace_data/Simulation/Output/mouse cerebellum/mixture_patch/probs_0.5_0.9FALSE/"
spaGCN(path=path, n_clusters = 4)

path1 = "./DESpace_data/Simulation/Output/mouse cerebellum/mixture_reverse_patch/probs_0.5_0.9FALSE/"
spaGCN(path=path, n_clusters = 4)
