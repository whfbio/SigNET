# ------------------------------------------------------------------------------
# Title: Indirect Spatial Communication Analysis (Commot - P0)
# Author: Yiran Song
# Date: March 18, 2025
# Description:
# This script runs indirect ligand-receptor communication analysis on the Xenium dataset
# at time point P0 using Commot with a larger interaction threshold.
#
# Key Functions:
# - Filters Xenium data for time point P0
# - Loads ligand-receptor interactions from CellChat (mouse)
# - Runs spatial communication with a distance threshold of 200
#
# Dependencies:
# - scanpy, anndata, pandas, numpy, commot
#
# Output:
# - Processed AnnData object with indirect ligand-receptor interactions (H5AD format)
#
# Usage:
# - Run after `2.1_commot_p0_direct.py`
# ------------------------------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import commot as ct

# Load Xenium dataset
adata = ad.read_h5ad('./xenium_total.h5ad')
adata.obsm["spatial"] = adata.obsm["spatial_fov"]

# Select P0 time point
adata_p0 = adata[adata.obs['time_point'] == 'p0'].copy()
del adata.obsm["spatial_fov"]
adata = None  # Free memory

# Load ligand-receptor interactions 
df_ligrec = ct.pp.ligand_receptor_database(database='CellChat', species='mouse')

# Perform spatial communication analysis
ct.tl.spatial_communication(adata_p0, database_name='user_database', df_ligrec=df_ligrec, dis_thr=200, heteromeric=True, pathway_sum=True)

# Save output
adata_p0.write('./adata_p0_spatial_communication_v2.h5ad')
