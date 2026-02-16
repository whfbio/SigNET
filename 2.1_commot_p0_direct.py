# ------------------------------------------------------------------------------
# Title: Direct Spatial Communication Analysis (Commot - P0)
# Author: Yiran Song
# Date: March 18, 2025
# Description:
# This script runs direct ligand-receptor communication analysis on the Xenium dataset
# at time point P0 using the Commot package.
#
# Key Functions:
# - Filters Xenium data for time point P0
# - Loads ligand-receptor interactions from CellChat (mouse)
# - Runs spatial communication with a distance threshold of 50
#
# Dependencies:
# - scanpy, anndata, pandas, numpy, commot
#
# Output:
# - Processed AnnData object with ligand-receptor interactions (H5AD format)
#
# Usage:
# - Run after processing Xenium data with SCENIC integration.
# ------------------------------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import commot as ct

# Load the Xenium data
adata = ad.read_h5ad('./xenium_total.h5ad')
adata.obsm["spatial"] = adata.obsm["spatial_fov"]

# Filter for time point P0
time_point_of_interest = 'p0'
adata_p0 = adata[adata.obs['time_point'] == time_point_of_interest].copy()

del adata.obsm["spatial_fov"]
adata = None  

# Load ligand-receptor interactions (ECM-Receptor & Cell-Cell Contact)
df_ligrec=ct.pp.ligand_receptor_database(database='CellChat', species='mouse',signaling_type =None)
df_ligrec = df_ligrec[df_ligrec.iloc[:,3].isin(['ECM-Receptor','Cell-Cell Contact'])]

# Run spatial communication analysis
ct.tl.spatial_communication(adata_p0, database_name='user_database', df_ligrec=df_ligrec, dis_thr=50, heteromeric=True, pathway_sum=True)

# Save the processed AnnData object
adata_p0.write('./adata_p0_spatial_communication_cell_contact.h5ad')