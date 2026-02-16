# ------------------------------------------------------------------------------
# Title: Gene Regulatory Network Inference (GRNBoost - P0)
# Author: Yiran Song
# Date: March 18, 2025
# Description:
# This script runs GRNBoost2 to infer gene regulatory networks from Commot spatial communication data.
#
# Key Functions:
# - Loads spatial communication results
# - Extracts ligand-receptor interactions
# - Runs GRNBoost2 for regulatory network inference
#
# Dependencies:
# - scanpy, anndata, pandas, numpy, arboreto, distributed
#
# Output:
# - Inferred gene regulatory network (GRN)
#
# Usage:
# - Run after `2.2_commot_p0_indirect.py`
# ------------------------------------------------------------------------------

import pandas as pd
import scanpy as sc
from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names
from distributed import Client, LocalCluster
import numpy as np
import anndata as ad

if __name__ == '__main__':
    
    adata = ad.read_h5ad('./adata_p0_spatial_communication_v2.h5ad')
    adata_2 = ad.read_h5ad('./cell-contact/adata_p0_spatial_communication_cell_contact.h5ad')
    
    receptor_genes = set()
    interaction_keys = [key for key in adata.obsp.keys() if "commot-user_database" in key]
    for key in interaction_keys:
        # Remove the prefix to get the interaction name
        interaction = key.split("commot-user_database-")[1]

        # Check if the interaction is a ligand-receptor pair (denoted by '-')
        if '-' in interaction:
            # Split into ligand and receptors
            parts = interaction.split('-')
            ligand = parts[0]  # The first part is the ligand
            receptors = parts[1].replace('_', ' ').split()  # Receptors can be separated by '_'

            # Add all receptor genes to the set
            receptor_genes.update(receptors)

    receptor_genes = {gene for gene in receptor_genes if gene}
    
    interaction_keys = [key for key in adata_2.obsp.keys() if "commot-user_database" in key]
    for key in interaction_keys:
        # Remove the prefix to get the interaction name
        interaction = key.split("commot-user_database-")[1]

        # Check if the interaction is a ligand-receptor pair (denoted by '-')
        if '-' in interaction:
            # Split into ligand and receptors
            parts = interaction.split('-')
            ligand = parts[0]  # The first part is the ligand
            receptors = parts[1].replace('_', ' ').split()  # Receptors can be separated by '_'

            # Add all receptor genes to the set
            receptor_genes.update(receptors)

    receptor_genes = {gene for gene in receptor_genes if gene}

    adata = ad.read_h5ad("./xenium_cytospace.h5ad")
    time_point_of_interest = 'p0'
    adata = adata[adata.obs['timepoint'] == time_point_of_interest].copy()
    print(adata)
    
    adata_genes = set(adata.var_names)
    # Find intersection of interaction receptor genes and genes in adata
    genes_in_data = receptor_genes.intersection(adata_genes)
    print(genes_in_data)
    gene_names = adata.var_names
    
    if isinstance(adata.X, np.ndarray):
        expression_matrix = adata.X # X[:, gene_indices]
    else:
        expression_matrix = adata.X.toarray()
        
    df_expression = pd.DataFrame(expression_matrix, columns=adata.var_names)
    print(df_expression.head())
    

    
    out_file = './grn/commot_grn_output_p0_cytospace_data.tsv'


    
    # tf_names is read using a utility function included in Arboreto
    #tf_names = load_tf_names(tf_file)
    
    # instantiate a custom Dask distributed Client
    client = Client(LocalCluster())

    # compute the GRN
    network = grnboost2(expression_data=df_expression, tf_names=genes_in_data, client_or_address=client)

    # write the GRN to file
    network.to_csv(out_file, sep='\t', index=False, header=False)