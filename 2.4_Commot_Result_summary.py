# ------------------------------------------------------------------------------
# Title: Summary of Commot Spatial Communication Results
# Author: Yiran Song
# Date: March 18, 2025
# Description:
# This script extracts and summarizes pathway-level ligand-receptor interactions
# from Commot spatial communication results.
#
# Key Functions:
# - Loads spatial communication results from Commot
# - Extracts pathway-specific interactions
# - Saves summary statistics by cell type and niche
# - Normalizes the results by cell type counts
#
# Dependencies:
# - scanpy, anndata, pandas, numpy, commot, plotly
#
# Output:
# - Summary table of interactions
#
# Usage:
# - Run after `2.3_grn_p0.py`
# ------------------------------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import commot as ct
import plotly

plotly_colors = plotly.colors.qualitative.Alphabet


adata = ad.read_h5ad('./cell-contact/adata_p0_spatial_communication_cell_contact.h5ad') # direct communication
adata_2 = ad.read_h5ad('./adata_p0_spatial_communication_v2.h5ad') # indirect communication
# ------------------------------------------------------------------------------
# Summary of cell type-level interactions
# ------------------------------------------------------------------------------

adata.obs['celltype_annotation'] = adata.obs['cytoSPACE.final.anno']
cell_types = adata.obs['celltype_annotation'].cat.categories
results = pd.DataFrame(columns=['Pathway', 'Sender_Celltype', 'Receiver_Celltype', 'Sum_Sender', 'Sum_Receiver'])
pathways = [key for key in adata.obsp.keys() if key.startswith('commot-user_database-')]
pathways

for pathway in pathways:

    comm_matrix = adata.obsp[pathway]
    
    for sender_type in cell_types:
        for receiver_type in cell_types:
            # Get indices for the sender and receiver cell types
            sender_cells = adata.obs.index[adata.obs['celltype_annotation'] == sender_type]
            receiver_cells = adata.obs.index[adata.obs['celltype_annotation'] == receiver_type]
            
            sender_indices = adata.obs_names.isin(sender_cells)
            receiver_indices = adata.obs_names.isin(receiver_cells)
            
            # Calculate the sum of outgoing (rows) and incoming (columns) communication
            sum_sender = comm_matrix[sender_indices, :][:, receiver_indices].sum(axis=1).sum()
            #sum_receiver = comm_matrix.loc[sender_indices, receiver_indices].sum(axis=0).sum()
            sum_receiver = sum_sender
            
            # Append the result to the DataFrame
            results = results.append({
                'Pathway': pathway,
                'Sender_Celltype': sender_type,
                'Receiver_Celltype': receiver_type,
                'Sum_Sender': sum_sender,
                'Sum_Receiver': sum_receiver
            }, ignore_index=True)
            
# Normalize the results by cell type counts
cell_counts = adata.obs['celltype_annotation'].value_counts()
results['Sender_Cell_Count'] = results['Sender_Celltype'].apply(lambda x: cell_counts.get(x, 1))
results['Receiver_Cell_Count'] = results['Receiver_Celltype'].apply(lambda x: cell_counts.get(x, 1))
results['Normalized_Sum_Sender_by_Sender_Count'] = results['Sum_Sender'] / results['Sender_Cell_Count']
results['Normalized_Sum_Receiver_by_Receiver_Count'] = results['Sum_Receiver'] / results['Receiver_Cell_Count']
results['Cell_Product'] = results['Sender_Cell_Count'] * results['Receiver_Cell_Count']
results['Normalized_Sum_Product'] = results['Sum_Sender'] / results['Cell_Product']  # Or 'Sum_Receiver' as needed
results['Normalized_Sum_Product_time107'] = results['Normalized_Sum_Product']  *1e7
results.to_csv("./cell-contact/Commot_results_by_celltype_normalized_p21.csv", index=False)





# ------------------------------------------------------------------------------
# Summary of niche-level interactions (Ligand-Receptor or Pathway)
# ------------------------------------------------------------------------------

results = []

# Process adata
for ligrec_key in adata.obsp.keys():
    if ligrec_key.startswith("commot-user_database-"):
        ligrec_cleaned = ligrec_key.replace("commot-user_database-", "")
        pathway = ligrec_pathway_map.get(ligrec_cleaned, "Unknown")  # Get pathway or "Unknown"

        # Extract the interaction matrix for this ligand-receptor pair
        interaction_matrix = adata.obsp[ligrec_key]

        # Convert the sparse matrix to CSR format if not already in that format
        if not isinstance(interaction_matrix, csr_matrix):
            interaction_matrix = interaction_matrix.tocsr()

        # Compute interaction metrics per niche
        for niche in adata.obs['niches'].unique():
            # Get cells in the current niche
            cells_in_niche = adata.obs[adata.obs['niches'] == niche].index

            # Convert string-based indices to integer positions
            cell_positions = np.array([adata.obs.index.get_loc(cell) for cell in cells_in_niche])

            if len(cell_positions) > 0:
                # Subset the sparse matrix using the integer indices
                interaction_scores = interaction_matrix[cell_positions, :][:, cell_positions].data
                sum_scores = np.sum(interaction_scores) if len(interaction_scores) > 0 else 0
                num_interactions = len(interaction_scores)
                avg_score = sum_scores / len(cell_positions) if num_interactions > 0 else np.nan
                avg_score_divided_by_two = sum_scores / (len(cell_positions) * len(cell_positions)) if num_interactions > 0 else np.nan
            else:
                avg_score = np.nan
                avg_score_divided_by_two = np.nan
                sum_scores = 0
                num_interactions = 0

            # Append the results
            results.append({
                'ligand_receptor': ligrec_cleaned,
                'pathway': pathway,
                'niche': niche,
                'num_cells_in_niche': len(cell_positions),
                'sum_interaction_scores': sum_scores,
                'num_interactions': num_interactions,
                'average_interaction_score_by_one': avg_score,
                'average_interaction_score_by_two': avg_score_divided_by_two
            })

# Repeat the same process for adata_2
for ligrec_key in adata_2.obsp.keys():
    if ligrec_key.startswith("commot-user_database-"):
        ligrec_cleaned = ligrec_key.replace("commot-user_database-", "")
        pathway = ligrec_pathway_map.get(ligrec_cleaned, "Unknown")  # Get pathway or "Unknown"

        # Extract the interaction matrix for this ligand-receptor pair
        interaction_matrix = adata_2.obsp[ligrec_key]

        # Convert the sparse matrix to CSR format if not already in that format
        if not isinstance(interaction_matrix, csr_matrix):
            interaction_matrix = interaction_matrix.tocsr()

        # Compute interaction metrics per niche
        for niche in adata_2.obs['niches'].unique():
            # Get cells in the current niche
            cells_in_niche = adata_2.obs[adata_2.obs['niches'] == niche].index

            # Convert string-based indices to integer positions
            cell_positions = np.array([adata_2.obs.index.get_loc(cell) for cell in cells_in_niche])

            if len(cell_positions) > 0:
                # Subset the sparse matrix using the integer indices
                interaction_scores = interaction_matrix[cell_positions, :][:, cell_positions].data
                sum_scores = np.sum(interaction_scores) if len(interaction_scores) > 0 else 0
                num_interactions = len(interaction_scores)
                avg_score = sum_scores / len(cell_positions) if num_interactions > 0 else np.nan
                avg_score_divided_by_two = sum_scores / (len(cell_positions) * len(cell_positions)) if num_interactions > 0 else np.nan
            else:
                avg_score = np.nan
                avg_score_divided_by_two = np.nan
                sum_scores = 0
                num_interactions = 0

            # Append the results
            results.append({
                'ligand_receptor': ligrec_cleaned,
                'pathway': pathway,
                'niche': niche,
                'num_cells_in_niche': len(cell_positions),
                'sum_interaction_scores': sum_scores,
                'num_interactions': num_interactions,
                'average_interaction_score_by_one': avg_score,
                'average_interaction_score_by_two': avg_score_divided_by_two
            })

df_results = pd.DataFrame(results)
df_results = df_results[df_results['pathway'] != 'Unknown']
df_results.to_csv('average_niche_interaction_scores_p0_v2_sumbynumber.csv', index=False)

