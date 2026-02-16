# SigNET: Spatial Signaling Network Analysis

## Overview

**SigNET** (Signaling-to-transcription regulatory NETwork) is a computational framework for reconstructing ligand-receptor-transcription factor regulatory cascades from spatially resolved transcriptomics data. SigNET identifies co-expression patterns between the receptors and their putative direct/indirect downstream TFs and connects these to the corresponding gene regulatory network inferred by SCENIC. The algorithm integrates spatial cell-cell communication inference with gene regulatory network analysis to identify how extracellular signals are transduced into transcriptional programs within tissue microenvironments.

### Key Features

- **Distance-stratified communication analysis**: Distinguishes direct (juxtacrine) from indirect (paracrine) signaling
- **Receptor-centric GRN inference**: Links cell surface receptors to downstream transcriptional cascades
- **Multi-timepoint compatible**: Applies uniformly across developmental stages or experimental conditions

### Conceptual Framework

```
Spatial Transcriptomics Data
         ↓
    [SigNET Algorithm]
         ↓
    ┌────┴────┬──────────┬─────────────┐
    ↓         ↓          ↓             ↓
 Ligand → Receptor → TF → Target Genes
    ↓         ↓          ↓             ↓
Sender    Receiver   Regulon      Program
 Cell      Cell      Network      Activation
```

**Note**: Throughout this documentation, postnatal day 0 (P0) is used as a representative example. The same workflow applies to all timepoints (P0, P7, P14, P21) or experimental conditions in your dataset.

---

## SigNET Algorithm Components

### Component 1: Distance-Stratified Spatial Communication

**Objective**: Infer ligand-receptor interactions accounting for spatial constraints of different signaling modalities.

**Method**: 
SigNET employs COMMOT (Cang et al., 2023) to model cell-cell communication using optimal transport theory. The algorithm operates at two spatial scales:

- **Direct signaling (≤50μm)**: Captures juxtacrine interactions including cell-cell contact and ECM-receptor binding. Limited to "Cell-Cell Contact" and "ECM-Receptor" categories from the CellChat ligand-receptor database.

- **Indirect signaling (≤200μm)**: Models paracrine and secreted factor signaling. Includes all ligand-receptor categories representing diffusible signals.

For each cell pair (*i*, *j*) and ligand-receptor pair (*L*, *R*), the algorithm computes:
- Spatial proximity weight based on Euclidean distance
- Expression-weighted interaction probability: *P*(*i*→*j*) ∝ *L*<sub>*i*</sub> × *R*<sub>*j*</sub> × *f*(distance)
- Pathway-level aggregation for related ligand-receptor families

**Output**: Sparse sender×receiver matrices stored in `.obsp` of AnnData objects, with keys formatted as `commot-user_database-[Ligand]-[Receptor]`.

---

### Component 2: Receptor-to-Transcription Factor Network Inference

**Objective**: Link extracellular signals (receptors) to intracellular gene regulatory programs (transcription factors).

**Method**:
SigNET extracts all receptor genes participating in identified spatial interactions and treats them as candidate transcriptional regulators. This receptor-centric approach reflects biological reality: receptor activation triggers signaling cascades that culminate in transcription factor activity.

The algorithm applies GRNBoost2 (Moerman et al., 2019), an ensemble gradient boosting method, to infer regulatory edges:

1. **Feature matrix**: Gene expression across all cells at the timepoint of interest
2. **Candidate TFs**: Union of all receptor genes from direct + indirect communication
3. **Target genes**: Genome-wide expression profiles
4. **Edge inference**: XGBoost regression identifies TF→Target relationships with importance weights

**Computational notes**: 
- Uses Dask distributed computing for scalability to thousands of genes
- Deterministic when random seed is fixed
- Outputs weighted edge list (TF, Target, Weight)

**Biological interpretation**: High-weight edges represent putative regulatory relationships where receptor gene expression predicts target gene expression across the cellular landscape.

---

### Component 3: Transcription Factor Module Construction

**Objective**: Organize TF-target relationships into coherent gene programs using multiple filtering strategies.

**Method**:
Raw GRN inference produces thousands of edges with varying confidence. SigNET constructs TF modules using five complementary filtering approaches:

1. **Global weight threshold** (*w* ≥ 1): Retains high-confidence edges across all TFs
2. **Top *N* targets per TF** (default *N*=50): Identifies core regulon for each TF
3. **Top *M* TFs per target** (default *M*=5): Finds primary regulators of each gene
4. **Standard deviation filtering** (μ + σ, μ + 3σ): Identifies outlier edges within each TF's distribution

Each filtering method captures different aspects of regulatory architecture:
- Global threshold: Universally strong relationships
- Top targets: TF-specific core programs
- Top TFs: Gene-specific regulatory logic
- SD-based: TF-specific regulatory specialization

**Output**: Multiple TF module sets, each representing a different perspective on the regulatory network. Consensus edges appearing across methods represent the most robust TF-target relationships.

---

### Component 4: Multi-Scale Communication Summarization

**Objective**: Aggregate single-cell interaction data to biologically interpretable scales (cell types and spatial niches).

**Method**:

**Cell-type level aggregation**:
For each ligand-receptor pair and sender-receiver cell type pair, SigNET computes:
- Raw interaction sum: Σ*S*<sub>*i*→*j*</sub> for all sender cells *i* of type *A* and receiver cells *j* of type *B*
- Normalized metrics accounting for population size:
  - Per-sender normalization: Σ*S* / *N*<sub>sender</sub>
  - Per-receiver normalization: Σ*S* / *N*<sub>receiver</sub>
  - Population product normalization: Σ*S* / (*N*<sub>sender</sub> × *N*<sub>receiver</sub>)

This normalization suite distinguishes genuine signaling preferences from statistical artifacts of cell abundance.

**Niche level aggregation**:
Within spatially-defined niches (identified by CellCharter or similar methods):
- Average interaction score per cell: Mean(*S*<sub>*i*→*j*</sub>) for all cell pairs within niche
- Total active interactions: Count of non-zero *S*<sub>*i*→*j*</sub> values
- Niche-specific pathway enrichment

**Biological interpretation**: Cell-type summaries reveal which cell populations communicate; niche summaries reveal where communication occurs in tissue space.

---

### Component 5: Visualization & Integration

**Objective**: Generate publication-quality visualizations linking spatial organization, communication networks, and gene regulation.

**Visualization types**:
- **Spatial distribution maps**: Overlay interaction scores on tissue coordinates
- **Communication heatmaps**: Cell type × cell type or niche × ligand-receptor matrices
- **Dimensionality reduction**: PCA of niche communication profiles to identify signaling archetypes
- **Network diagrams**: Ligand→Receptor→TF→Target cascade visualization
- **Temporal dynamics**: Comparison across developmental timepoints or conditions

**Key outputs**:
- Niche composition and enrichment analysis
- Pathway-specific communication patterns (e.g., SEMA3, Notch, Kit-Kitl signaling)
- Spatially-resolved interaction maps at single-cell resolution
- Integrated signaling network graphs

---

## Implementation Scripts

The SigNET framework is implemented across six scripts that execute the algorithm components sequentially:

| Script | Component | Input | Output |
|--------|-----------|-------|--------|
| `2.1_commot_direct.py` | Direct communication | Xenium spatial data | Direct interaction matrices (≤50μm) |
| `2.2_commot_indirect.py` | Indirect communication | Xenium spatial data | Indirect interaction matrices (≤200μm) |
| `2.3_grn.py` | GRN inference | Communication results + expression | TF-target edge list |
| `2.4_Commot_Result_summary.py` | Multi-scale aggregation | Interaction matrices | Cell type & niche summaries |
| `2.5_SigNET_network.R` | TF module construction | GRN edge list | Filtered TF modules |
| `2.6_main_figure_generation.R` | Visualization | All previous outputs | Publication figures |

**Timepoint generalization**: These scripts use P0 (postnatal day 0) as the working example. To analyze other timepoints, simply modify the `time_point_of_interest` variable in scripts 2.1-2.4 (e.g., change `'p0'` to `'p7'`, `'p14'`, or `'p21'`). The algorithm logic remains identical across timepoints.

---

## Software Requirements

### Computational Environment
- **RAM**: 32+ GB recommended
- **Python**: ≥3.8
- **R**: ≥4.4.0

#### Python Environment (≥3.8)
```bash
scanpy==1.11.5
pandas==2.2.3
numpy==2.2.4
anndata==0.11.1
commot==0.0.3
arboreto==0.1.6
distributed>=2024.1.0
plotly==5.18.0
scipy>=1.11.0
```

#### R Environment (≥4.4.0)
```r
Seurat (5.1.0)
ggplot2 (3.5.0)
dplyr (1.1.4)
tidyr (1.3.1)
reshape2 (1.4.4)
readr (2.1.5)
patchwork (1.3.0)
gridExtra (2.3)
cowplot (1.1.3)
pheatmap (1.0.12)
ggrepel (0.9.6)
viridis (0.6.5)
tibble (3.2.1)
ComplexHeatmap (Bioconductor)
```

---

## Installation

Complete installation instructions are provided in separate files:
- **Python packages**: `requirements.txt` 
- **R packages**: `install_packages.R`

Quick install:
```bash
# Python
pip install -r requirements.txt

# R
Rscript install_packages.R
```

---

## Usage

Execute scripts sequentially for complete SigNET analysis:

```bash
python 2.1_commot_direct.py
python 2.2_commot_indirect.py
python 2.3_grn.py
python 2.4_Commot_Result_summary.py
Rscript 2.5_SigNET_network.R
Rscript 2.6_main_figure_generation.R
```

To analyze different timepoints, modify the `time_point_of_interest` variable in scripts 2.1-2.4 (e.g., `'p0'` → `'p7'`, `'p14'`, or `'p21'`).

---

## Input Data Format

SigNET requires spatially resolved transcriptomics data in AnnData (`.h5ad`) or Seurat (`.rds`) format with the following annotations:

**Required `.obs` metadata**:
- `time_point`: Experimental timepoint or condition (e.g., 'p0', 'p7', 'p14', 'p21')
- `cytoSPACE.final.anno` or equivalent: Cell type annotations
- `niches`: Spatial niche assignments (from CellCharter or similar)

**Required `.obsm` keys**:
- `spatial` or `spatial_fov`: X,Y coordinates for each cell

**Primary input files**:
- `xenium_total.h5ad` - Integrated spatial transcriptomics dataset
- `xenium_cytospace.h5ad` - Gene expression matrix for GRN inference
- `ref.xenium.with12cellcharter.rds` - Seurat object with niche annotations (for visualization)

---

## Output Files & Interpretation

### Spatial Communication Matrices
- **Format**: AnnData `.obsp` keys: `commot-user_database-[Ligand]-[Receptor]`
- **Structure**: Sparse CSR matrices (sender cells × receiver cells)
- **Values**: Interaction probability scores (0-∞, higher = stronger predicted interaction)

### Gene Regulatory Networks
- **Format**: TSV edge list with columns: TF (receptor gene), Target, Weight
- **Interpretation**: Weight reflects regulatory importance; typically threshold at weight ≥ 1 for high-confidence edges

### TF Modules
- **Format**: Multiple CSV files, one per filtering method
- **Content**: TF-Target pairs organized into coherent gene programs
- **Use**: Pathway enrichment, network visualization, mechanistic hypothesis generation

### Summary Statistics
- **Cell-type summaries**: Normalized sender-receiver communication intensities
- **Niche summaries**: Average interaction scores per cell within spatial domains
- **Use**: Identify cell type communication preferences and niche-specific signaling signatures

---

## Citation

If you use the SigNET framework, please cite:

**This work**:  
Charting Postnatal Heart Development Using In Vivo Single-Cell Functional Genomics
Haofei Wang, Yanhan Dong, Yiran Song, Marazzano Colon, Nicholas Yapundich, Shea Ricketts, Xingyan Liu, Gregory Farber, Yunzhe Qian, Li Qian, Jiandong Liu
bioRxiv 2025.03.10.642473; doi: https://doi.org/10.1101/2025.03.10.642473

**Underlying methods**:
1. **COMMOT**: Cang, Z., & Nie, Q. (2023). Inferring spatial and signaling relationships between cells from single cell transcriptomic data. *Nature Communications*, 14(1), 1-13.

2. **GRNBoost2**: Moerman, T., et al. (2019). GRNBoost2 and Arboreto: efficient and scalable inference of gene regulatory networks. *Bioinformatics*, 35(12), 2159-2161.

3. **Seurat**: Hao, Y., et al. (2021). Integrated analysis of multimodal single-cell data. *Cell*, 184(13), 3573-3587.

4. **Scanpy**: Wolf, F. A., et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19(1), 1-5.

---

## License & Availability

- **Code**: Available at https://github.com/whfbio/SigNET
---

## Contact

For questions regarding the SigNET algorithm:
- **Implementation**: Yiran Song, Haofei Wang
- **Repository**: https://github.com/whfbio/SigNET
---

**Algorithm Version**: 1.0  
**Last Updated**: Feb 2026
