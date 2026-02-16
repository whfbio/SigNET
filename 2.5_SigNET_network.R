# ------------------------------------------------------------------------------
# Title: Construction of Co-Expression Network (SigNET)
# Author: Yiran Song
# Date: March 18, 2025
# Description:
# This script constructs a co-expression network using transcription factor (TF)
# modules derived from Commot gene regulatory network (GRN) data.
#
# Key Functions:
# - Reads ligand-receptor interaction data
# - Filters interactions based on weight thresholds
# - Constructs TF modules using multiple criteria:
#   - Global weight threshold
#   - Top targets per TF
#   - Top TFs per target
#   - Standard deviation-based filtering
# - Saves the resulting co-expression modules as CSV files
#
# Dependencies:
# - reshape2, dplyr, readr
#
# Output:
# - CSV files containing TF modules for different filtering methods
#
# Usage:
# - Run after GRN inference (`2.3_grn_p0.py`)
# ------------------------------------------------------------------------------

# Load required libraries
library(reshape2)
library(dplyr)
library(readr)

# ------------------------------------------------------------------------------
# Step 1: Set Parameters and Load Data
# ------------------------------------------------------------------------------

# Define threshold parameters
weightThreshold <- 1    # Global weight threshold
topThr <- c(0.1)        # Additional weight thresholds
nTopTfs <- c(5)         # Top TFs per target
nTopTargets <- 50       # Top targets per TF
weightCol <- "weight"   # Name of the weight column

# Load GRN interaction data
linkList <- read_delim("./grn/commot_grn_output_p0_cytospace_data.tsv", delim = "\t", col_names = FALSE)
colnames(linkList) <- c("TF", "Target", weightCol)

# Ensure correct data types
linkList$TF <- as.character(linkList$TF)
linkList$Target <- as.character(linkList$Target)

# ------------------------------------------------------------------------------
# Step 2: Construct TF Modules Based on Weight Thresholds
# ------------------------------------------------------------------------------

# Filter interactions based on weight threshold
linkList <- linkList %>% filter(!!sym(weightCol) >= weightThreshold)

# Initialize list to store TF modules
tfModules <- list()
allName <- paste0("w", format(weightThreshold, scientific=FALSE))
tfModules[[allName]] <- split(linkList$Target, linkList$TF)

# Create modules using additional weight thresholds
if (!is.null(topThr)) {
  topThr <- setNames(topThr, paste0("w", format(topThr, scientific=FALSE)))
  for (i in seq_along(topThr)) {
    llminW <- linkList %>% filter(!!sym(weightCol) > topThr[i])
    tfModules[[names(topThr)[i]]] <- split(llminW$Target, llminW$TF)
  }
}

# Create modules using the top N targets per TF
if (!is.null(nTopTargets)) {
  nTopTargets <- setNames(nTopTargets, paste0("top", nTopTargets))
  for (i in seq_along(nTopTargets)) {
    tfModules[[names(nTopTargets)[i]]] <- lapply(tfModules[[allName]], function(x) x[1:min(length(x), nTopTargets[i])])
  }
}

# ------------------------------------------------------------------------------
# Step 3: Construct TF Modules Based on Top TFs Per Target
# ------------------------------------------------------------------------------

if (!is.null(nTopTfs)) {
  linkList_byTarget <- split(linkList, linkList$Target)
  nTopTfs <- setNames(nTopTfs, paste0("top", nTopTfs, "perTarget"))
  
  topTFsperTarget <- lapply(linkList_byTarget, function(llt) {
    nTFs <- nTopTfs[nTopTfs <= nrow(llt)]
    if (length(nTFs) == 0) return(NULL)
    lapply(nTFs, function(x) llt$TF[1:x])
  })
  
  # Remove NULL entries
  topTFsperTarget <- Filter(Negate(is.null), topTFsperTarget)
  
  # Convert to a data frame
  topTFsperTargetDf <- do.call(rbind, lapply(names(topTFsperTarget), function(target) {
    methods <- names(topTFsperTarget[[target]])
    data.frame(Target = target, TF = unlist(topTFsperTarget[[target]]), method = rep(methods, sapply(topTFsperTarget[[target]], length)))
  }))
}

# ------------------------------------------------------------------------------
# Step 4: Construct Modules Based on Standard Deviation Filtering
# ------------------------------------------------------------------------------

topPerTf <- function(ll, weightCol) {
  ll_split <- split(ll, ll$TF)
  tptf <- lapply(names(ll_split), function(tf) {
    tfMean <- mean(ll_split[[tf]][, weightCol])
    tfSd <- sd(ll_split[[tf]][, weightCol])
    list(
      top1sd = ll_split[[tf]][which(ll_split[[tf]][, weightCol] >= tfMean + tfSd), "Target"],
      top3sd = ll_split[[tf]][which(ll_split[[tf]][, weightCol] >= tfMean + (3 * tfSd)), "Target"]
    )
  })
  names(tptf) <- names(ll_split)
  
  # Convert to data frame
  tptf_melted <- melt(tptf)
  colnames(tptf_melted) <- c("Target", "method", "TF")
  tptf_melted <- tptf_melted[, c("Target", "TF", "method")]
  return(tptf_melted)
}

# Apply the function
byFun <- topPerTf(linkList, weightCol)

# ------------------------------------------------------------------------------
# Step 5: Save Results as CSV Files
# ------------------------------------------------------------------------------

# Convert all TF modules to a single data frame
tfModules_melted <- melt(tfModules)
colnames(tfModules_melted) <- c("Target", "TF", "method")

# Combine all module data frames
tfModulesDf <- bind_rows(tfModules_melted, topTFsperTargetDf, byFun)
tfModulesDf <- unique(tfModulesDf)  # Remove duplicates

# Save full dataset
write_csv(tfModulesDf, "./grn/tfModules_p0_cytospace.csv")

# Save modules based on individual methods
methods <- unique(tfModulesDf$method)

for (method in methods) {
  df_method <- filter(tfModulesDf, method == !!method)
  file_name <- paste0("./grn/tfModules_p0_cytospace_", method, ".csv")
  write_csv(df_method, file_name)
  cat("Saved:", file_name, "\n")
}

