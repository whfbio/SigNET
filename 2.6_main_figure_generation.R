# ------------------------------------------------------------------------------
# Title: Figure 2 - Spatial Niche Composition and Cellchat Analysis
# Author: Haofei Wang
# Date: March 24, 2025
# Description:
# This script generates visualizations for spatial clustering results across timepoints 
# using Seurat and CellCharter annotations. Includes:
# - Spatial niche UMAPs
# - Cell type vs. niche dot plots
# - Barplots of niche proportions
# - Cell-cell communication summary (Commot)
#
# Dependencies:
# - Seurat
# - ggplot2
# - patchwork
# - reshape2
#
# Usage:
# Load `ref.xenium` object before running.
# ------------------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(gridExtra)
library(cowplot)
library(grid)
library(pheatmap)
library(ggrepel)
wd<-"./PS_38/analysis/"
setwd(wd)
source("./PS_16/Reference_data_analysis/fig_colors.R")
source("./PS_38/analysis/scripts/cellchat_subset_functions.R")
source("./PS_38/analysis/scripts/functions.R")

###Figure 2A
ref.xenium<-readRDS("./robjs/ref.xenium.with12cellcharter.rds")
fov.list<-list("p0"="fov","p7"="fov.2","p14"="fov.3","p21"="fov.4")
p<- list()
timepoint <-"p0"
xenium.temp <- subset(ref.xenium,time_point==timepoint)
p[[timepoint]]<- ImageDimPlot(xenium.temp, fov = fov.list[[timepoint]], nmols = 20000, group.by = "X12_cluster_cellcharter", cols = xenium_12niche_color,
                              boundaries = "segmentation", border.color = NA, dark.background = F, crop = T, axes = F) +
  ggtitle(paste0(timepoint,"_all"))
print(paste0("finish plotting source ", timepoint))
timepoint <-"p7"
xenium.temp <- subset(ref.xenium,time_point==timepoint)
p[[timepoint]]<- ImageDimPlot(xenium.temp, fov = fov.list[[timepoint]], nmols = 20000, group.by = "X12_cluster_cellcharter", cols = xenium_12niche_color,
                              boundaries = "segmentation", border.color = NA, dark.background = F, crop = T, axes = F) +
  ggtitle(paste0(timepoint,"_all"))
print(paste0("finish plotting source ", timepoint))
timepoint <-"p14"
xenium.temp <- subset(ref.xenium,time_point==timepoint)
p[[timepoint]]<- ImageDimPlot(xenium.temp, fov = fov.list[[timepoint]], nmols = 20000, group.by = "X12_cluster_cellcharter", cols = xenium_12niche_color,
                              boundaries = "segmentation", border.color = NA, dark.background = F, crop = T, axes = F) +
  ggtitle(paste0(timepoint,"_all"))
print(paste0("finish plotting source ", timepoint))
timepoint <-"p21"
xenium.temp <- subset(ref.xenium,time_point==timepoint)
p[[timepoint]]<- ImageDimPlot(xenium.temp, fov = fov.list[[timepoint]], nmols = 20000, group.by = "X12_cluster_cellcharter", cols = xenium_12niche_color,
                              boundaries = "segmentation", border.color = NA, dark.background = F, crop = T, axes = F) +
  ggtitle(paste0(timepoint,"_all"))
print(paste0("finish plotting source ", timepoint))

pdf("./plots/cellcharter_niche_plots_12_recolor.pdf", width = 10, height= 10)
grid.arrange(rotate_image(p[["p0"]],90), rotate_image(p[["p7"]],180), rotate_image(p[["p14"]],90), rotate_image(p[["p21"]],90),nrow =2, ncol = 2)
dev.off()

###Figure 2B
xenium.celltype<-c(xenium_cm_lv2_pal,xenium_ec_lv2_pal,xenium_fib_lv2_pal,xenium_peri_lv2_pal,xenium_other_lv2_pal) 
xenium.celltype.name<-names(xenium.celltype)
timepoint <-"p14"
##examine the count matrix and remove the one showing limited cells numbers
xenium.temp <- subset(ref.xenium,time_point==timepoint)
#data<-as.data.frame(table(xenium.temp$cytoSPACE.final.anno, xenium.temp$X18_cluster_cellcharter))
data<-as.data.frame(table(ref.xenium$cytoSPACE.final.anno, ref.xenium$X12_cluster_cellcharter))

colnames(data)<-c("CellType","Niche","Count")

# Calculate the total counts for each column
totals <- data %>%
  group_by(Niche) %>%
  summarise(Total = sum(Count))

# Merge the totals back with the original data
data <- data %>%
  left_join(totals, by = "Niche") %>%
  mutate(Percentage = (Count / Total) * 100)
#data<-filter(data, data$Niche != niche_remove)
library(tibble)  # Required for rownames_to_column and column_to_rownames
count_matrix <- data %>%
  select(CellType, Niche, Count) %>% # Keep only necessary columns
  pivot_wider(names_from = Niche, values_from = Count) %>% # Convert to wide format
  column_to_rownames("CellType") # Set CellType as row names
count_matrix <- as.matrix(count_matrix)

fc.matrix<-count_matrix
for (i in 1:nrow(count_matrix)){
  for (j in 1:ncol(count_matrix)){
    precentage.1<-count_matrix[i,j]/sum(count_matrix[,j])
    precentage.2<-sum(count_matrix[i,-j])/sum(count_matrix[,-j])
    fc<-precentage.1/precentage.2
    fc.matrix[i,j]<-fc
  }
}
library(reshape2)
# Convert the fold change matrix to long format
fold_change_data <- melt(fc.matrix)
colnames(fold_change_data) <- c("CellType", "Niche", "FoldChange")

# Convert Niche in both data frames to character
data$Niche <- as.character(data$Niche)
fold_change_data$Niche <- as.character(fold_change_data$Niche)

# Merge the fold change data with the original data
data <- data %>%
  left_join(fold_change_data, by = c("CellType", "Niche"))

# Define the pseudocount
pseudocount <- 1e-6

# Add the pseudocount to FoldChange and calculate the log2 fold change
data <- data %>%
  mutate(FoldChange_Adjusted = FoldChange + pseudocount,
         Log2FoldChange = log2(FoldChange_Adjusted))

#cap the fold change at 5
data[, "Log2FoldChange"] <- pmin(data[, "Log2FoldChange"], 2)
data[, "Log2FoldChange"] <- pmax(data[, "Log2FoldChange"], -2)


# Apply the factor with the desired order
data$CellType <- factor(data$CellType, levels = rev(xenium.celltype.name))
data$Niche <- factor(data$Niche, levels = c("11","10","9","8","7","6","5","4","3","2","1","0"))

# Create dot plot
# Create the dot plot with color mapped to fold change
pdf(paste0("./plots/total_cellcharter_12niche_cellcomposition.pdf"), height=5, width =10)
ggplot(data, aes(y = Niche, x = CellType)) +
  geom_point(aes(size = Percentage, color = Log2FoldChange)) +
  scale_size_continuous(range = c(0, 8)) +  # Adjust range for better visibility
  scale_color_gradient2(low = "#3c5b54", mid="#efe6d5", high = "#e73213", midpoint=0) +  # Adjust color scheme as needed
  theme_minimal() +
  labs(title = paste0("Proportion of Cell Types in Different Niches"),
       x = "Niche",
       y = "Cell Type",
       size = "Proportion",
       color = "log2 Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines
dev.off()

###Figure 2C
Idents(ref.xenium)<-"time_point"
pt <- table(Idents(ref.xenium), ref.xenium$cluster_cellcharter)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
pt$Var1<- factor(pt$Var1, levels = c("p0","p7","p14","p21"))

pdf("./plots/cellcharter_7niche_cellcomposition_barplot.pdf", height=6, width =4)
ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = xenium_niche_color_p0_p7) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines
dev.off()

Idents(ref.xenium)<-"time_point"
pt <- table(Idents(ref.xenium), ref.xenium$X12_cluster_cellcharter)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
pt$Var1<- factor(pt$Var1, levels = c("p0","p7","p14","p21"))
pt$Var2 <- factor(pt$Var2, levels = c("0", "1", "2","3", "4", "5","6", "7", "8","9", "10", "11"))

pdf("./plots/cellcharter_12niche_cellcomposition_barplot.pdf", height=6, width =4)
ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = xenium_12niche_color) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines
dev.off()

# Create the bar plot
pdf(paste0("./plots/", timepoint, "_niche_cellcomposition.pdf"))
ggplot(data, aes(x = Niche, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = CellType), position = position_stack(vjust = 0.5), size =3) +
  theme_minimal() +
  labs(title = "Percentage of Each Cell Type in Each Niche",
       x = "Niche",
       y = "Percentage",
       fill = "Cell Type")
dev.off()

###Figure 2E
result.direct.total<-data.frame(celltype = as.character(), sender_signal_intensity=as.numeric(),receiver_signal_intensity=as.numeric(), time_point =as.character())
result.secreted.total<-data.frame(celltype = as.character(), sender_signal_intensity=as.numeric(),receiver_signal_intensity=as.numeric(), time_point =as.character())

p0.direct<-read.csv("./COMMOT/results/Direct_contact/Commot_results_by_celltype_normalized_p0.csv")
p0.secret<-read.csv("./COMMOT/results/secreted/Commot_results_by_celltype_normalized_p0.csv")
p0.celltype<-c("Ankrd1.CM","Ddc.CM","Epic","Peri","Postn.Fib","SMC","Slit2.CM",
               "a.CM","arteriole.EC","cap.EC","cap.EC.High.Myo.","endocardial.EC",
               "p.CM","p.EC","p.Fib","v.CM","valve.Fib")

p0.direct<-filter(p0.direct,p0.direct$Sender_Celltype %in% p0.celltype & p0.direct$Receiver_Celltype %in% p0.celltype)
p0.direct <- p0.direct %>%
  mutate(pathway_name = gsub("commot-user_database-", "", Pathway))
p0.direct.pathway<- filter(p0.direct, pathway_name %in% unique(LR_database_mouse$pathway_name))
p0.secret<-filter(p0.secret,p0.secret$Sender_Celltype %in% p0.celltype & p0.secret$Receiver_Celltype %in% p0.celltype)
p0.secret <- p0.secret %>%
  mutate(pathway_name = gsub("commot-user_database-", "", Pathway))
p0.secret.pathway<- filter(p0.secret, pathway_name %in% unique(LR_database_mouse$pathway_name))

result.direct<-data.frame(celltype = as.character(), sender_signal_intensity=as.numeric(),receiver_signal_intensity=as.numeric())
i<-1
for (celltype in p0.celltype){
  df<-filter(p0.direct.pathway, p0.direct.pathway$Sender_Celltype == celltype)
  sender.intensity<-sum(df$Normalized_Sum_Receiver_by_Receiver_Count)
  df<-filter(p0.direct.pathway, p0.direct.pathway$Receiver_Celltype == celltype)
  receiver.intensity<-sum(df$Normalized_Sum_Receiver_by_Receiver_Count)
  result.direct[i,]=c(celltype, sender.intensity, receiver.intensity)
  i <- i+1
}

result.secreted<-data.frame(celltype = as.character(), sender_signal_intensity=as.numeric(),receiver_signal_intensity=as.numeric())
i<-1
for (celltype in p0.celltype){
  df<-filter(p0.secret.pathway, p0.secret.pathway$Sender_Celltype == celltype)
  sender.intensity<-sum(df$Normalized_Sum_Receiver_by_Receiver_Count)
  df<-filter(p0.secret.pathway, p0.secret.pathway$Receiver_Celltype == celltype)
  receiver.intensity<-sum(df$Normalized_Sum_Receiver_by_Receiver_Count)
  result.secreted[i,]=c(celltype, sender.intensity, receiver.intensity)
  i <- i+1
}
result.direct$time_point<-"p0"
result.secreted$time_point<-"p0"
result.direct.total<-rbind(result.direct.total,result.direct)
result.secreted.total<-rbind(result.secreted.total,result.secreted)

p7.direct<-read.csv("./COMMOT/results/Direct_contact/Commot_results_by_celltype_normalized_p7.csv")
p7.secret<-read.csv("./COMMOT/results/secreted/Commot_results_by_celltype_normalized_p7.csv")
p7.celltype<-c("Ankrd1.CM","Ddc.CM","Epic","Peri","Postn.Fib","SMC","Slit2.CM",
               "a.CM","arteriole.EC","cap.EC","cap.EC.High.Myo.","endocardial.EC",
               "p.CM","p.EC","p.Fib","v.CM","valve.Fib","Gfpt2.Fib")
p7.direct<-filter(p7.direct,p7.direct$Sender_Celltype %in% p7.celltype & p7.direct$Receiver_Celltype %in% p7.celltype)
p7.direct <- p7.direct %>%
  mutate(pathway_name = gsub("commot-user_database-", "", Pathway))
p7.direct.pathway<- filter(p7.direct, pathway_name %in% unique(LR_database_mouse$pathway_name))
p7.secret<-filter(p7.secret,p7.secret$Sender_Celltype %in% p7.celltype & p7.secret$Receiver_Celltype %in% p7.celltype)
p7.secret <- p7.secret %>%
  mutate(pathway_name = gsub("commot-user_database-", "", Pathway))
p7.secret.pathway<- filter(p7.secret, pathway_name %in% unique(LR_database_mouse$pathway_name))

result.direct<-data.frame(celltype = as.character(), sender_signal_intensity=as.numeric(),receiver_signal_intensity=as.numeric())
i<-1
for (celltype in p7.celltype){
  df<-filter(p7.direct.pathway, p7.direct.pathway$Sender_Celltype == celltype)
  sender.intensity<-sum(df$Normalized_Sum_Receiver_by_Receiver_Count)
  df<-filter(p7.direct.pathway, p7.direct.pathway$Receiver_Celltype == celltype)
  receiver.intensity<-sum(df$Normalized_Sum_Receiver_by_Receiver_Count)
  result.direct[i,]=c(celltype, sender.intensity, receiver.intensity)
  i <- i+1
}
result.secreted<-data.frame(celltype = as.character(), sender_signal_intensity=as.numeric(),receiver_signal_intensity=as.numeric())
i<-1
for (celltype in p7.celltype){
  df<-filter(p7.secret.pathway, p7.secret.pathway$Sender_Celltype == celltype)
  sender.intensity<-sum(df$Normalized_Sum_Receiver_by_Receiver_Count)
  df<-filter(p7.secret.pathway, p7.secret.pathway$Receiver_Celltype == celltype)
  receiver.intensity<-sum(df$Normalized_Sum_Receiver_by_Receiver_Count)
  result.secreted[i,]=c(celltype, sender.intensity, receiver.intensity)
  i <- i+1
}
result.direct$time_point<-"p7"
result.secreted$time_point<-"p7"
result.direct.total<-rbind(result.direct.total,result.direct)
result.secreted.total<-rbind(result.secreted.total,result.secreted)

p14.direct<-read.csv("./COMMOT/results/Direct_contact/Commot_results_by_celltype_normalized_p14.csv")
p14.secret<-read.csv("./COMMOT/results/secreted/Commot_results_by_celltype_normalized_p14.csv")
p14.celltype<-c("Ankrd1.CM","Epic","Peri","Postn.Fib","SMC","Slit2.CM",
                "a.CM","arteriole.EC","cap.EC","cap.EC.High.Myo.","endocardial.EC",
                "v.CM","valve.Fib","Gfpt2.Fib","Col8a1.Fib","Lsamp.Fib","p.CM","p.EC","p.Fib")
p14.direct<-filter(p14.direct,p14.direct$Sender_Celltype %in% p14.celltype & p14.direct$Receiver_Celltype %in% p14.celltype)
p14.direct <- p14.direct %>%
  mutate(pathway_name = gsub("commot-user_database-", "", Pathway))
p14.direct.pathway<- filter(p14.direct, pathway_name %in% unique(LR_database_mouse$pathway_name))
p14.secret<-filter(p14.secret,p14.secret$Sender_Celltype %in% p14.celltype & p14.secret$Receiver_Celltype %in% p14.celltype)
p14.secret <- p14.secret %>%
  mutate(pathway_name = gsub("commot-user_database-", "", Pathway))
p14.secret.pathway<- filter(p14.secret, pathway_name %in% unique(LR_database_mouse$pathway_name))

result.direct<-data.frame(celltype = as.character(), sender_signal_intensity=as.numeric(),receiver_signal_intensity=as.numeric())
i<-1
for (celltype in p14.celltype){
  df<-filter(p14.direct.pathway, p14.direct.pathway$Sender_Celltype == celltype)
  sender.intensity<-sum(df$Normalized_Sum_Receiver_by_Receiver_Count)
  df<-filter(p14.direct.pathway, p14.direct.pathway$Receiver_Celltype == celltype)
  receiver.intensity<-sum(df$Normalized_Sum_Receiver_by_Receiver_Count)
  result.direct[i,]=c(celltype, sender.intensity, receiver.intensity)
  i <- i+1
}
result.secreted<-data.frame(celltype = as.character(), sender_signal_intensity=as.numeric(),receiver_signal_intensity=as.numeric())
i<-1
for (celltype in p14.celltype){
  df<-filter(p14.secret.pathway, p14.secret.pathway$Sender_Celltype == celltype)
  sender.intensity<-sum(df$Normalized_Sum_Receiver_by_Receiver_Count)
  df<-filter(p14.secret.pathway, p14.secret.pathway$Receiver_Celltype == celltype)
  receiver.intensity<-sum(df$Normalized_Sum_Receiver_by_Receiver_Count)
  result.secreted[i,]=c(celltype, sender.intensity, receiver.intensity)
  i <- i+1
}
result.direct$time_point<-"p14"
result.secreted$time_point<-"p14"
result.direct.total<-rbind(result.direct.total,result.direct)
result.secreted.total<-rbind(result.secreted.total,result.secreted)

p21.direct<-read.csv("./COMMOT/results/Direct_contact/Commot_results_by_celltype_normalized_p21.csv")
p21.secret<-read.csv("./COMMOT/results/secreted/Commot_results_by_celltype_normalized_p21.csv")
p21.celltype<-c("Ankrd1.CM","Epic","Peri","Postn.Fib","SMC","Slit2.CM",
                "a.CM","arteriole.EC","cap.EC","cap.EC.High.Myo.","endocardial.EC",
                "v.CM","valve.Fib","Gfpt2.Fib","Col8a1.Fib","Lsamp.Fib")
p21.direct<-filter(p21.direct,p21.direct$Sender_Celltype %in% p21.celltype & p21.direct$Receiver_Celltype %in% p21.celltype)
p21.direct <- p21.direct %>%
  mutate(pathway_name = gsub("commot-user_database-", "", Pathway))
p21.direct.pathway<- filter(p21.direct, pathway_name %in% unique(LR_database_mouse$pathway_name))
p21.secret<-filter(p21.secret,p21.secret$Sender_Celltype %in% p21.celltype & p21.secret$Receiver_Celltype %in% p21.celltype)
p21.secret <- p21.secret %>%
  mutate(pathway_name = gsub("commot-user_database-", "", Pathway))
p21.secret.pathway<- filter(p21.secret, pathway_name %in% unique(LR_database_mouse$pathway_name))

result.direct<-data.frame(celltype = as.character(), sender_signal_intensity=as.numeric(),receiver_signal_intensity=as.numeric())
i<-1
for (celltype in p21.celltype){
  df<-filter(p21.direct.pathway, p21.direct.pathway$Sender_Celltype == celltype)
  sender.intensity<-sum(df$Normalized_Sum_Receiver_by_Receiver_Count)
  df<-filter(p21.direct.pathway, p21.direct.pathway$Receiver_Celltype == celltype)
  receiver.intensity<-sum(df$Normalized_Sum_Receiver_by_Receiver_Count)
  result.direct[i,]=c(celltype, sender.intensity, receiver.intensity)
  i <- i+1
}
result.secreted<-data.frame(celltype = as.character(), sender_signal_intensity=as.numeric(),receiver_signal_intensity=as.numeric())
i<-1
for (celltype in p21.celltype){
  df<-filter(p21.secret.pathway, p21.secret.pathway$Sender_Celltype == celltype)
  sender.intensity<-sum(df$Normalized_Sum_Receiver_by_Receiver_Count)
  df<-filter(p21.secret.pathway, p21.secret.pathway$Receiver_Celltype == celltype)
  receiver.intensity<-sum(df$Normalized_Sum_Receiver_by_Receiver_Count)
  result.secreted[i,]=c(celltype, sender.intensity, receiver.intensity)
  i <- i+1
}

result.direct$time_point<-"p21"
result.secreted$time_point<-"p21"
result.direct.total<-rbind(result.direct.total,result.direct)
result.secreted.total<-rbind(result.secreted.total,result.secreted)


xenium_other_lv2_pal <- c("Blood"="#f4895f", "Epic"="#95cf92", "SMC"="#ffd92f", "Neur"="#AAF6D2")
xenium_cm_lv2_pal<-c("Ddc.CM"="#651251", "v.CM"= '#1E0012', "Ankrd1.CM"='#D29ACA', "p.CM" = "#CFC8DA","Slit2.CM"="#B9A2BE" ,"a.CM" ="#9F7498","av.CM"="#72476D")
xenium_ec_lv2_pal<-c("cap.EC" ="#183642","arteriole.EC"="#6DB2BF","endocardial.EC"="#4699B7", "p.EC"="#CEE9F4","lym.EC"="#94CFC9","cap.EC.High.Myo."="#20566E")
xenium_fib_lv2_pal<-c("Gfpt2.Fib"="#cb4154","Col8a1.Fib"="#d55d6c","Lsamp.Fib"="#e07984","Postn.Fib"="#ea949b", "valve.Fib"="#f5b0b3", "p.Fib"="#ffcccb")
xenium_peri_lv2_pal<-c("p.Peri"="#b58188","Peri"="#6f1926")
color.scheme<-c(xenium_peri_lv2_pal,xenium_fib_lv2_pal,xenium_ec_lv2_pal,xenium_cm_lv2_pal,xenium_other_lv2_pal)
ordered_celltypes <- c("v.CM", "Ddc.CM", "av.CM", "a.CM", "Ankrd1.CM", "Slit2.CM", 
                       "p.CM", "cap.EC", "cap.EC.High.Myo.", "endocardial.EC", "arteriole.EC", 
                       "lym.EC", "p.EC", "Gfpt2.Fib", "Col8a1.Fib", "Lsamp.Fib", "Postn.Fib", "valve.Fib","p.Fib",
                       "Peri","Epic","SMC")

# Ensure the celltype column is a factor with the desired order
result.direct.total$celltype <- factor(result.direct.total$celltype, levels = ordered_celltypes)
result.direct.total$time_point <- factor(result.direct.total$time_point, levels = c("p0", "p7", "p14", "p21"))
result.direct.total$sender_signal_intensity <- as.numeric(result.direct.total$sender_signal_intensity)
result.direct.total$receiver_signal_intensity <- as.numeric(result.direct.total$receiver_signal_intensity)
pdf("./COMMOT/plots/direct_signal_intensity_over_time.pdf", onefile = T)
ggplot(result.direct.total, aes(x = time_point, y = sender_signal_intensity, group = celltype, color = celltype)) +                         # Line plot
  geom_smooth(se = FALSE, method = "loess") +  # Smoothed lines without standard error ribbon
  labs(title = "Sending Signal Intensity Across Time Points",
       x = "Time Point",
       y = "Sender Signal Intensity") +
  scale_color_manual(values = color.scheme) +  # Use the custom colors
  theme_minimal() +                       # Clean theme
  theme(
    legend.position = "right",            # Adjust legend position
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank()    # Remove minor grid lines
  )
ggplot(result.direct.total, aes(x = time_point, y = receiver_signal_intensity, group = celltype, color = celltype)) +                         # Line plot
  geom_smooth(se = FALSE, method = "loess") +  # Smoothed lines without standard error ribbon
  labs(title = "Receiving Signal Intensity Across Time Points",
       x = "Time Point",
       y = "Receiver Signal Intensity") +
  scale_color_manual(values = color.scheme) +  # Use the custom colors
  theme_minimal() +                       # Clean theme
  theme(
    legend.position = "right",            # Adjust legend position
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank()    # Remove minor grid lines
  )
dev.off()

# Ensure the celltype column is a factor with the desired order
result.secreted.total$celltype <- factor(result.secreted.total$celltype, levels = ordered_celltypes)
result.secreted.total$time_point <- factor(result.secreted.total$time_point, levels = c("p0", "p7", "p14", "p21"))
result.secreted.total$sender_signal_intensity <- as.numeric(result.secreted.total$sender_signal_intensity)
result.secreted.total$receiver_signal_intensity <- as.numeric(result.secreted.total$receiver_signal_intensity)
pdf("./COMMOT/plots/indirect_signal_intensity_over_time.pdf", onefile = T)
ggplot(result.secreted.total, aes(x = time_point, y = sender_signal_intensity, group = celltype, color = celltype)) +                         # Line plot
  geom_smooth(se = FALSE, method = "loess") +  # Smoothed lines without standard error ribbon
  labs(title = "Sending Signal Intensity Across Time Points",
       x = "Time Point",
       y = "Sender Signal Intensity") +
  scale_color_manual(values = color.scheme) +  # Use the custom colors
  theme_minimal() +                       # Clean theme
  theme(
    legend.position = "right",            # Adjust legend position
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank()    # Remove minor grid lines
  )
ggplot(result.secreted.total, aes(x = time_point, y = receiver_signal_intensity, group = celltype, color = celltype)) +                         # Line plot
  geom_smooth(se = FALSE, method = "loess") +  # Smoothed lines without standard error ribbon
  labs(title = "Receiving Signal Intensity Across Time Points",
       x = "Time Point",
       y = "Receiver Signal Intensity") +
  scale_color_manual(values = color.scheme) +  # Use the custom colors
  theme_minimal() +                       # Clean theme
  theme(
    legend.position = "right",            # Adjust legend position
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank()    # Remove minor grid lines
  )
dev.off()


###Figure 2F
library(tidyr)      # for pivot_wider
library(dplyr)      # for pipe %>%
library(pheatmap)
library(tibble)
commot_files <- list(
  p0 = "./COMMOT/results/commot_niches/average_niche_interaction_scores_p0_v2.csv",
  p7 = "./COMMOT/results/commot_niches/average_niche_interaction_scores_p7_v2.csv",
  p14 = "./COMMOT/results/commot_niches/average_niche_interaction_scores_p14_v2.csv",
  p21 = "./COMMOT/results/commot_niches/average_niche_interaction_scores_p21_v2.csv"
)

# Read and process data into a list
commot_niche_pathways <- lapply(commot_files, function(file) {
  df <- read.csv(file)
#  df %>%
#    mutate(Pathway = gsub("commot-user_database-", "", Pathway))
})

# Load CellChatDB.mouse data
data(CellChatDB.mouse)
LR_database_mouse <- CellChatDB.mouse$interaction

# Define categories
direct.cat <- c("ECM-Receptor", "Cell-Cell Contact")
indirect.cat <- c("Secreted Signaling")

# Process pathways
process_pathways <- function(df, categories) {
  df %>%
    filter(annotation %in% categories) %>%
    mutate(
      receptor.symbol = gsub(", ", "_", receptor.symbol),
      LRpair = paste(ligand.symbol, receptor.symbol, sep = "-")
    )
}

direct.pathway <- process_pathways(LR_database_mouse, direct.cat)
indirect.pathway <- process_pathways(LR_database_mouse, indirect.cat)

# Filter commot_niche_pathways for direct pathways
commot_niche_LR_direct <- lapply(commot_niche_pathways, function(df) {
  df %>%
    filter(ligand_receptor %in% direct.pathway$LRpair)
})
# Filter commot_niche_pathways for direct pathways
commot_niche_LR_indirect <- lapply(commot_niche_pathways, function(df) {
  df %>%
    filter(ligand_receptor %in% indirect.pathway$LRpair)
})
p0.niche<-c(1,3,4,5,6,7,8,10,9)
p7.niche<-c(0,1,3,4,5,6,7,8,10,2)
p14.niche<-c(0,1,3,4,5,6,7,8,10,11)
p21.niche<-c(0,1,3,4,5,6,7,8,10,11)
niche.time<-list("p0" = p0.niche,"p7" = p7.niche,"p14" = p14.niche,"p21" = p21.niche)

timepoint <- "p21"
df<-commot_niche_LR_indirect[[timepoint]][,c("ligand_receptor","niche","average_interaction_score_by_one")]
df<-df[df$niche %in% niche.time[[timepoint]],]
heatmap_data <- df %>%
  pivot_wider(names_from = niche, values_from = average_interaction_score_by_one, values_fill = list(average_interaction_score_by_one = 0)) %>%
  column_to_rownames("ligand_receptor")  # convert the LigandReceptor column to row names
heatmap_data[is.na(heatmap_data)] <- 0

# Convert the resulting data frame to a matrix
activity_matrix <- as.matrix(heatmap_data)

# Run PCA (center and scale data)
pca_result <- prcomp(t(activity_matrix), center = TRUE, scale. = TRUE)

# You can view a summary of the PCA results:
summary(pca_result)
pdf(paste0("../../../Manuscript/Figure 3/",timepoint,"indirect_niche_LR_PCA.pdf"), height = 5.5, width = 5)
# Basic scatter plot of PC1 vs PC2
plot(pca_result$x[, 1:2],
     xlab = "PC1",
     ylab = "PC2",
     main = "PCA Plot: Ligand-Receptor Activity by Niche",
     pch = 19, 
     col = xenium_12niche_color)
# Add labels for each dot. Here, we use the row names as the IDs.
text(pca_result$x[, 1], pca_result$x[, 2],
     labels = rownames(pca_result$x),  # or use another vector with IDs
     pos = 3,   # positions text above the point
     cex = 0.7) # adjust text size if needed
dev.off()

###Figure 2G
commot.p7.niche<-read.csv("./COMMOT/results/commot_niches/average_niche_interaction_scores_p7_v2.csv")

df<-dplyr::filter(commot.p7.niche,commot.p7.niche$ligand_receptor==unique(commot.p7.niche$ligand_receptor)[1])
df_transposed <- as.data.frame(t(df))
colnames(df_transposed)<-df_transposed["niche",]
df.final<-df_transposed["average_interaction_score_by_one",]
df.result<-df.final[FALSE,]
heatmap.list<-list()

pt <-"SEMA3"
df.result<-df.final[FALSE,]
df.pathway<-dplyr::filter(commot.p7.niche, commot.p7.niche$pathway %in% pt)
for (lr in unique(df.pathway$ligand_receptor)){
  df<-dplyr::filter(df.pathway,df.pathway$ligand_receptor==lr)
  df_transposed <- as.data.frame(t(df))
  colnames(df_transposed)<-df_transposed["niche",]
  df.final<-df_transposed["average_interaction_score_by_one",]
  rownames(df.final)<-lr
  df.result<-rbind(df.result,df.final)
}
df.result.sort<-df.result[,c(" 1", " 2"," 3", " 4", " 5"," 6", " 7"," 8", "10")]
df.result.sort[is.na(df.result.sort)] <- 0

df.result.num <- as.data.frame(lapply(df.result.sort, as.numeric))
rownames(df.result.num)<-rownames(df.result.sort)
colnames(df.result.num)<-colnames(df.result.sort)

# Scale data by row
scaled_df <- t(apply(as.matrix(df.result.num), 1, scale))
colnames(scaled_df)<-c("1","2","3", "4", "5","6", "7","8", "10")

heatmap.list[[pt]]<-Heatmap(scaled_df, 
                            name = "communication_Intensity",  # Label for the color scale
                            row_title = "Rows",  # Title for rows
                            column_title = paste0("communication_",pt),  # Title for columns
                            show_row_dend = TRUE,  # Show row dendrogram
                            show_column_dend = TRUE,  # Show column dendrogram
                            cluster_rows = TRUE,  # Cluster rows
                            cluster_columns = FALSE,  # Cluster columns
                            row_labels = rownames(scaled_df),
                            column_labels = colnames(scaled_df),
                            col = colorRampPalette(c("white", "#bdbdbd", "black"))(50))  # Color scale
heatmap.list[[pt]] = draw(heatmap.list[[pt]])
pdf("../../../Manuscript/Figure 3/p7_SEMA3_LR_communication_heatmap.pdf", onefile = T)

draw(heatmap.list[[pt]])

dev.off()

###Figure 2H
# Generate a color palette
color_palette <- viridis(10)

# Load Xenium datasets
timepoints <- c("p0", "p7", "p14", "p21")
xenium.obj <- lapply(timepoints, function(tp) readRDS(paste0("./robjs/ref.xenium.", tp, ".1011.rds")))
names(xenium.obj) <- timepoints
heatmap.list<-list()
timepoint<-"p7"

# Extract pathway names starting with "commot"
pathway <- grep("Sema3", colnames(xenium.obj[[timepoint]]@meta.data), value = TRUE)

# Generate heatmaps only if pathways are found
if (length(pathway) > 0) {
  heatmap.list <- lapply(pathway, function(pathwayname) {
    ImageFeaturePlot(xenium.obj[[timepoint]], features = pathwayname, max.cutoff = 1, size = 1, cols = color_palette, border.size = NA)
  })
  
  # # Save all plots to a PDF
  # pdf(file.path("./plots", paste0(timepoint, "_cell_indirect_communication_spatial_distribution.pdf")), onefile = TRUE)
  # invisible(lapply(heatmap.list, print))  # Use invisible to suppress unnecessary output
  # dev.off()
}
cropped.coords <- Crop(xenium.obj[[timepoint]][["fov.2"]], x = c(2750, 3250), y = c(2750, 3250), coords = "plot")
xenium.obj[[timepoint]][["zoom"]] <- cropped.coords
DefaultBoundary(xenium.obj[[timepoint]][["zoom"]]) <- "segmentation"
ImageFeaturePlot(xenium.obj[[timepoint]], fov= "zoom", features = "commot.user_database.Sema3a.Nrp1_Plxna4", max.cutoff = 1,size = 0.75, cols = color_palette,border.color = NA)
ImageFeaturePlot(xenium.obj[[timepoint]], fov= "zoom", features = "commot.user_database.Sema3c.Plxnd1", max.cutoff = 1,size = 0.75, cols = color_palette,border.color = NA)
ImageFeaturePlot(xenium.obj[[timepoint]], fov= "zoom", features = "commot.user_database.SEMA3", max.cutoff = 1,size = 0.75, cols = color_palette,border.color = NA)

###Figure 2J

library(viridis)
p0.xenium<-readRDS("./robjs/ref.xenium.p0.direct.1011.rds")
p7.xenium<-readRDS("./robjs/ref.xenium.p7.direct.1011.rds")
p14.xenium<-readRDS("./robjs/ref.xenium.p14.direct.1011.rds")
p21.xenium<-readRDS("./robjs/ref.xenium.p21.direct.1011.rds")

color_palette <- viridis(20)

pdf("../../../Manuscript/Figure 3/p7_cap_LR_visual.pdf")
ImageFeaturePlot(p7.xenium, fov= "fov.4", features = c("commot.user_database.Ptprm.Ptprm"), max.cutoff = 2,size = 0.75, cols = color_palette, border.size = NA,dark.background =F)
ImageFeaturePlot(p7.xenium, fov= "fov.4", features = c("commot.user_database.Mpzl1.Mpzl1"), max.cutoff = 2,size = 0.75, cols = color_palette, border.size = NA,dark.background =F)
ImageFeaturePlot(p7.xenium, fov= "fov.4", features = c("commot.user_database.Sema6a.Plxna4"), max.cutoff = 1,size = 0.75, cols = color_palette, border.size = NA,dark.background =F)
ImageFeaturePlot(p7.xenium, fov= "fov.4", features = c("commot.user_database.Dll4.Notch1"), max.cutoff = 1,size = 0.75, cols = color_palette, border.size = NA,dark.background =F)
ImageFeaturePlot(p7.xenium, fov= "fov.4", features = c("commot.user_database.Esam.Esam"), max.cutoff = 1,size = 0.75, cols = color_palette, border.size = NA,dark.background =F)
ImageFeaturePlot(p7.xenium, fov= "fov.4", features = c("commot.user_database.Efnb2.Ephb4"), max.cutoff = 2,size = 0.75, cols = color_palette, border.size = NA,dark.background =F)
dev.off()

###Figure 2L
p0.xenium<-readRDS("./robjs/ref.xenium.p0.1011.rds")
p7.xenium<-readRDS("./robjs/ref.xenium.p7.1011.rds")
p14.xenium<-readRDS("./robjs/ref.xenium.p14.1011.rds")
p21.xenium<-readRDS("./robjs/ref.xenium.p21.1011.rds")
pdf("../../../Manuscript/Figure 3/Kit_visual.pdf")
ImageFeaturePlot(p0.xenium, fov= "fov", features = c("commot.user_database.Kitl.Kit"), max.cutoff = 1,size = 1, cols = color_palette, border.size = NA,dark.background =F)
ImageFeaturePlot(p7.xenium, fov= "fov.2", features = c("commot.user_database.Kitl.Kit"), max.cutoff = 1,size = 1, cols = color_palette, border.size = NA,dark.background =F)
ImageFeaturePlot(p14.xenium, fov= "fov.3", features = c("commot.user_database.Kitl.Kit"), max.cutoff = 1,size = 1, cols = color_palette, border.size = NA,dark.background =F)
ImageFeaturePlot(p21.xenium, fov= "fov.4", features = c("commot.user_database.Kitl.Kit"), max.cutoff = 1,size = 0.75, cols = color_palette, border.size = NA,dark.background =F)
dev.off()


