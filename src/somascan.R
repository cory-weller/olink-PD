#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(ggplot2)
library(WGCNA)
library(foreach)
library(doMC)
registerDoMC(cores=4)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

DATADIR <- '/gpfs/gsfs9/users/CARD/projects/proteomicsXprogression_PD/data'
ANNO_FN <- 'PPMI_Project_151_pqtl_Analysis_Annotations_20210210.csv'
SOMASCAN_FN <- 'inputDF_rm_GBA.csv'

annotations <- fread(paste0(DATADIR,'/',ANNO_FN))
annotations[, c(3,8:61) := NULL]

ss_data <- fread(paste0(DATADIR,'/',SOMASCAN_FN))
ss_data[DUMMYREC=='HC_all', ARM := 'HC']
ss_data[DUMMYREC=='PD_all', ARM := 'PD']



annotations[, SOMA_SEQ_ID := gsub('-', '_', SOMA_SEQ_ID)]
annotations[, SOMA_SEQ_ID := paste0('anti_', SOMA_SEQ_ID)]
annotations[, SOMA_SEQ_ID := gsub('_3$', '', SOMA_SEQ_ID)]

ids <- intersect(annotations$SOMA_SEQ_ID, colnames(ss_data))
annotations <- annotations[SOMA_SEQ_ID %in% ids]
setkey(annotations, SOMA_SEQ_ID)

o <- foreach(gid=unique(annotations$SOMA_SEQ_ID), .combine='rbind') %do% {
    gene_symbols <- paste(annotations[gid]$TARGET_GENE_SYMBOL, collapse=';')
    data.table('SOMA_SEQ_ID'=gid, 'GENE_SYMBOLS'=gene_symbols)
}
setkey(o, SOMA_SEQ_ID)

dat <- ss_data[, .SD, .SDcols=c('PATNO','ARM',ids)]
setnames(dat, ids, )

patient_ids <- dat$PATNO
patient_arm <- dat$ARM

dat.t <- as.data.table(t(dat[, -c(1:2)]))
setnames(dat.t, as.character(dat$PATNO))

dat.t <- cbind(o, dat.t)

### good to go from here


# Assuming information columns come first, observation columns come after
# Define the # of information columns
n_info_cols <- 2
sample_id_col <- 1

info.dt <- dat.t[,.SD, .SDcols=(1:n_info_cols)]
col_names <- info.dt[[sample_id_col]]

obs.dt <- dat.t[,.SD, .SDcols=!(1:n_info_cols)]

row_names <- colnames(obs.dt)
datExpr0 <- t(as.data.frame(obs.dt))


# Organize data such that rows = samples; columns = genes
rownames(datExpr0) <- row_names
colnames(datExpr0) <- col_names


gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

if (!gsg$allOK) {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
        printFlush(paste("Removing genes:", paste(colnames(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
        printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}



sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


# Plot a line to show the cut
abline(h = 40, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 40, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================
# traits <- fread('reports/mds_scores.tsv')
# traits <- traits[participant_id %like% "^PP-"]
# traits[, ID := gsub('-', '', participant_id)]
# traits[,participant_id := NULL]
# traits <- traits[visit_month == 0]
# traits <- dcast(traits, visit_month + ID ~ variable, value.var='summary_score')
# traits <- traits[ID %in% rownames(datExpr)]
# datTraits <- as.data.frame(traits[, .SD, .SDcols=c('mds_i','mds_ii','mds_iii')])
# rownames(datTraits) <- traits[['ID']]


# collectGarbage();


# #=====================================================================================
# #
# #  Code chunk 8
# #
# #=====================================================================================


# # Re-cluster samples
# sampleTree2 = hclust(dist(datExpr), method = "average")
# # Convert traits to a color representation: white means low, red means high, grey means missing entry
# traitColors = numbers2colors(datTraits, signed = FALSE);
# # Plot the sample dendrogram and the colors underneath.
# plotDendroAndColors(sampleTree2, traitColors,
#                     groupLabels = names(datTraits), 
#                     main = "Sample dendrogram and trait heatmap")


### CLUSTERING


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
# Choose lowest power that hits 90% R-squared
chosen_power <- sft$fitIndices[which(sft$fitIndices$SFT.R.sq > 0.9),][1,1]
# got NA so manually choose 10, best option


if(is.na(chosen_power)) {
    cat('No power thresholds exceed R-squared of 90%!\n')
}

# chosen_power <- 3

# Read the docs for this function and optimize!
net = blockwiseModules(datExpr, power = chosen_power,
                       TOMType = "unsigned", minModuleSize = 15,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM", 
                       verbose = 3)


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)

# module Eigengenes
MEs = net$MEs;

modules <- data.frame('group'=names(net$colors), 'module'=as.vector(net$colors))
modules <- modules[order(modules$module),]

cardio_genes <- unique(dat[olink_panel=='cardiometabolic'][['UniProt']])
neuro_genes <- unique(dat[olink_panel=='neurology'][['UniProt']])
onco_genes <- unique(dat[olink_panel=='oncology'][['UniProt']])
infla_genes <- unique(dat[olink_panel=='inflammation'][['UniProt']])

setDT(modules)
sum(modules[module==0][['group']] %in% )  
sum(modules[module==1][['group']] %in% neuro_genes)
sum(modules[module==2][['group']] %in% cardio_genes)
sum(modules[module==3][['group']] %in% )

geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "olink-PPMI.RData")