#Set working directory to folder containing counts and design data
setwd("C:/Users/desqz2/Desktop/ARPM")
#Load required packages
library(DESeq2)
library(edgeR)
library(WGCNA)
library(dendextend)
library(heatmap3)
library(limma)
library(dplyr)
library(tidyverse)
library(scales)
library(sva)
library(ggrepel)
library(patchwork)
library(BIGverse)

### Load raw OPAR count matrix
raw_counts <- read.csv("combined_count_matrix_1_modified.csv", header=T, row.names=1)
head(raw_counts)
design<-read.csv("DesignMatrix.csv", row.names = 1)
#Normalize Counts with DESeq using variance stabilizingg transfomation
dds <- DESeqDataSetFromMatrix(dat ,dt , design = ~Type)
dim(dds) 
dds <- DESeq(dds)
a <- varianceStabilizingTransformation(dds)
b <- getVarianceStabilizedData(dds)
b1 <- rowVars(b)
summary(b1)
q00_wpn <- quantile( rowVars(b)) 
expr_normalized_Data <- b[ b1 > q00_wpn, ]
xx<-expr_normalized_Data
dim(xx)
write.csv(xx, "Normalized_data_matrix.csv")
dt2$CV<-apply(dt2,1, function(x) sd(x) / mean(x) * 100) # CV of each gene

# Select hypervariable genes based CV greater than 4 % 
dat3<-dat2[dat2$CV > 4,]

#WGCNA analysis 
dat<-read.csv("dat3.csv", row.names = 1) # construct WGCNA based Hypervariable genes
dim(dat3)
input_mat = t(dat3)
# Choose a set of soft threshold parameters
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(input_mat, powerVector = powers, verbose = 5) 
# Scale-free topology fit index as a function of the soft-thresholding power
#pdf(file = "2-n-sft.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


picked_power = 5
temp_cor <- cor       
cor <- WGCNA::cor         
netwk <- blockwiseModules(input_mat,
                          
                          power = picked_power,                
                          networkType = "signed",
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          minModuleSize = 50,
                          maxBlockSize = 5000,
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          numericLabels = T,
                          verbose = 3)




module_eigengenes <- netwk$MEs
netwk$MEs
# Print out a preview
head(module_eigengenes)

df<-module_eigengenes
write.csv(df,"module_eigengeneMerge.csv")

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

netwk$colors[netwk$blockGenes[[1]]]
table(netwk$colors)
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]
module_df
write.csv(module_df, "modGSE107361.csv")





# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

module_order


MEs0
write.csv(MEs0, "MEs0_GSEMerge.csv")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

MEList = moduleEigengenes(input_mat, colors = mergedColors)
MEs = MEList$eigengenes
plotEigengeneNetworks(MEs0, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))

#**************************************************
# Define numbers of genes and samples
nGenes = ncol(input_mat);
nSamples = nrow(input_mat);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(input_mat, mergedColors)$eigengenes
MEs = orderMEs(MEs0)
dim(MEs)
MEs$MEgreen

# Read clincial traits of asthma and control subjects
bac_traits = read.csv("GSE130588_clinicalData2.csv", row.names = 1)
bac_traits[,-1] # for CSV file
bac_traits
rownames(bac_traits)
rownames(MEs)
dd<-rownames(MEs)
write.csv(dd, "GSE121212_MEs.csv")

# sample names should be consistent in eigen genes and traits !!!!
bac_traits = bac_traits[match(rownames(MEs), rownames(bac_traits)), ]

table(rownames(MEs) == rownames(bac_traits))
# Calculate pearson correlation coefficients between module eigen-genes and traits
moduleTraitCor = cor(MEs, bac_traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
moduleTraitCor 
moduleTraitPvalue

xx2<-moduleTraitPvalue
moduleTraitCor
xx3<-moduleTraitCor

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 1), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(bac_traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))




# Gene Significance and Module Membership

# Define variable Lesional containing from input-dat
Lesional = as.data.frame(bac_traits$Lesional);
names(Lesional) = "Lesional"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
modNames
Lesional
geneModuleMembership = as.data.frame(cor(input_mat, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(input_mat, Lesional, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Lesional), sep="");
names(GSPvalue) = paste("p.GS.", names(Lesional), sep="");

MEs$MEblue
module = "blue"
column = match(module, modNames);
column

moduleGenes = mergedColors==module;
moduleGenes

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Asthma Status",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 0.9, cex.lab = 1, cex.axis = 1, col = module) 



# For count data, we used DESeq2 package for DEGs analysis
# DEGs analysis
dat<-read.csv("Count_GEOData_data_matrix.csv", row.names = 1)
dt<-read.csv("designmatrix.csv", row.names = 1)
#Create a Datasheet from count matrix and labels
dds <- DESeqDataSetFromMatrix(dat ,dt , design = ~Type)
dim(dds)
dds$Type <- relevel(dds$Type, ref = "C") # C->set reference control subjects
dds <- DESeq(dds)
res <- results(dds)
head(res)
class(res)
write.csv(res, "DEGs_GSE157194__all.csv")

# Limma package was used for DEG analysis for normalized GEO dataset 

dat<-read.csv("Normalized_data_matrix.csv", row.names = 1)

design.mat<-read.csv("design matrix.csv")
design.mat<-design.mat[,-1]
design.mat
dim(design.mat)
contrast.mat<-matrix(c(1,-1), ncol = 1)
dimnames(contrast.mat)<-list(c('A', 'C'), "Diff")
sample<-factor(rep(c("A", "C"), c(21,410)))
design.mat<-model.matrix(~0+sample)
colnames(design.mat)<-levels(sample)
design.mat
dim(design.mat)
contrast.mat<-makeContrasts(diff = A - C, levels = design.mat)
contrast.mat
fit<-lmFit(M, design.mat)
fit2<-contrasts.fit(fit, contrast.mat)
fit3<-eBayes(fit2)
deg<-topTable(fit3) # top ranked DEGs
deg1 <- topTable(fit3, n=Inf, coef=1,adjust.method="BH" )
DEG1 <- as.data.frame(deg1)
View(DEG1)
write.csv(DEG1, "DEGs_based on_limma.csv")

# batch effect adjustment based on Surrogate Variable Analysis (SVA)
datx<-read.csv("data with batch_effect.csv", row.names = 1)
b<-read.csv("Design_matrix.csv")
b$batch<-as.factor(b$batch)
dim(b)
dim(datx)
batch = b$batch
batch
# parametric adjustment
combat_gdata = ComBat(dat= datx, batch= batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
combat_gdata
write.csv(combat_gdata, "Batch_adjusted_by_sva.csv")
# non-parametric adjustment, mean-only version
combat_edata2 = ComBat(dat=datx, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TR)
combat_edata2 
write.csv(combat_gdata2, "Batch_adjusted_by_sva.csv")
















