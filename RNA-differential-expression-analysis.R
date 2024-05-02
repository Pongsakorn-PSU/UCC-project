#----------R script for differential expression analysis----------#
# install the package
install.packages("matrixStats")
BiocManager::install("IsoformSwitchAnalyzeR")
BiocManager::install("DEXSeq")
BiocManager::install("SummarizedExperiment")
BiocManager::install("MatrixGenerics",force = TRUE)

library(matrixStats)
library(DEXSeq)
library(SummarizedExperiment)
library(MatrixGenerics)
library(IsoformSwitchAnalyzeR)

# Import gene expression files
stringTieQuant <- importIsoformExpression(
  parentDir = "/path/to/gene/expr/directory", 
  addIsofomIdAsColumn = FALSE,
  readLength=150 
)

# Create design matrix
myDesign <- data.frame(
  sampleID = colnames(stringTieQuant$abundance),
  condition = rep(c('UCC','Control','Tumor'), times=c(9,7,15))
)

# Create switchAnalyzeRlist
switchAnalyzeRlist <- importRdata(
  isoformCountMatrix   = stringTieQuant$counts,
  isoformRepExpression = stringTieQuant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "path/to/file/merge.gtf",
)

## Extract gene count matrix
geneCount <- extractGeneExpression(
  switchAnalyzeRlist,
  extractCounts = TRUE # set to FALSE for abundances
)
head(geneCount)

## Filter out lncRNA
V38.gtf <- rtracklayer::import("path/to/ref/gencode.v38.annotation.gtf")
head(V38.gtf)
table(V38.gtf$gene_type)
lncrna <- unique(V38.gtf[V38.gtf$gene_type=="lncRNA"]$gene_name)
lncrna
geneCount <- geneCount[!(geneCount$gene_name %in% lncrna), ]
dim(geneCount)
head(geneCount)

# Combine duplicate row
head(geneCount)
length(geneCount$gene_name)
length(unique(geneCount$gene_name))
geneCount[1] <- list(NULL)

geneCount_uq <- geneCount %>%
  group_by(gene_name) %>%
  summarise(across(colnames(geneCount)[-(1)], sum)) %>%
  as.data.frame()

gene_name_uq <- geneCount_uq[,1]
head(gene_name_uq)
geneCount_uq[1] <- list(NULL)

## Delete gene id and gene name colume on counts
counts_mat <- apply(as.matrix.noquote(geneCount_uq),
                    2,
                    as.numeric)
class(counts_mat)
row.names(counts_mat) <- gene_name_uq #assign gene name as row name
head(counts_mat)

#----------Differential expression analysis with edgeR----------#

##Install and load edgeR software and rtracklayer
BiocManager::install('edgeR')
BiocManager::install("rtracklayer")
library(edgeR)
library(rtracklayer)
library(dplyr)

## create design matrix
group <- myDesign$condition
group
design <- model.matrix(~0+group)
design

#Single factor analysis with edgeR
y <- DGEList(counts=counts_mat, group=group)
head(y)

# log transform and normalization
logcpm <- cpm(y, prior.count = 2, log=TRUE)
class(logcpm) #Should be numeric matrix
head(logcpm)
colnames(logcpm)

## Export count matrix as .csv file
write.csv(logcpm, file='location/to/save/gcount.csv', row.names = TRUE )

# Calculate normalization factor
genexp <- calcNormFactors(y)
head(genexp)

#GLM Common dispersion
genexp <- estimateGLMCommonDisp(genexp, design)
head(genexp)

#Estimate GLM trended dispersions
genexp <- estimateGLMTrendedDisp(genexp, design)
head(genexp)

#Tagwise dispersion of each gene
genexp <- estimateGLMTagwiseDisp(genexp, design)
head(genexp)

#Testing for DE genes (BBz, CD28, CD28_40, CD79A_40)
fit <- glmQLFit(genexp, design)
fit
genediff <- glmQLFTest(fit, contrast=c(-1,0,1))
topTags(genediff)
head(genediff)
class(genediff)

#Constructing a table showing differential expression with edgeR, trimming out those p-value >= 0.03
genediff$table <- cbind(genediff$table, FDR=p.adjust(genediff$table$PValue, method = 'fdr'))
head(genediff$table)

#Select gene with p-value >=0.03 and sort table by p-value
gsign <- genediff$table[genediff$table$PValue<0.1,]
gsign <- gsign[order(gsign$PValue),]

dim(gsign)
head(gsign, 10)

# Extract Up and Down regulated gene
GUp <- gsign %>% dplyr::filter(logFC>0)
GDown <- gsign %>% dplyr::filter(logFC<(0))
head(GUp)
dim(GUp)
dim(GDown)

write.table(rownames(GUp), file = "location/to/save/RNAsign_up.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE)
write.table(rownames(GDown), file = "location/to/save/RNAsign_down.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE)

##Write table .txt
write.table(gsign, file = "location/to/save/gsign_urine_vs_control.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

#------------------------------ End ------------------------------#
