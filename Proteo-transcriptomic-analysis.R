##---------------Integrated proteo-transcriptomic analysis---------------##
# 1) Matching Protein and mRNA annotation
# Retrieve gene expression data
head(logcpm)
dim(logcpm)
# Retrieve protein expression data
head(ProtCounts_mat)
ncol(ProtCounts_mat)
# Select only data from shared sample in both data sets
RiP <- colnames(logcpm[,colnames(logcpm) %in% colnames(ProtCounts_mat)])
RiP
PiR <- colnames(ProtCounts_mat[,colnames(ProtCounts_mat) %in% colnames(logcpm)])
PiR
sample <- intersect(RiP,PiR) #find sample present in RNA and Protein data
sample
# Subset gene and prot count table which presented in shared samples
PC <- ProtCounts_mat[,sample]
GC <- logmRNA[,sample]
# Retrieve protein name and change to gene name
write.table(Prot_name, file = "location/to/save/prot_name.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE)
# Go search ID from Uniprot Database (https://www.uniprot.org/id-mapping)
# Save mapping gene name in prot_to_genename.tsv
# Retrieve mapping gene name
MappedG <- read.table(file='location/to/prot_to_genename.tsv')
head(MappedG)
MappedG <- MappedG[-1,]
row.names(MappedG) <- MappedG[,1]
colnames(MappedG) <- c("Prot_name","Gene_name")
head(MappedG)

# 2) Convert protein name to gene name
Name <- data.frame(Prot_name = rownames(PC),
                   Gene_name = 1:nrow(PC))

for (i in 1:nrow(Name)) {
  Name$Gene_name[i] <- MappedG$Gene_name[which(MappedG$Prot_name == Name$Prot_name[i])]
}
rownames(PC) <- Name$Gene_name

# 3) Urine-tissue mRNA correlation
# Retrieve gene expression data of urine and cancer tissue
dim(logcpm)
mRNAcount <- logcpm[(rownames(logcpm) %in% mrna), ]
mRNAcount <- mRNAcount %>% as.data.frame()

# Combine duplicate gene by merging expression level of their transcript
mRNAcount <- mRNAcount %>%
  group_by(gene_name) %>%
  summarise(across(colnames(geneCount)[-(1)], sum)) %>%
  as.data.frame()

mrna_uq <- mRNAcount[,1]
mrna_uq
row.names(mRNAcount) <- mrna_uq
mRNAcount[1] <- list(NULL)
head(mRNAcount)

# Log transformation and as.matrix numeric
logmRNA <- log2(mRNAcount+1)
head(logmRNA)
colnames(logmRNA)
logmRNA <- apply(as.matrix.noquote(logmRNA),
                 2,
                 as.numeric)
row.names(logmRNA) <- mrna_uq

# Select gene expression of urinary cell and tissue 
logmRNA_u <- subset(logmRNA, select = c("sample ID (column name) of urinary sample"))
logmRNA_t <- subset(logmRNA, select = c("sample ID (column name) of tissue sample")
logmRNA_c <- subset(logmRNA, select = c("sample ID (column name) of control sample")

# subtract with median of normal
logmRNA_u <- logmRNA_u - apply(logmRNA_c, 1, median)
logmRNA_t <- logmRNA_t - apply(logmRNA_c, 1, median)

# Perform Spearman correlation
res <- c()
x <- c()
pval <- c()

for (i in 1:nrow(logmRNA_t)) {
  res <- c(res, cor(logmRNA_u[i,], logmRNA_t[i,], method = "spearman"))
  x <- cor.test(logmRNA_u[i,], logmRNA_t[i,], method = "spearman")
  pval <- c(pval, x$p.value)
}

Cor_coef <- cbind(gene_name = rownames(logmRNA_t), Spearman_coef = res, pval)
Cor_coef <- data.frame(Cor_coef)
Coefsign_ut <- Cor_coef %>% dplyr::filter(pval<0.1)
dim(Coefsign_ut)
head(Coefsign_ut)

# Save data
write.table(Cor_coef, file = "F:/Pongsakorn/BCa/Urine_tissue_RNA_corr.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE)

# 4) Identify gene with dysregulation in both RNA and protein level
head(genediff$table) #Gene expression profile
head(toptable) #Protein expression profile
toptable <- cbind(Gene_name = 1:nrow(toptable),toptable)

for (i in 1:nrow(toptable)) {
  toptable$Gene_name[i] <- MappedG$Gene_name[which(MappedG$Prot_name == toptable$Prot_name[i])]
}
# Create CoSigntable
CoSigntable <- data.frame(Prot_name = toptable$Prot_name,
                          Gene_name = toptable$Gene_name,
                          logFC_prot = c(1:nrow(res)),
                          Pval_prot = c(1:nrow(res)),
                          logFC_rna = c(1:nrow(res)),
                          Pval_rna = c(1:nrow(res)))
head(CoSigntable)
# Fill the data
for (i in 1:nrow(CoSigntable)) {
  CoSigntable$Pval_prot[i] <- toptable$P.Value[which(toptable$Prot_name == CoSigntable$Prot_name[i])]
  CoSigntable$logFC_prot[i] <- toptable$logFC[which(toptable$Prot_name == CoSigntable$Prot_name[i])]
  CoSigntable$Pval_rna[i] <- genediff$table$PValue[which(rownames(genediff$table) == CoSigntable$Gene_name[i])]
  CoSigntable$logFC_rna[i] <- genediff$table$logFC[which(rownames(genediff$table) == CoSigntable$Gene_name[i])]
}
# Identify gene with significant expression in both RNA and protein level (p<0.1)
CoSign_final <- CoSigntable %>% dplyr::filter(Pval_prot>1) %>% dplyr::filter(Pval_rna>1)
CoSign_final
