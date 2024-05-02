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

# 3) mRNA-Protein correlation
# Subset GC by shared gene name in both data set
final_p <- PC[(rownames(PC) %in% intersect(rownames(PC),rownames(GC))),]
final_g <- GC[(rownames(GC) %in% intersect(rownames(PC),rownames(GC))),]
nrow(final_g)
nrow(final_p)
# Log2 transform of both data set (if retrieve from logmRNA, no need to transform final_g)
final_g <- apply(as.matrix.noquote(final_g), 2, as.numeric)
final_p <- log2(final_p+1)
final_g <- log2(final_g+1)
class(final_g)
class(final_p)
# Rearrange row name in the same order
head(row.names(final_g))
head(row.names(final_p))
head(final_p)
final_g <- final_g[rownames(final_p),]
head(final_g)
# subtract with median of normal
final_g <- final_g[,1:6] - apply(final_g[,7:11], 1, median) #sample number 7-11 is normal sample
final_p <- final_p[,1:6] - apply(final_p[,7:11], 1, median)
# Spearman correlation
res <- c()
x <- c()
pval <- c()

for (i in 1:nrow(final_g)) {
  res <- c(res, cor(final_g[i,], final_p[i,], method = "spearman"))
  x <- cor.test(final_g[i,], final_p[i,], method = "spearman")
  pval <- c(pval, x$p.value)
}

Cor_coef <- cbind(gene_name = rownames(secret_g), Spearman_coef = res, pval)
Cor_coef <- data.frame(Cor_coef)

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
