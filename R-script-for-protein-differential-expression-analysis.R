#------------------- Protein Expression data analysis -------------------#

# 1) Install required package
install.packages("BiocManager")
BiocManager::install("Biobase")
BiocManager::install("limma")

library(Biobase)
library(dplyr)

# 2) import expression data
# Protein expression table was obtained from Protein Pilot v.5.0.2.0 software
ProtCounts <- read.table('path/to/protein/expression/table/file.txt',
                         header=TRUE, sep="\t",na.strings = "NULL")
# Modifying row name and column name of the table
ProtCounts <- ProtCounts[-1,]
head(ProtCounts)
colnames(ProtCounts)
row.names(ProtCounts)
rname <- ProtCounts[,1]
head(rname)
Uniprot_ID <- sapply(strsplit(rname, split='|', fixed=TRUE), 
                     function(x) (x[2]))
Uniprot_ID
Prot_name <- sapply(strsplit(rname, split='|', fixed=TRUE), 
                    function(x) (x[3]))
Prot_name <- sapply(strsplit(Prot_name, split='_', fixed=TRUE), 
                    function(x) (x[1]))
Prot_name

# Assign row name
ProtCounts[,1] <- NULL
ProtCounts_mat <- apply(as.matrix.noquote(ProtCounts), 2,
                        as.numeric)
rownames(ProtCounts_mat) <- Prot_name
head(ProtCounts_mat) # Check the row name
class(ProtCounts_mat) # check the class of data

# Create phenotype data
condition <- gsub('BC.*','UCC', colnames(ProtCounts_mat)) # BC = bladder cancer
condition <- gsub('NC.*','Control', condition) # NC = Normal control
condition
table(condition)
length(condition)

Phenotype <- data.frame(sample_ID = colnames(ProtCounts_mat),
                        condition = condition,
                        stage = rep("early", each=length(condition)),
                        recurrent = rep("no_recur", each=length(condition)),
                        type = rep("MIBC", each=length(condition)),
                        UD = rep("UCC", each=length(condition)),
                        row.names = colnames(ProtCounts_mat))
# Assign sample group
latestage <- c("BC29","BC31","BC45","BC46","BC49","BC51","BC53","BC61")
recur <- c("BC13","BC24","BC42","BC46","BC48","BC50","BC53")
NMIBC <- c("BC13","BC19","BC23","BC30","BC31","BC34","BC39","BC40","BC42","BC43","BC44","BC45","BC46","BC48","BC49","BC50","BC51","BC52","BC53","BC55","BC56","BC57","BC60")
NC <- grep("NC", colnames(ProtCounts_mat), value=TRUE)

# Fill the phenotype data of patients
for (i in 1:nrow(Phenotype)) {
  if(Phenotype$sample_ID[i] %in% latestage){
    Phenotype$stage[i] <- c("late")
  } else if (Phenotype$condition[i] == "Control") {
    Phenotype$stage[i] <- c("Control")
  }
  if(Phenotype$sample_ID[i] %in% recur){
    Phenotype$recurrent[i] <- c("recur")
  }
  if(Phenotype$sample_ID[i] %in% NMIBC){
    Phenotype$type[i] <- c("NMIBC")
  }
  if(Phenotype$sample_ID[i] %in% NC){
    Phenotype$UD[i] <- sampleID$V3[which(Phenotype$sample_ID[i] == sampleID$V2)]
  }
}
Phenotype

# Create phenotype data of NMIBC group
Phenotype_NMIBC <- Phenotype %>% dplyr::filter(type == "NMIBC")
Phenotype_NMIBC

# Creat feature table
library(stringr)
library(dplyr)
library(tidyr)
feature <- data.frame(Prot_name = Prot_name,
                      Uniprot_ID = Uniprot_ID)
row.names(feature) <- Prot_name
head(feature)
dim(feature)

# 3) Create ExpressionSet object
library(limma)
eset <- ExpressionSet(assayData = ProtCounts_mat, # x is expression table
                      phenoData = AnnotatedDataFrame(Phenotype), # p is phenotype table
                      featureData = AnnotatedDataFrame(feature)) # f is feature

eset_NMIBC <- ExpressionSet(assayData = ProtCounts_mat[,Phenotype_NMIBC$sample_ID], # x is expression table
                      phenoData = AnnotatedDataFrame(Phenotype_NMIBC), # p is phenotype table
                      featureData = AnnotatedDataFrame(feature))
# View the number of feature (rows) and sample (columns)
dim(eset)
pData(eset)

# 4) Create design matrix
condition
pData(eset)

design <- model.matrix(~0+condition, data = pData(eset))
colnames(design) <- c("Control","UCC")
design_NMIBC <- model.matrix(~0+recurrent, data = pData(eset_NMIBC))
colnames(design) <- c("no_recur","recur")

# 5) Expression data normalization
# Create new ExpressionSet to store normalized data
eset_norm <- eset
eset_norm_NMIBC <- eset_NMIBC

# Log2 transform normalization
exprs(eset_norm) <- log2(exprs(eset_norm)+1)
exprs(eset_norm)
exprs(eset_norm_NMIBC) <- log2(exprs(eset_norm_NMIBC)+1)
exprs(eset_norm_NMIBC)

# 6) DEP with limma
library(limma)

# Fit the model
fit <- lmFit(eset_norm, design)
fit_NMIBC <- lmFit(eset_norm_NMIBC, design)
fit$design
fit_NMIBC$design

cont.matrix <- makeContrasts(UCCvsControl=UCC-Control, levels = design)
cont.matrix_NMIBC <- makeContrasts(Recurvsno_recur=recur-no_recur, levels = design)
cont.matrix
cont.matrix_NMIBC

fit2 <- contrasts.fit(fit, cont.matrix)
fit2_NMIBC <- contrasts.fit(fit_NMIBC, cont.matrix_NMIBC)

# Calculate the t-statistics
fit2 <- eBayes(fit2)
fit2_NMIBC <- eBayes(fit2_NMIBC)

toptable <- topTable(fit2, adjust="BH", number="inf", sort.by = "P")
toptable_NMIBC <- topTable(fit2_NMIBC, adjust="BH", number="inf", sort.by = "P")

head(toptable)
head(toptable_NMIBC)

# Filter differently expressed protein (p.value < 0.1)
SigBCa_vs_NORM <- toptable %>% dplyr::filter(P.Value<0.1)
SigBCavsNORM
Recur_vs_No_recur <- toptable_NMIBC %>% dplyr::filter(P.Value<0.1)
Recur_vs_No_recur

# Select upregulated and downregulated protein
Up <- SigBCavsNORM %>% dplyr::filter(logFC>1)
Down <- SigBCavsNORM %>% dplyr::filter(logFC<(-1))
dim(Up)
dim(Down)

#------------------- End -------------------#