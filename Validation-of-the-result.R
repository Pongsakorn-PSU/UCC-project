###------------------- Regression model prediction -------------------###

# Retrieve expression data of this study
exprs(eset_norm)
head(exprs(eset_norm))

# Create expression matrix
x <- t(exprs(eset_norm))
x <- data.frame(cbind(x,condition = pData(eset_norm)$condition))
head(x)

# Assign condition (UCC = 1, Control = 0)
for (i in 1:length(x$condition)) {
  if(x$condition[i] == "UCC"){
    x$condition[i] <- 1
  } else if (x$condition[i] == "Control") {
    x$condition[i] <- 0
  }
}

# Change all column to numeric
for (i in 1:ncol(x)) {
  if(class(x[,i]) == "character"){
    x[,i] <- as.numeric(x[,i])
  }
}
x <- as.data.frame(x)

# Retrieve validation dataset from public repository
Expr_validate <- read.table('path/to/project/PXD010260/Expression_table.txt',
                            header=TRUE, sep="\t",na.strings = "NULL")
head(Expr_validate)
rownames(Expr_validate) <- Expr_validate$Accession
Expr_validate <- Expr_validate[,-c(1:2)]
head(Expr_validate)

Expr_validate <- log2(Expr_validate+1) #log2 transform
Expr_validate <- data.frame(t(Expr_validate)) # transpose

for (i in 1:ncol(Expr_validate)) {
  if(class(Expr_validate[,i]) == "character"){
    Expr_validate[,i] <- as.numeric(Expr_validate[,i])
  }
}

Expr_validate <- scale(Expr_validate)

# Retrieve data of Co-sign gene
head(CoSign_final)

# Create model performance for each specific protein
logit_fit <- glm(condition ~ Q05823, family = "binomial", # Q05823 is the uniprot ID of interested protein, which is changeable
                 data = x)
summary(logit_fit)

logit_performance <- x |> select(condition, Q05823) |> 
  mutate(predicted_group = predict(logit_fit, newdata = x,
                                   type = "response")) |> 
  mutate(binomial_response = 
           (predicted_group^condition)*(1-predicted_group)^(1-condition)) |> 
  mutate(saturated_response = 1) |> 
  mutate(null_response = 0.5)

logit_likelihood <- logit_performance |> 
  summarize(across(contains("response"), ~sum(log(.x))))
logit_likelihood

# Prediction
new_group <- data.frame(Q05823 = Expr_validate[,"Q05823"])
new_group
new_group <- new_group |> 
  mutate(predicted_response = predict(logit_fit, 
                                      newdata = new_group, type = "response")) %>% 
  mutate(predicted_response2 = ifelse(predicted_response <= 0.5, 0, 1)) %>% 
  mutate(real_response = c(rep(1, 8), rep(0,8)))
table(new_group$predicted_response2, new_group$real_response)
new_group

# ROC analysis
library("pROC")
roc <- plot.roc(new_group$real_response,new_group$predicted_response)
roc
plot(roc)
