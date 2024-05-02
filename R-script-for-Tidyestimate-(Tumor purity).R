##--------------- ESTIMATE tumor purity assessment (via tidyestimate)---------------##
## install tidyestimate
devtools::install_github("KaiAragaki/tidy_estimate")
library(tidyestimate)

## Retrieve input for tidyestimate (It can be data.frame, tibble, or matrix) - gene expression table
head(geneCount_uq)

## Run tidyestimate
Tumorscores <- geneCount_uq |>
  filter_common_genes(id = "hgnc_symbol", tell_missing = FALSE, find_alias = TRUE) |> 
  estimate_score(is_affymetrix = TRUE)
Tumorscores

## Save data in .txt file
write.table(Tumorscores, file = "location/to/save/TumorPure_ESTIMATE.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)

## Cellular component by xCell can perform online at http://xCell.ucsf.edu/
## Create input for xCell by following command
write.table(geneCount_uq, file = "location/to/save/geneCount_uq.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)

#------------------------------ End ------------------------------#