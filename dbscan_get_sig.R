# ---------------------------------------------------------------------------------------------
# Function for generating gene signatures based on the t.test output of "dbscan_programs"
# --------------------------------------------------------------------------------------------- 

# - dbscan_programs_output = outuput from "dbscan_programs"
# - cell_line = cell line name
# - max_genes = maximum number of genes in the signature
# - p_val = t.test p val cutoff
# - log2_fc = differential expression fold change cutoff
# - max_size = maximum percentage of cells in each cluster


# Returns a list indicating the cells in each clusters and a dataframe with t.test results comparing gene expression of cell inside vs. outside the respective cluster.
# If no clusters were identified, returns NA

dbscan_get_sig <- function(dbscan_programs_output, cell_line, max_genes = 50, p_val = 0.001, log2_fc = 1, max_size=0.9) {
    a <- dbscan_programs_output[["clusters_sig"]]
    a <- lapply(a, function(x) x[x[,"log2.FC."] >= log2_fc & x[,"ttest_p"] < p_val,])
    a <- lapply(a, function(x) rownames(x)[order(x[,"log2.FC."], decreasing = T)][1:max_genes])
    a <- lapply(a, function(x) x[!is.na(x)])
    
    meta <-  readRDS("CCLE_heterogeneity_Rfiles/CCLE_metadata.RDS")           
    b <-  dbscan_programs_output[["clusters_cells"]]
    b <-  sapply(b, function(x) (length(x)/meta[cell_line,"n_cells"]) < max_size)
     
    return(a[which(b)])             
}
