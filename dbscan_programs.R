# ---------------------------------------------------------------------------------------------
# Function for getting heterogeneity programs using tSNE and density-based clustering (DBSCAN)
# --------------------------------------------------------------------------------------------- 

# - tpm = CPM expression matrix (rows = genes, columns = cells)
# - initial_dims = the number of dimensions that should be retained in the initial PCA step of the tSNE algorithm
# - perplexity = perplexity parameter to be used by the tSNE algorithm (should not be bigger than 3 * perplexity < ncol(cpm) - 1)
# - max_iter = number of iterations to be used by the tSNE algorithm
# - eps = size of the epsilon neighborhood to be used by the DBSCAN algorithm
# - minpts = number of minimum points in the eps region (for core points) to be used by the DBSCAN algorithm
# - min_size = minimum percentage of cells in each cluster
# - max_size = maximum percentage of cells in each cluster


# Returns a list indicating the cells in each clusters and the expression program of each cluster
# If no clusters were identified, returns NA

library(Rtsne)
library(dbscan)

dbscan_programs <- function(cpm, initial_dims=50, perplexity = 35, max_iter=7000, eps = 1.8, minpts=10, min_size= -Inf, max_size = Inf, seed=1) {
  
  # Processing cpm data
  CP100K_log <- log2((cpm/10) + 1)
  CP100K_log <- CP100K_log[apply(CP100K_log, 1, function(x) length(which(x > 3.5)) > ncol(CP100K_log)*0.02),]
  CP100K_log <- CP100K_log - rowMeans(CP100K_log)
    
  if(3*perplexity >= ncol(CP100K_log) - 1) {
    perplexity <- (ncol(CP100K_log) - 2)/3
    warning(paste("perplexity is too large - using ", (ncol(CP100K_log) - 2)/3, sep = ""))
  } 
  
  # Performing tSNE                                       
  set.seed(seed = seed)
  tsne <-  Rtsne(t(CP100K_log), dims = 2, pca = T, max_iter = max_iter, perplexity = perplexity, initial_dims = initial_dims)$Y
  rm(.Random.seed, envir=.GlobalEnv)
  
  # Perfoming density-based clustering 
  ds <-  dbscan(tsne, eps = eps, minPts = minpts)$cluster
  
  # Filtering clusters based on size 
  clusters <- split(colnames(CP100K_log), as.character(ds))
  clusters <- clusters[names(clusters)!="0"]
  
  if(length(clusters) == 0) return(NA) 
  
  clusters_final <- clusters[sapply(clusters, function(x) length(x) > ncol(CP100K_log)*min_size & length(x) < ncol(CP100K_log)*max_size )] 
  
  if(length(clusters_final) == 1) return(NA) 
  
  # Getting differentially expressed genes in each cluster 
  clusters_sig <- list()
  
  for(i in names(clusters_final)) {
    a <- clusters_final[[i]] # gets cells in the cluster
    if(length(a) > (ncol(CP100K_log)-2)) {
      clusters_sig[[i]] <- NA
    } else {
      fc <- rowMeans(CP100K_log[,is.element(colnames(CP100K_log), a)]) - rowMeans(CP100K_log[,!is.element(colnames(CP100K_log), a)])
      p <- apply(CP100K_log, 1, function(x) t.test(x[is.element(colnames(CP100K_log), a)], x[!is.element(colnames(CP100K_log), a)])$p.value)
      clusters_sig[[i]] <- data.frame("log2(FC)" = fc, "ttest_p" =p)
    }
  }
  
  return(list("clusters_cells" = clusters_final, "clusters_sig" = clusters_sig))
}




