# --------------------------------------------------
# Function for vizualizing density-based clustering
# --------------------------------------------------

# - tpm = CPM expression matrix (rows = genes, columns = cells)
# - initial_dims = the number of dimensions that should be retained in the initial PCA step of the tSNE algorithm
# - perplexity = perplexity parameter to be used by the tSNE algorithm (should not be bigger than 3 * perplexity < ncol(cpm) - 1)
# - max_iter = number of iterations to be used by the tSNE algorithm
# - eps = size of the epsilon neighborhood to be used by the DBSCAN algorithm
# - minpts = number of minimum points in the eps region (for core points) to be used by the DBSCAN algorithm
# - min_size = minimum percentage of cells in each cluster
# - max_size = maximum percentage of cells in each cluster


# Returns a 2D tSNE plot with cells colored according to density-based clustering (DBSCAN)

library(Rtsne)
library(dbscan)
library(ggplot2)

dbscan_tsne_plot <- function(cpm, initial_dims=50, perplexity = 35, max_iter=7000, eps = 1.8, minpts=10, min_size= -Inf, max_size = Inf, seed=1, save.plot=FALSE, dir=NULL, plot.title=NULL)
 {
  
  ## Processing cpm data
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
  ds <- as.character(ds)
  ds[ds=="0"] <- NA
  
  # Filtering clusters based on size 
  size <- table(ds)
  ds[is.element(ds, names(size)[size > max_size*ncol(CP100K_log) | size < min_size*ncol(CP100K_log)])] <- NA
  
  # Plotting clusters
  tsne_plot <- data.frame("tSNE1"= tsne[,1], "tSNE2"=tsne[,2], "clusters"=ds)
  
  p1 <- ggplot(tsne_plot, aes(x = tSNE1, y = tSNE2, colour = clusters)) + 
    geom_point(alpha = 0.9, size=2.2) + 
    theme(panel.background = element_blank(), panel.border=element_rect(fill=F), axis.line = element_blank(),  axis.text = element_text(size=11), axis.title = element_text(size=12),  plot.title = element_text(size = 13), legend.key=element_blank()) +
    scale_color_brewer(palette = "Set2", name="Discrete\nClusters", na.value = "gray80", breaks=unique(na.omit(tsne_plot$clusters))) +
    labs(x="tSNE 1", y="tSNE 2", title=plot.title) 

  if(save.plot == TRUE) {  
    if(is.null(dir))  {
      pdf(paste0(plot.title, "tsne_clusters.pdf"), width = 3.8, height = 3.2)
      print(p1)
      dev.off()   
    } else {
        pdf(paste0(dir, "/", plot.title, "tsne_clusters.pdf"), width = 3.8, height = 3.2)
        print(p1)
        dev.off()   
    }
  }
  
  return(p1)
}











