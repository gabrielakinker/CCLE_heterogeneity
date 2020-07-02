############################################################################################################################################################################################### DISCRETE PROGRAMS OF VARIABILITY ########################################################


######################## Prepearing the data ########################## 
expr <- readRDS("/home/labs/tirosh/kinker/analyses/CCLE/data processing/scCCLE&pool1_tpm_by_cell_line.rds") 
expr <- lapply(expr, function(x) log2((x/10) + 1)) # log transformation correcting for average library size
expr <- lapply(expr, function(x) x[apply(x, 1, function(y) length(which(y > 3.5)) > ncol(x)*0.02),]) # selects the top genes
expr <- lapply(expr, function(x) x-rowMeans(x)) # centers the data gene-wise


####################### Getting cell lines with more than one cluster #######################

# eps 1.8
discr_clusters_minpt5_1.8 <- lapply(expr, function(x) get_discrete_programs_tsne(expression_data =  x, distance_cutoff = 1.8, points_cutoff=5, min_size_cutoff= 0.02, max_size_cutoff = Inf,  perplexity=35, iteractions = 7000, seed=1))
discr_clusters_minpt5_1.8 <- discr_clusters_minpt5_1.8[sapply(discr_clusters_minpt5_1.8, is.list)]

  ## save output
  saveRDS(discr_clusters_minpt5_1.8, "discr_clusters_minpt5_eps1.8.RDS")
  
# eps 1.6
discr_clusters_minpt5_1.5 <- lapply(expr, function(x) get_discrete_programs_tsne(expression_data =  x, distance_cutoff = 1.5, points_cutoff=5, min_size_cutoff= 0.02, max_size_cutoff = Inf,  perplexity=35, iteractions = 7000, seed=1))
discr_clusters_minpt5_1.5 <- discr_clusters_minpt5_1.5[sapply(discr_clusters_minpt5_1.5, is.list)]
  
  ## save output
  saveRDS(discr_clusters_minpt5_1.5, "discr_clusters_minpt5_eps1.5.RDS")
  

# eps 1.2
discr_clusters_minpt5_1.2 <- lapply(expr, function(x) get_discrete_programs_tsne(expression_data =  x, distance_cutoff = 1.2, points_cutoff=5, min_size_cutoff= 0.02, max_size_cutoff = Inf,  perplexity=35, iteractions = 7000, seed=1))
discr_clusters_minpt5_1.2 <- discr_clusters_minpt5_1.2[sapply(discr_clusters_minpt5_1.2, is.list)]

  ## save output
  saveRDS(discr_clusters_minpt5_1.2, "discr_clusters_minpt5_eps1.2.RDS")

###################### Plotting cell lines with more than one cluster eps 1.8 ##########################
lapply(names(discr_clusters_minpt5_1.8), function(x) plot_tsne_clusters(expression_data =  expr[[x]], distance_cutoff = 1.8,  points_cutoff=5, min_size_cutoff= 0.02, max_size_cutoff=Inf, perplexity=35, iteractions = 7000, seed=1, save.plot = T, plot.title = paste(x, "_minpt5_", sep = "")))


################# Plotting cell lines with only one cluster eps 1.8 ############################
plot_tsne_clusters(expression_data =  expr$COLO741_SKIN, distance_cutoff = 1.8,  points_cutoff=5, min_size_cutoff= 0.02, max_size_cutoff=Inf, perplexity=35, iteractions = 7000, seed=1, save.plot = T, plot.title = paste("COLO741_SKIN", "_minpt5_", sep = ""))

plot_tsne_clusters(expression_data =  expr$JHU006_UPPER_AERODIGESTIVE_TRACT, distance_cutoff = 1.8,  points_cutoff=5, min_size_cutoff= 0.02, max_size_cutoff=Inf, perplexity=35, iteractions = 7000, seed=1, save.plot = T, plot.title = paste("JHU006_UPPER_AERODIGESTIVE_TRACT", "_minpt5_", sep = ""))

plot_tsne_clusters(expression_data =  expr$HCC38_BREAST, distance_cutoff = 1.8,  points_cutoff=5, min_size_cutoff= 0.02, max_size_cutoff=Inf, perplexity=35, iteractions = 7000, seed=1, save.plot = T, plot.title = paste("HCC38_BREAST", "_minpt5_", sep = ""))

plot_tsne_clusters(expression_data =  expr$MFE319_ENDOMETRIUM, distance_cutoff = 1.8,  points_cutoff=5, min_size_cutoff= 0.02, max_size_cutoff=Inf, perplexity=35, iteractions = 7000, seed=1, save.plot = T, plot.title = paste("MFE319_ENDOMETRIUM", "_minpt5_", sep = ""))
  
  
  








