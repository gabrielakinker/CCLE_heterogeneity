# ----------------------------------------------------------------------------------------------------
# Module 6. Evaluating the association between CNA subclones and expression heterogeneity in cell lines.
# ----------------------------------------------------------------------------------------------------

# **************************************************************************
# Basic setup

# load necessary R packages and functions
library(reshape2)
library(ggplot2)
source("nmf_cell_class.R")
source("robust_nmf_programs.R")

# read scRNA-seq data (CPM) from cell lines 
expr_ccle <- readRDS("CCLE_heterogeneity_Rfiles/CCLE_scRNAseq_CPM.RDS")     

# read inferred cna clones 
clones <- readRDS("Expected_results/module5/clone_assignment_ccle.RDS")

# read discrete programs 
dbscan_programs <- readRDS("Expected_results/module1/discr_clusters_minpt5_eps1.8_ccle.RDS")       

# read nmf programs (all - nmf ranks 6-9)
nmf_programs_genes <- readRDS("Expected_results/module1/nmf_w_basis_ccle.RDS") # nmf gene scores
nmf_programs_cells <- readRDS("Expected_results/module1/nmf_h_coef_ccle.RDS") # nmf cell scores

# **************************************************************************
# Select continuous nmf programs

# definy nmf cell program classification 
nmf_cell_scores_class <- unlist(unname(lapply(nmf_programs_cells, nmf_cell_class)), recursive = F)

# identify robust nmf programs and remove redanduncy
nmf_programs_sig  <- lapply(nmf_programs_genes, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
nmf_robust <- robust_nmf_programs(nmf_programs_sig, intra_min = 35, intra_max = 10, inter_filter=F)

# select only robust programs from nmf cell scores/class
nmf_programs_cells <-  unlist(unname(lapply(nmf_programs_cells, function(x) apply(x, 2, as.list))), recursive = F)
nmf_programs_cells <- lapply(nmf_programs_cells, unlist)
nmf_programs_cells <- nmf_programs_cells[nmf_robust]

nmf_cell_scores_class <- nmf_cell_scores_class[nmf_robust]

# indentify discrete nmf programs by comparing nmf and dbscan cell classification
nmf_vs_dbscan <- list()

for(i in names(nmf_cell_scores_class)) {
  if(!is.element(gsub(".{4}$", "", i), names(dbscan_programs))) next
  a <- colnames(expr_ccle[[gsub(".{4}$", "", i)]]) 
  b <- as.character(is.element(a,nmf_cell_scores_class[[i]])) 
  c <- dbscan_programs[[gsub(".{4}$", "", i)]][["clusters_cells"]]
  d <- c() 
  for(j in names(c)) {
    e <- table(b, as.character(is.element(a,c[[j]])))
    d[j] <- fisher.test(e)$p.value
  }
  nmf_vs_dbscan[[i]] <- d
}

nmf_vs_dbscan <- nmf_vs_dbscan[sapply(nmf_vs_dbscan, function(x) length(which(x<0.001)) != 0)]                                       
nmf_programs_cells <- nmf_programs_cells[!is.element(names(nmf_programs_cells), names(nmf_vs_dbscan))]                                     
                                               
# **************************************************************************
# Association between continuous programs and cna clones
                                        
# compare nmf and cna cell assigments
cna_conti <- list()

for(i in names(nmf_programs_cells)[is.element(gsub(".{4}$", "", names(nmf_programs_cells)), names(clones))]) {
  a <- nmf_programs_cells[[i]] # gets nmf scores for program i  
  b <- clones[[gsub(".{4}$", "", i)]] # gets clones
  a <- a[b$cell]  
  if(length(unique(b$final_class))==2) { # two clones
   cna_conti[[i]] <- t.test(a~b[,"final_class"])$p.value  
  }
  if(length(unique(b$final_class))>2) { # more than two clones
   cna_conti[[i]] <- summary(aov(a~b[,"final_class"]))[[1]][1,5]     
  }
}

cna_conti <- cna_conti[sapply(cna_conti, function(x) x<0.001)]

# **************************************************************************
# Association between discrete programs and cna clones
                                                                      
# unlist discret programs                              
dbscan_programs <- lapply(dbscan_programs, function(x) x[["clusters_cells"]])
dbscan_programs <- unlist(dbscan_programs, recursive = F)                       
                          
cna_discr <- list()

for(i in names(dbscan_programs)[is.element(gsub("\\..*", "", names(dbscan_programs)), names(clones))]) {
  a <- colnames(expr_ccle[[sub("\\..*", "", i)]]) # gets all cells from the selected cell line
  b <- as.character(is.element(a,dbscan_programs[[i]])) # checks which cells are part of the discrete cluster
  c <- clones[[sub("\\..*", "", i)]] # gets gmm test for the selected cell line for each chromosome arm
  e <- data.frame("clone" = c[,"final_class"], "discr" = b[match(c[,"cell"], a)]) # combines gmm test results and discrete program classification
  cna_discr[[i]] <- chisq.test(table(e))$p.value # fisher test - classification vs discrete program classification
}
                         
cna_discr <- cna_discr[sapply(cna_discr, function(x) x<0.001)]


# **************************************************************************
# Plot results              
                              
# pie chart
plot_cna <- data.frame("Discrete\nprograms" = c(length(dbscan_programs)-length(cna_discr), length(cna_discr)), "Continuous\nprograms" = c(length(nmf_cell_scores)-length(cna_conti), length(cna_conti)), check.names = F)
rownames(plot_cna) <- c("Not linked to CNA sub-clone", "Linked to CNA sub-clone")
plot_cna <- round(t(100*t(plot_cna)/colSums(plot_cna)),0)
plot_cna <- melt(as.matrix(plot_cna))

pdf("Output/module6/subclones_vs_expr_heterog.pdf", width = 6, height = 5)                    
ggplot(plot_cna, aes(x="", y=value, fill=Var1)) +
  geom_bar(stat="identity", alpha=0.8) +
  coord_polar("y") +
  scale_fill_manual(values=c("gray88", "maroon"), name="") +
  labs(x="", y="") +
  facet_wrap(~ Var2)  +
  theme(panel.background=element_blank() ,axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), strip.background = element_blank(), legend.position = "bottom", legend.justification = "center", strip.text.x=element_text(margin=margin(b=5)))
dev.off()

                              