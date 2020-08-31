# ----------------------------------------------------------------------------------------------------
# Module 5. Inferring CNV subclones in cell lines using scRNA-seq data
# ----------------------------------------------------------------------------------------------------

# **************************************************************************
# Basic setup

# load necessary R packages and functions
library(mclust)

# read gene locus annotation
gene_locus <- readRDS("CCLE_heterogeneity_Rfiles/gene_locus.RDS")

# read scRNA-seq data (CPM) from cell lines 
expr_ccle <- readRDS("CCLE_heterogeneity_Rfiles/CCLE_scRNAseq_CPM.RDS")      
                     
# select common genes between datasets
common_genes <- Reduce(intersect, c(lapply(expr_ccle, rownames), list(gene_locus$HGNC.symbol)))

gene_locus <- gene_locus[match(common_genes, gene_locus$HGNC.symbol),]
expr_ccle <- lapply(expr_ccle, function(x) x[common_genes,])
                   
# select top expressed genes
ave_expr <- rowMeans(sapply(expr_ccle, rowMeans))                  
expr_ccle <- lapply(expr_ccle, function(x) x[order(ave_expr, decreasing = T)[1:7000],])     
                    
gene_locus <- gene_locus[match(row.names(expr_ccle$NCIH2126_LUNG), gene_locus$HGNC.symbol),]
                    
# process data   
expr_ccle <- lapply(expr_ccle, function(x) log2((x/10) + 1))
expr_ccle <- lapply(expr_ccle, function(x) t(t(x)-colMeans(x)))
ave_expr_log <- rowMeans(sapply(expr_ccle, rowMeans))
expr_ccle <- lapply(expr_ccle, function(x) x - ave_expr_log)

# **************************************************************************
# Infer large scale copy number aberrations in chromosome arms
                    
# truncate expression data
for(i in names(expr_ccle)) {
  expr_ccle[[i]][expr_ccle[[i]] > 3] <- 3
  expr_ccle[[i]][expr_ccle[[i]] < -3] <- -3
}

# genes by chromosome location                     
gene_order <- gene_locus$HGNC.symbol[order(gene_locus$Chromosome.scaffold.name, gene_locus$Karyotype.band)]        
                    
# calculate running average expression in windows of 100 genes 
cna_infer <- list()

for(i in names(expr_ccle)) {
  a <- expr_ccle[[i]][gene_order,]
  b <- data.frame(matrix(ncol = nrow(a) - 99, nrow = ncol(a)), row.names = colnames(a))
  for(j in 1:ncol(b)) {
    b[,j] <- colMeans(a[j:(j+99),])
  }
  cna_infer[[i]] <- b
}                    
                    
                    
# save
saveRDS(cna_infer, "Output/module5/cna_infer_ccle.RDS")  
                 
# determine the limits of each chromosome arm in the inferred cna matrix
chr_arms_size <- table(gene_locus$chr.arm) # get the number of genes in each chr arm 
chr_arms_size <- chr_arms_size[order(as.numeric(gsub('.{1}$', '', names(chr_arms_size))), names(chr_arms_size))] # order by chr location

window_vs_arm <- c(rep(names(chr_arms_size)[1],chr_arms_size[1]-50)) # annotate the chormosome arm of each window in the inferred cna matrix
for(i in 2:length(chr_arms_size)) {
  window_vs_arm <- c(window_vs_arm, rep(names(chr_arms_size)[i], (chr_arms_size)[i]))
}
window_vs_arm <- window_vs_arm[1:6901]
                                        
# average the inferred cna matrix by chromossome arm                                        
cna_infer_arm <- lapply(cna_infer, function(x) data.frame(aggregate(t(x), by=list(window_vs_arm), mean), row.names = 1))

# save
saveRDS(cna_infer_arm, "Output/module5/cna_infer_arm_ccle.RDS" )

# **************************************************************************
# Detection of arm-level CNA subclones                 
                        
# fit a bimodal Gaussian mixture (via EM algorithm) for each arm in each cell line
gmm_test <- lapply(cna_infer_arm, function(x) apply(x, 1, function(y) Mclust(y, G=1:2)))

# confindently define cell lines that have subclones (i.e. those with > 20 cells classified into a second mode with > 99% confidence)                                                   
cna_subclones <- lapply(gmm_test, function(x) lapply(x, function(y) data.frame("class"=y$classification, "uncer" = y$uncertainty)))

for(i in names(cna_subclones)) {
  for(j in names(cna_subclones[[i]])) {
    a <- cna_subclones[[i]][[j]]
    a <- a[a$uncer<0.01,]
    if(length(which(table(a$class) > 20)) <2) cna_subclones[[i]][[j]] <- NULL
  }
}

cna_subclones <- cna_subclones[sapply(cna_subclones, function(x) length(x)!=0)]
  
# save 
saveRDS(cna_subclones, "Output/module5/cna_subclones_ccle.RDS")                    
                    
# confindently assign cells to clones (only consider cells classified with >90% confidence)
clone_assignment <- list()     
                                      
for(i in names(cna_subclones)) {
  if(length(cna_subclones[[i]]) > 1) { # more than one clone - we considered all combinations of modes with at least 5 cells.
    a <- colnames(expr_ccle[[i]])
    b <- do.call(paste, lapply(cna_subclones[[i]], function(x) x[["class"]]))
    b <- data.frame("cell"=a, "class"=b, stringsAsFactors = F)
    c <- sapply(cna_subclones[[i]], function(x) x[["uncer"]])
    b <- b[apply(c,1,function(x) length(which(x<=0.1))==length(x)),]
    b <- b[is.element(b$class, names(which(table(b$class)>=5))),]
    d <- data.frame(names(sort(table(b$class), decreasing = T)), 1:length(unique(b$class)), stringsAsFactors = F)
    b$final_class <- d[,2][match(b$class, d[,1])]
    b$final_class <- paste(i, "_", b$final_class, sep="")
    clone_assignment[[i]] <- b
  }
  if(length(cna_subclones[[i]]) == 1) { # only one clone -  we defined two clones corresponding to the two modes
    a <- colnames(expr_ccle[[i]])
    b <- data.frame("cell"=a, "class"=cna_subclones[[i]][[1]][["class"]], stringsAsFactors = F)
    b <- b[cna_subclones[[i]][[1]][["uncer"]] <= 0.1,]
    b$final_class <- paste(i, "_", b$class, sep="")
    clone_assignment[[i]] <- b
  }
}

# save 
saveRDS(clone_assignment, "Output/module5/clone_assignment_ccle.RDS")                                        
