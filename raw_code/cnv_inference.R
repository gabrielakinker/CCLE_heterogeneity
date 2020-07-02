###################################################################################################################################  CNV INFERENCE  - scRNA-seq CCLE ############################################
library(diptest)
library(mclust)
library(ggplot2)
library(scales)

################ processing data ##################


### reading gene locus annotation
gene_annotation <- read.table(file.choose(), sep = "\t", header = T, stringsAsFactors = F)
gene_annotation <- gene_annotation[is.element(gene_annotation$Chromosome.scaffold.name, c(1:22, "X")),]
gene_annotation$Chromosome.scaffold.name[gene_annotation$Chromosome.scaffold.name =="X"] <- 23
gene_annotation$Chromosome.scaffold.name <- as.numeric(gene_annotation$Chromosome.scaffold.name)
gene_annotation$chr.arm <- paste(gene_annotation$Chromosome.scaffold.name, substring(gene_annotation$Karyotype.band,1,1), sep = "")

### reading TPM expression data
expr_tpm <- readRDS("/home/labs/tirosh/kinker/analyses/CCLE/data processing/scCCLE&pool1_tpm_by_cell_line.rds") 

### reading CNV annotation from CCLE
cnv <- read.table(file.choose(), sep = "\t", header = T, stringsAsFactors = F)
rownames(cnv) <- cnv$SYMBOL
cnv <- cnv[,-(1:5)]
cnv <- cnv[,is.element(colnames(cnv), names(expr_tpm))]

### getting common genes 
common_genes <- Reduce(intersect, c(lapply(expr_tpm, rownames), list(gene_annotation$HGNC.symbol), list(rownames(cnv))))

### selecting common gene from matrix
gene_annotation <- gene_annotation[match(common_genes, gene_annotation$HGNC.symbol),]
expr_tpm <- lapply(expr_tpm, function(x) x[common_genes,])
cnv <- cnv[common_genes,]

### calculating average TPM values
ave_tpm <- rowMeans(sapply(expr_tpm, rowMeans))

### selecting top genes
expr <- lapply(expr_tpm, function(x) x[order(ave_tpm, decreasing = T)[1:7000],])

### log tranformation and library size correction 
expr <- lapply(expr, function(x) log2((x/10) + 1))
expr <- lapply(expr, function(x) t(t(x)-colMeans(x)))

### centering data gene-wise
ave_log_tpm <- rowMeans(sapply(expr, rowMeans))
expr_final <- lapply(expr, function(x) x - ave_log_tpm)


### truncating the data
for(i in names(expr_final)) {
  expr_final[[i]][expr_final[[i]] > 3] <- 3
  expr_final[[i]][expr_final[[i]] < -3] <- -3
}

### selecting genes in the annotation matrix
gene_annotation_sel <- gene_annotation[match(row.names(expr_final$NCIH2126_LUNG), gene_annotation$HGNC.symbol),]

### using expression data to calculate running average of 100 genes - genes ordered by chr localion
gene_order <- gene_annotation_sel$HGNC.symbol[order(gene_annotation_sel$Chromosome.scaffold.name, gene_annotation_sel$Karyotype.band)] # orders the genes by chromossome location

cnv_smooth_infer <- list()

for(i in names(expr_final)) {
  a <- t(expr_final[[i]]) 
  a <- a[,gene_order]
  b <- data.frame(matrix(ncol = ncol(a) - 99, nrow = nrow(a)), row.names = rownames(a))
  for(j in 1:ncol(b)) {
    b[,j] <- rowMeans(a[,j:(j+99)])
  }
  cnv_smooth_infer[[i]] <- b
}

  # save
  saveRDS(cnv_smooth_infer, "cnv_smooth_infer.RDS" )
  
### calculating average cnv (CCLE) per chromossome arm 
chr_arm <- split(as.character(gene_annotation_sel$HGNC.symbol), gene_annotation_sel$chr.arm) # divides genes by chromossome arm
cnv_by_arm <- sapply(chr_arm, function(y) if(length(y) > 1) colMeans(cnv[y,]) else cnv[y,]) # averages cnvs values by chromossome arm

### plotting average cnv (CCLE) per chromossome arm
cnv_by_arm_plot <- melt(cnv_by_arm)
cnv_by_arm_plot$status <- ""
cnv_by_arm_plot$status[cnv_by_arm_plot$value > -0.2 & cnv_by_arm_plot$value < 0.2] <- "Reference\ncell lines"
cnv_by_arm_plot <- cnv_by_arm_plot[order(as.numeric(gsub('.{1}$', '', cnv_by_arm_plot$Var2)), cnv_by_arm_plot$Var2),]
cnv_by_arm_plot$Var2 <- gsub("23", "X", cnv_by_arm_plot$Var2)

ggplot(cnv_by_arm_plot, aes(x=factor(Var2, levels = unique((Var2))), y=value, color=status)) +
  geom_point(size=1) +
  labs(x="Chromossome arm", y="CCLE CNV (log2)", color="") +
  scale_color_manual(values=c("gray70", "indianred2"), breaks = "Reference\ncell lines") +
  theme(axis.text.x = element_text(size=10, angle = 60, vjust = 0.8), axis.title = element_text(size=12), legend.text = element_text(size=12))
  
pdf("fig.2_ref_CNV.pdf", width = 10, height = 4)

### selecting model cell lines for each chromossome arm
model_by_arm <- apply(cnv_by_arm, 2, function(x) rownames(cnv_by_arm)[x > -0.2 & x < 0.2])
  
### determining the limits of each arm in the smoothed running average matrix
chr_arms_size <- table(gene_annotation_sel$chr.arm) # gets the number of genes in each chr arm 
chr_arms_size <- chr_arms_size[order(as.numeric(gsub('.{1}$', '', names(chr_arms_size))), names(chr_arms_size))] # orders by chr location

window_vs_arm <- c(rep(names(chr_arms_size)[1],chr_arms_size[1]-50)) # creates a annotation of chr arm for each of the 100 genes window

for(i in 2:length(chr_arms_size)) {
  window_vs_arm <- c(window_vs_arm, rep(names(chr_arms_size)[i], (chr_arms_size)[i]))
}

window_vs_arm <- window_vs_arm[1:6901]

### calculating average cnv (inferred) by  chromossome arm
expr_by_arm <- lapply(cnv_smooth_infer, function(x) data.frame(aggregate(t(x), by=list(window_vs_arm), mean), row.names = 1))

  # save
  saveRDS(expr_by_arm, "cnv_smooth_infer_by_arm.RDS" )

### testing the unimodality of cell scores by chromossome arm
# Gaussian finite mixture model fitted by EM algorithm 
gmm_test <- lapply(expr_by_arm, function(x) apply(x, 1, function(y) Mclust(y, G=1:2)))

gmm_final <- lapply(gmm_test, function(x) lapply(x, function(y) data.frame("class"=y$classification, "uncer" = y$uncertainty)))

for(i in names(gmm_final)) {
  for(j in names(gmm_final[[i]])) {
    a <- gmm_final[[i]][[j]]
    a <- a[a$uncer<0.01,]
    if(length(which(table(a$class) > 20)) <2) gmm_final[[i]][[j]] <- NULL
  }
}

gmm_final <- gmm_final[sapply(gmm_final, function(x) length(x)!=0)]
  
    # save 
    saveRDS(gmm_test, "gmm_test_complete.RDS")
    saveRDS(gmm_final, "gmm_test_filtered.RDS")

### creating a reference for the smoothed inferred cnv using the model cell lines for each arm
ref_cnv_infer_smooth <- list()

for(i in 1:ncol(cnv_smooth_infer$NCIH2126_LUNG)) {
  a <- model_by_arm[[window_vs_arm[i]]] # gets model cell lines for that window
  ref_cnv_infer_smooth[[i]] <- sapply(cnv_smooth_infer[a], function(x) mean(x[,i])) # averages the cnv inferrence of that window for each cell line
}

  # save
  saveRDS(ref_cnv_infer_smooth, "ref_cnv_smooth_infer.RDS" )
  

### for plotting: getting chromosome limits 
vlines_cnv <- table(gene_annotation_sel$Chromosome.scaffold.name[order(gene_annotation_sel$Chromosome.scaffold.name)])[1] - 50
for(i in 2:22) {
  vlines_cnv <- c(vlines_cnv, table(gene_annotation_sel$Chromosome.scaffold.name[order(gene_annotation_sel$Chromosome.scaffold.name)])[i] + vlines_cnv[i-1])
}

  # save
  saveRDS(vlines_cnv, "vlines_cnv_plot.RDS" )

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Hartigans' dip test
  dip_test <- lapply(expr_by_arm, function(x) apply(x, 1, function(y) dip.test(y, simulate.p.value = T)$p.value))
  dip_test_final <- lapply(dip_test, function(x) names(x)[x<0.05])
  dip_test_final <- dip_test_final[sapply(dip_test_final, function(x) length(x) !=0)]
  
  # save
  saveRDS(dip_test_final, "CNVclones_dip_test.RDS" )
  
  
### plotting selected cell lines
test <- cnv_smooth_infer$SCC25_UPPER_AERODIGESTIVE_TRACT
rownames(test) <- paste("c", rownames(test))
min <- sapply(ref_cnv_infer_smooth, function(x) quantile(x,0.25))
max <- sapply(ref_cnv_infer_smooth, function(x) quantile(x,0.75))
mean <- sapply(ref_cnv_infer_smooth, mean)
sd <- sapply(ref_cnv_infer_smooth, sd)

test <- t(t(test)-mean)

test_final <- test

for(i in 1:ncol(test_final)) {
  a <- test_final[,i]
  a[test[,i] > max[i] + 0.2] <-  a[test[,i] > max[i] + 0.2] - max[i]
  a[test[,i] < min[i] - 0.2] <-  a[test[,i] < min[i] - 0.2] - min[i]
  a[test[,i] >= min[i] - 0.2 & test[,i] <= max[i] + 0.2] <- 0
  test_final[,i] <- a
}


cnv_ind_hc <- hclust(as.dist(1-cor(t(test))), method = "average")
cnv_ind_smooth_plot <- melt(as.matrix(test_final[order(expr_by_arm$SCC9_UPPER_AERODIGESTIVE_TRACT["1q",]),]))
cnv_ind_smooth_plot <- melt(as.matrix(test_final))


ggplot(data = cnv_ind_smooth_plot, aes(y=Var1, x=Var2)) + 
  geom_raster(aes(fill=value), interpolate = T) + 
  scale_fill_gradient2(limits=c(-0.8, 0.8), midpoint = 0, low=c("dodgerblue4"), mid= c("white"), high = c("darkred"),   oob=squish, name="Inferred CNV (log2)") +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=NA, color="black"), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text.y  = element_blank(), legend.title = element_text(size=11, vjust = 0.8), axis.title= element_text(size=12), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.position = "bottom", legend.justification = "right")+ 
  labs(x="\nChromosome", y="") +
  scale_x_discrete(breaks= colnames(test)[c(1, vlines_cnv_all)], labels= c(names(vlines_cnv_all), "X")) +
  geom_vline(xintercept = vlines_cnv_all[1]) +
  geom_vline(xintercept = vlines_cnv_all[2]) +
  geom_vline(xintercept = vlines_cnv_all[3]) +
  geom_vline(xintercept = vlines_cnv_all[4]) +
  geom_vline(xintercept = vlines_cnv_all[5]) +
  geom_vline(xintercept = vlines_cnv_all[6]) +
  geom_vline(xintercept = vlines_cnv_all[7]) +
  geom_vline(xintercept = vlines_cnv_all[8]) +
  geom_vline(xintercept = vlines_cnv_all[9]) +
  geom_vline(xintercept = vlines_cnv_all[10]) +
  geom_vline(xintercept = vlines_cnv_all[11]) +
  geom_vline(xintercept = vlines_cnv_all[12]) +
  geom_vline(xintercept = vlines_cnv_all[13]) +
  geom_vline(xintercept = vlines_cnv_all[14]) +
  geom_vline(xintercept = vlines_cnv_all[15]) +
  geom_vline(xintercept = vlines_cnv_all[16]) +
  geom_vline(xintercept = vlines_cnv_all[17]) +
  geom_vline(xintercept = vlines_cnv_all[18]) +
  geom_vline(xintercept = vlines_cnv_all[19]) +
  geom_vline(xintercept = vlines_cnv_all[20]) +
  geom_vline(xintercept = vlines_cnv_all[21]) +
  geom_vline(xintercept = vlines_cnv_all[22])









expr_by_arm_ref <- list()
for(i in names(chr_arm)) {
  a <- model_by_arm[[i]]
  b <- chr_arm[[i]]
  if(length(chr_arm[[i]]) > 1) expr_by_arm_ref[[i]] <-  sapply(expr_final[a], function(x) mean(x[b,]))
  else expr_by_arm_ref[[i]] <-  sapply(expr_final[a], function(x) x[b,])
}


### normalizing expr by arm
expr_by_arm_norm <- list()

for(i in names(expr_by_arm)) {
  a <- t( (t(expr_by_arm[[i]])-sapply(expr_by_arm_ref, mean))/sapply(expr_by_arm_ref, sd))
  a[a>-2 & a < 2] <- 0
  expr_by_arm_norm[[i]] <- a
}













gmm_test_sig <- list()
for(i in names(gmm_test)) {
  a <- list()
  for(j in names(gmm_test[[i]])) {
    if(gmm_test[[i]][[j]]$G ==1 | min(table(gmm_test[[i]][[j]]$classification)) <= 0.02*length(gmm_test[[i]][[j]]$classification) )  a[[j]] <- 1
    else a[[j]] <- mclustBootstrapLRT(t(expr_by_arm[[i]][j,]), modelName =gmm_test[[i]][[j]]$modelName, maxG = max(gmm_test[[i]][[j]]$G))
  }
  gmm_test_sig[[i]] <- a
} 



  
  lapply(expr_by_arm, function(x) apply(x, 2, function(y) mclustBootstrapLRT(y, modelNames ="E" )$p.value))

gmm_testI <- lapply(expr_by_arm, function(x) apply(x, 2, function(y) Mclust(y, modelNames ="V")))
gmm_test_sigI <- lapply(expr_by_arm, function(x) apply(x, 2, function(y) mclustBootstrapLRT(y, modelNames ="V" )$p.value))

gmm_test_final <-unlist(lapply(gmm_test, function(x) lapply(x, function(y) aggregate(y$data, list(y$classification), mean))), recursive = F)
gmm_test_final <- gmm_test_final[sapply(gmm_test_final, function(x) length(x$V1) > 1)]
gmm_test_final <- gmm_test_final[sapply(gmm_test_final, function(x) length(which(x$V1 > 0.25 | x$V1 < -0.25)) != 0)]

### 
median_expr_by_arm <- apply(sapply(expr_by_arm, function(x) apply(x, 2, mean)), 1, median)
expr_by_arm_cent <- lapply(expr_final, function(x) sapply(chr_arm, function(y) if(length(y) > 1) colMeans(x[y,]) else x[y,]))



