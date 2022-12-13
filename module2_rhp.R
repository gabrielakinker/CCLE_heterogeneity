# ----------------------------------------------------------------------------------------------------
# Module 2. Definying heterogeneity patterns that are shared between multiple cell lines  and between multiple tumors
# (i.e. recurrent heterogeneous programs, RHPs)
# ----------------------------------------------------------------------------------------------------

# **************************************************************************
# Basic setup

# load necessary R packages and functions
library(reshape2)
library(ggplot2)
library(scales)
library(egg)
source("custom_magma.R")
source("cancer_color_ccle.R")
source("dbscan_get_sig.R")
source("robust_nmf_programs.R")
source("nmf_cell_class.R")


# read scRNA-seq data from cell lines and tumors
expr_ccle <- readRDS("CCLE_heterogeneity_Rfiles/CCLE_scRNAseq_CPM.RDS") # CCLE cell lines
expr_tumor <- readRDS("CCLE_heterogeneity_Rfiles/tumors_scRNAseq_logTPM.RDS") # human tumors

# read metadata
meta_ccle <- readRDS("CCLE_heterogeneity_Rfiles/CCLE_metadata.RDS")

# read heterogeneity programs identified using tSNE and DBSCAN (i.e. discrete heterogeneity)
dbscan_programs_ccle <- readRDS("Expected_results/module1/discr_clusters_minpt5_eps1.8_ccle.RDS")

# read heterogeneity programs identified using NMF (i.e. discrete + continuous heterogeneity)
nmf_programs_genes_ccle <- readRDS("Expected_results/module1/nmf_w_basis_ccle.RDS") # nmf gene scores
nmf_programs_cells_ccle <- readRDS("Expected_results/module1/nmf_h_coef_ccle.RDS") # nmf cell scores

nmf_programs_genes_tumor <- readRDS("Expected_results/module1/nmf_w_basis_tumor.RDS") # nmf gene scores
nmf_programs_cells_tumor <- readRDS("Expected_results/module1/nmf_h_coef_tumor.RDS") # nmf cell scores


# ************************************************************************** 
# Identifying recurrent heterogeneous programs (RHPs) in cell lines - tSNE + DBSCAN

# get gene programs for dbscan clusters (ignore cluster composed of more than 90% of cells)
dbscan_programs_sig_ccle <- lapply(names(dbscan_programs_ccle), function(x) dbscan_get_sig(dbscan_programs_ccle[[x]], cell_line=x, max_size=0.9))
names(dbscan_programs_sig_ccle) <- names(dbscan_programs_ccle)
dbscan_programs_sig_ccle <- unlist(dbscan_programs_sig_ccle, recursive = F)

# manually remove 2 redundant clusters                              
intersect(dbscan_programs_sig_ccle$`RERFLCAI_LUNG.2`, dbscan_programs_sig_ccle$`RERFLCAI_LUNG.4`)
intersect(dbscan_programs_sig_ccle$`RERFLCAI_LUNG.1`, dbscan_programs_sig_ccle$`RERFLCAI_LUNG.4`)
dbscan_programs_sig_ccle[c("RERFLCAI_LUNG.1", "RERFLCAI_LUNG.2")] <- NULL
                              
# calculate similarity between programs
dbscan_intersect_ccle <- sapply(dbscan_programs_sig_ccle, function(x) sapply(dbscan_programs_sig_ccle,function(y) length(intersect(y,x))/length(union(x,y))))*100
 
# hierarchical clustering of the similarity matrix                                                                   
dbscan_intersect_hc_ccle <- hclust(as.dist(100-dbscan_intersect_ccle), method="average") 
dbscan_intersect_hc_ccle <- reorder(as.dendrogram(dbscan_intersect_hc_ccle), colMeans(dbscan_intersect_ccle))
dbscan_intersect_ccle <- dbscan_intersect_ccle[order.dendrogram(dbscan_intersect_hc_ccle), order.dendrogram(dbscan_intersect_hc_ccle)]
                                                                   
# plot similarity matrix heatmap
dbscan_intersect_melt_ccle <- reshape2::melt(dbscan_intersect_ccle) 
                                                                   
p1.1 <- ggplot(data = dbscan_intersect_melt_ccle, aes(x=Var1, y=Var2, fill=value, color=value)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  scale_color_gradient2(limits=c(2, 25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  theme( axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.border = element_rect(fill=F),  axis.line = element_blank(), axis.text = element_text(size = 13), axis.title = element_text(size = 15), legend.title = element_text(size=12), legend.text = element_text(size = 11), legend.text.align = 0.5, legend.direction = "horizontal", legend.position = c(1.20,-0.1), plot.margin = unit(c(2,6,2,2), "cm")) + 
  scale_x_discrete(name="Discrete programs", breaks=unique(dbscan_intersect_melt_ccle$Var1)[seq(5, length(unique(dbscan_intersect_melt_ccle$Var1)), by=5)], labels= seq(5, length(unique(dbscan_intersect_melt_ccle$Var1)), by=5)) + 
  scale_y_discrete(name="", breaks=unique(dbscan_intersect_melt_ccle$Var1)[seq(5, length(unique(dbscan_intersect_melt_ccle$Var1)), by=5)], labels= seq(5, length(unique(dbscan_intersect_melt_ccle$Var1)), by=5)) +
  guides(fill = guide_colourbar(barheight = 0.8, barwidth = 4.8, title.position = "top", title.hjust = 0.5))                                                                   
                       
# plot similarity matrix cancer type annotation
dbscan_annot_type_ccle <- data.frame("cell_lines" =sub("\\..*", "", rownames(dbscan_intersect_ccle)), stringsAsFactors = F)
dbscan_annot_type_ccle$type <- meta_ccle$cancer_type_trunc[match(dbscan_annot_type_ccle$cell_lines, rownames(meta_ccle))]
dbscan_annot_type_ccle$type <- factor(dbscan_annot_type_ccle$type, levels=c(sort(unique(dbscan_annot_type_ccle$type))[sort(unique(dbscan_annot_type_ccle$type))!="Others"], "Others"))
                                                                   
p1.2 <- ggplot(dbscan_annot_type_ccle, aes(y="", x=1:nrow(dbscan_annot_type_ccle), fill=type, color=type)) +
  geom_tile() +
  scale_color_manual(values = cancer_color[match(levels(dbscan_annot_type_ccle$type), cancer_color$type), "color"] , name="") +
  scale_fill_manual(values = cancer_color[match(levels(dbscan_annot_type_ccle$type), cancer_color$type), "color"] , name="") +
   scale_x_continuous(expand = c(0,0)) +                                                           
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 11),  axis.title = element_text(size=8),  plot.margin = unit(c(2,4,-2,1), "cm"), legend.text.align = 0, legend.key.size = unit(0.4, "cm"), legend.key=element_blank()) +                                  
  labs(x="", y="") +
  guides(fill = guide_legend(title = "",  override.aes = list(size = 5), ncol=1), color=FALSE) 

# combine p1 plots                                                                   
pdf("Output/module2/dbscan_metaprograms_eps1.8_ccle.pdf", height = 6.5, width = 9, onefile = F)
egg::ggarrange(p1.2,p1.1, nrow=2, heights = c(0.5,9))
dev.off()
 
# manually define metaprograms (genes observed in at least 25% of programs composing the respective metacluster) 
dbscan_meta1_ccle <- sort(table(unlist(dbscan_programs_sig_ccle[colnames(dbscan_intersect_ccle)[4:9]])), decreasing=T)
dbscan_meta1_ccle[dbscan_meta1_ccle/length(4:9) > 0.25]
dbscan_meta2_ccle <- sort(table(unlist(dbscan_programs_sig_ccle[colnames(dbscan_intersect_ccle)[12:18]])), decreasing=T)
dbscan_meta2_ccle[dbscan_meta2_ccle/length(12:18) > 0.25]                                                                                                                         
# ************************************************************************** 
# Identifying recurrent heterogeneous programs (RHPs) in cell lines - NMF 
                                                                   
# get gene programs (top 50 genes by NMF score)
nmf_programs_sig_ccle <- lapply(nmf_programs_genes_ccle, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
                                             
# for each cell line, select robust NMF programs (i.e. obseved using different ranks in the same cell line), remove redundancy due to multiple ranks, and apply a filter based on the similarity to programs from other cell lines. 
nmf_filter_ccle <- robust_nmf_programs(nmf_programs_sig_ccle, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10)
nmf_programs_sig_ccle <- lapply(nmf_programs_sig_ccle, function(x) x[, is.element(colnames(x), nmf_filter_ccle),drop=F])
nmf_programs_sig_ccle <- do.call(cbind, nmf_programs_sig_ccle)
                            
# calculate similarity between programs
nmf_intersect_ccle <- apply(nmf_programs_sig_ccle , 2, function(x) apply(nmf_programs_sig_ccle , 2, function(y) length(intersect(x,y)))) 
                                                                
# hierarchical clustering of the similarity matrix 
nmf_intersect_hc_ccle <- hclust(as.dist(50-nmf_intersect_ccle), method="average") 
nmf_intersect_hc_ccle <- reorder(as.dendrogram(nmf_intersect_hc_ccle), colMeans(nmf_intersect_ccle))
nmf_intersect_ccle <- nmf_intersect_ccle[order.dendrogram(nmf_intersect_hc_ccle), order.dendrogram(nmf_intersect_hc_ccle)]

# plot similarity matrix heatmap     
nmf_intersect_meltI_ccle <- reshape2::melt(nmf_intersect_ccle) 

p2.1 <- ggplot(data = nmf_intersect_meltI_ccle, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  scale_x_discrete(name="\nPrograms", breaks=unique(nmf_intersect_meltI_ccle$Var1)[seq(100, length(unique(nmf_intersect_meltI_ccle$Var1)), by=100)], labels= seq(100, length(unique(nmf_intersect_meltI_ccle$Var1)), by=100)) + 
  scale_y_discrete(name="\nPrograms", breaks=unique(nmf_intersect_meltI_ccle$Var2)[seq(100, length(unique(nmf_intersect_meltI_ccle$Var2)), by=100)], labels= seq(100, length(unique(nmf_intersect_meltI_ccle$Var2)), by=100))  +
  geom_vline(xintercept=c(424,504), linetype="longdash", size=0.6)+
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))


# for each program, plot correlation between nmf cell scores and cell complexity
complexity_ccle <- lapply(expr_ccle, function(x) apply(x, 2, function(y) length(which(y != 0))))                                                               
nmf_corr_complexity_ccle <- data.frame(matrix(ncol = 1, nrow = ncol(nmf_intersect_ccle)), row.names=colnames(nmf_intersect_ccle))
colnames(nmf_corr_complexity_ccle) <- "corr"                                             
for(i in rownames(nmf_corr_complexity_ccle)) {
    a <- gsub(".{4}$", "", i)
    b <- nmf_programs_cells_ccle[[a]][,i]
    c <- complexity_ccle[[a]]
    nmf_corr_complexity_ccle[i,1] <- cor(b,c)
}

p2.2 <- ggplot(nmf_corr_complexity_ccle, aes(y=corr, x = 1:nrow(nmf_corr_complexity_ccle))) +
  geom_smooth(span=0.05, se = FALSE, color="gray36", size= 0.8, method = "loess") + 
  geom_point(alpha=0.3, size = 1.5) + 
  theme(axis.ticks = element_blank(), axis.line = element_blank(), panel.background = element_rect(fill = "gray97"), panel.grid = element_blank(), axis.text.x = element_blank(), axis.title = element_text(size = 12), axis.text.y = element_text(size = 10), plot.margin=unit(c(1,1,-0.6,1), "cm") ) + 
  labs(y="Cell complexity\ncorrelation", x="") + 
  scale_y_continuous(expand = c(0,0), limits = c(-0.9,0.9), breaks = seq(-0.8, 0.8, 0.4)) + 
  scale_x_continuous(expand = c(0,0))

# combine p2 plots
pdf("Output/module2/nmf_metaprograms_complexity_ccle.pdf", height = 6, width = 7, onefile=FALSE) 
egg::ggarrange(p2.2,p2.1, nrow=2, ncol = 1, heights = c(1,5))
dev.off()
                                           
# manually remove programs associated with cell complexity and save output
nmf_intersect_ccle <- nmf_intersect_ccle[-c(424:504, 509:514),-c(424:504, 509:514)]
saveRDS(nmf_programs_sig_ccle[,colnames(nmf_intersect_ccle)], "Output/module2/nmf_programs_sig_ccle.RDS")                                                                                          
# hierarchical clustering of the filtered similarity matrix 
nmf_intersect_hc_ccle <- hclust(as.dist(50-nmf_intersect_ccle), method="average") 
nmf_intersect_hc_ccle <- reorder(as.dendrogram(nmf_intersect_hc_ccle), colMeans(nmf_intersect_ccle))
nmf_intersect_ccle <- nmf_intersect_ccle[order.dendrogram(nmf_intersect_hc_ccle), order.dendrogram(nmf_intersect_hc_ccle)]

# plot similarity matrix heatmap     
nmf_intersect_meltII_ccle <- reshape2::melt(nmf_intersect_ccle) 

p3.1 <- ggplot(data = nmf_intersect_meltII_ccle, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  scale_x_discrete(name="\nPrograms", breaks=unique(nmf_intersect_meltII_ccle$Var1)[seq(100, length(unique(nmf_intersect_meltII_ccle$Var1)), by=100)], labels= seq(100, length(unique(nmf_intersect_meltII_ccle$Var1)), by=100)) + 
  scale_y_discrete(name="\nPrograms", breaks=unique(nmf_intersect_meltII_ccle$Var2)[seq(100, length(unique(nmf_intersect_meltII_ccle$Var2)), by=100)], labels= seq(100, length(unique(nmf_intersect_meltII_ccle$Var2)), by=100)) +
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
                                             
# annotate nmf programs that are discrete
nmf_programs_cells_class_ccle <- unlist(unname(lapply(nmf_programs_cells_ccle, nmf_cell_class)), recursive = F)

nmf_vs_dbscan <- list()
for(i in names(nmf_programs_cells_class_ccle)) {
  #if(!is.element(gsub(".{4}$", "", i), names(dbscan_programs_ccle))) dbscan_vs_nmf[[i]] <- NA
  a <- colnames(expr_ccle[[gsub(".{4}$", "", i)]]) 
  b <- as.character(is.element(a,nmf_programs_cells_class_ccle[[i]])) 
  c <- dbscan_programs_ccle[[gsub(".{4}$", "", i)]][["clusters_cells"]]
  d <- c() 
  for(j in names(c)) {
    e <- table(b, as.character(is.element(a,c[[j]])))
    d[j] <- fisher.test(e)$p.value
  }
  nmf_vs_dbscan[[i]] <- d
}


nmf_vs_dbscan <- nmf_vs_dbscan[sapply(nmf_vs_dbscan, function(x) length(which(x<0.001)) != 0)]
nmf_vs_dbscan <- data.frame("discrete"=is.element(colnames(nmf_intersect_ccle), names(nmf_vs_dbscan)))                                                       
p3.2 <- ggplot(nmf_vs_dbscan, aes(x=1:nrow(nmf_vs_dbscan_ccle), y=0, fill=discrete, colour=discrete))+
  geom_tile() + 
  scale_fill_manual(values = c("white", "black"), breaks="TRUE", labels="Discrete programs", name="") +
  scale_colour_manual(values = c("white", "black"), breaks="TRUE", labels="Discrete programs", name="") +
  theme(axis.ticks = element_blank(), panel.background = element_blank(),panel.border=element_rect(fill=F), axis.line = element_blank(), axis.text = element_blank(), plot.margin = unit(c(1,1,0,1), "cm")) + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +                                           
  labs(x="", y="") 


# combine p3 plots
pdf("Output/module2/nmf_metaprograms_discrete_ccle.pdf", height = 6, width = 7, onefile=FALSE)
egg::ggarrange(p3.2,p3.1, nrow=2, ncol = 1, heights = c(0.4,6))
dev.off()

# manually select non-cell cycle programs   
nmf_intersect_nc_ccle  <- nmf_intersect_ccle[-c(433:800), -c(433:800)]
                                                       
saveRDS(nmf_programs_sig_ccle[,colnames(nmf_intersect_nc_ccle)], "Output/module2/nmf_programs_sig_nc_ccle.RDS")           
                                                       
# hierarchical clustering of the filtered similarity matrix      
nmf_intersect_hc_nc_ccle <- hclust(as.dist(50-nmf_intersect_nc_ccle), method="average") 
nmf_intersect_hc_nc_ccle <- reorder(as.dendrogram(nmf_intersect_hc_nc_ccle), colMeans(nmf_intersect_nc_ccle))
nmf_intersect_nc_ccle <-  nmf_intersect_nc_ccle[order.dendrogram(nmf_intersect_hc_nc_ccle), order.dendrogram(nmf_intersect_hc_nc_ccle)]
  

# plot similarity matrix heatmap     
nmf_intersect_melt_nc_ccle <- melt(nmf_intersect_nc_ccle)

p4.1 <- ggplot(data = nmf_intersect_melt_nc_ccle, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() +
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")+  
  theme( axis.ticks = element_blank(), panel.background = element_blank(), panel.border=element_rect(colour = "gray40", size = 0.4, fill=F),  axis.line = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10),  legend.title = element_text(size=10), legend.text = element_text(size = 9), legend.text.align = 0.5, legend.direction = "horizontal", legend.justification="bottom", plot.margin = unit(c(0.5,3,0.5,0.5), "cm")) + 
  scale_y_discrete(name="\nPrograms", breaks=unique(nmf_intersect_melt_nc_ccle$Var2)[seq(100, length(unique(nmf_intersect_melt_nc_ccle$Var2)), by=100)], labels= seq(100, length(unique(nmf_intersect_melt_nc_ccle$Var2)), by=100)) +
   guides(color = guide_colourbar(barheight = 0.6, barwidth = 4.4, title.position = "top", title.hjust = 0.5)) 
                                                       
# plot similarity matrix cancer type annotation
nmf_annot_type_ccle <- data.frame("cell_lines" =gsub(".{4}$", "",rownames(nmf_intersect_nc_ccle)), stringsAsFactors = F)
nmf_annot_type_ccle$type <- meta_ccle$cancer_type_trunc[match(nmf_annot_type_ccle$cell_lines, rownames(meta_ccle))]
nmf_annot_type_ccle$type <- factor(nmf_annot_type_ccle$type, levels=c(sort(unique(nmf_annot_type_ccle$type))[sort(unique(nmf_annot_type_ccle$type))!="Others"], "Others"))
                                                                   
p4.2 <- ggplot(nmf_annot_type_ccle, aes(y="", x=1:nrow(nmf_annot_type_ccle), fill=type, color=type)) +
  geom_tile() +
  scale_color_manual(values = cancer_color[match(levels(nmf_annot_type_ccle$type), cancer_color$type), "color"] , name="") +
  scale_fill_manual(values = cancer_color[match(levels(nmf_annot_type_ccle$type), cancer_color$type), "color"] , name="") +                                                    
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 10),  axis.title = element_text(size=8),  legend.text.align = 0, legend.key.size = unit(0.4, "cm"), legend.key=element_blank(), legend.position=c(1.3,-4), plot.margin = unit(c(0.5,3,-0.6,0.5), "cm")) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x="", y="") +
  guides(fill = guide_legend(title = "", override.aes = list(size = 5), ncol = 1), color=FALSE) 
                         
                                                       
# manually define metaprograms (genes observed in at least 25% of programs composing the respective metacluster) 
nmf_meta1_ccle <- sort(table(nmf_programs_sig_ccle[,colnames(nmf_intersect_nc_ccle)[10:25]])/length(10:25), decreasing=T)
nmf_meta1_ccle_programs <- colnames(nmf_intersect_nc_ccle)[10:25]
                                                 
nmf_meta2_ccle <- sort(table(nmf_programs_sig_ccle[,colnames(nmf_intersect_nc_ccle)[46:53]])/length(46:53), decreasing=T) 
nmf_meta2_ccle_programs <- colnames(nmf_intersect_nc_ccle)[46:53]
                                                            
nmf_meta3_ccle <- sort(table(nmf_programs_sig_ccle[,colnames(nmf_intersect_nc_ccle)[57:75]])/length(57:75), decreasing=T)
nmf_meta3_ccle_programs <- colnames(nmf_intersect_nc_ccle)[57:75]                                     
nmf_meta4_ccle <- sort(table(nmf_programs_sig_ccle[,colnames(nmf_intersect_nc_ccle)[81:95]])/length(81:95), decreasing=T)
nmf_meta4_ccle_programs <- colnames(nmf_intersect_nc_ccle)[81:95]                    

nmf_meta5_ccle <- sort(table(nmf_programs_sig_ccle[,colnames(nmf_intersect_nc_ccle)[98:109]])/length(98:109), decreasing=T)
nmf_meta5_ccle_programs <- colnames(nmf_intersect_nc_ccle)[98:109]                                     
nmf_meta6_ccle <- sort(table(nmf_programs_sig_ccle[,colnames(nmf_intersect_nc_ccle)[116:132]])/length(116:132), decreasing = T)
nmf_meta6_ccle_programs <- colnames(nmf_intersect_nc_ccle)[116:132]
                                                           
nmf_meta7_ccle <- sort(table(nmf_programs_sig_ccle[,colnames(nmf_intersect_nc_ccle)[140:167]])/length(140:167), decreasing=T)
nmf_meta7_ccle_programs <- colnames(nmf_intersect_nc_ccle)[140:167]
                                               
nmf_meta8_ccle <- sort(table(nmf_programs_sig_ccle[,colnames(nmf_intersect_nc_ccle)[197:250]])/length(197:250), decreasing=T)
nmf_meta8_ccle_programs<-colnames(nmf_intersect_nc_ccle)[197:250]

nmf_meta9_ccle <- sort(table(nmf_programs_sig_ccle[,colnames(nmf_intersect_nc_ccle)[261:278]])/length(261:278), decreasing=T)
nmf_meta9_ccle_programs <- colnames(nmf_intersect_nc_ccle)[261:278]

nmf_meta10_ccle <- sort(table(nmf_programs_sig_ccle[,colnames(nmf_intersect_nc_ccle)[278:432]])/length(278:432), decreasing=T)
nmf_meta10_ccle_programs <- colnames(nmf_intersect_nc_ccle)[278:432]

nmf_meta_all_ccle <- list(skinpig=nmf_meta1_ccle, emtI=nmf_meta2_ccle, emtII=nmf_meta3_ccle, ifn=nmf_meta4_ccle, emtIII=nmf_meta5_ccle, p53snc=nmf_meta6_ccle, episen=nmf_meta7_ccle, stress=nmf_meta8_ccle, protmat=nmf_meta9_ccle, protdreg=nmf_meta10_ccle)       
                                      
nmf_meta_all_ccle_top25 <- lapply(nmf_meta_all_ccle, function(x) x[x>=0.25]) 
                                  
saveRDS(nmf_meta_all_ccle,"Output/module2/nmf_metaprograms_sig_nc_ccle.RDS")
                                     
saveRDS(nmf_meta_all_ccle_top25,"Output/module2/nmf_metaprograms_sigtop25_nc_ccle.RDS") 
                                      
saveRDS(list(skinpig=nmf_meta1_ccle_programs, emtI=nmf_meta2_ccle_programs, emtII=nmf_meta3_ccle_programs, ifn=nmf_meta4_ccle_programs, emtIII=nmf_meta5_ccle_programs, p53snc=nmf_meta6_ccle_programs, episen=nmf_meta7_ccle_programs, stress=nmf_meta8_ccle_programs, protmat=nmf_meta9_ccle_programs, protdreg=nmf_meta10_ccle_programs),"Output/module2/nmf_metaprograms_programs_nc_ccle.RDS")                                                                         

# plot NMF gene scores
nmf_nc_genes <- do.call(c,lapply(nmf_meta_all_ccle_top25, function(x) rev(names(x))))
nmf_nc_genes <- sapply(unlist(unname(lapply(nmf_programs_genes_ccle, function(x) apply(x, 2, list))), recursive=F)[colnames(nmf_intersect_nc_ccle)], function(x) unlist(x)[nmf_nc_genes])
rownames(nmf_nc_genes) <- NULL
nmf_nc_genes[is.na(nmf_nc_genes)] <- 0

nmf_nc_genes_melt <- melt(nmf_nc_genes)

p4.3 <- ggplot(data = nmf_nc_genes_melt, aes(x=Var2, y=Var1, fill=value*100)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(0,8), low= "#b5cce3", mid=  "white", high = "darkred", midpoint = 4,   oob=squish, name=expression(paste("NMF score (10"^"2", ")", sep = ""))) +
   scale_color_gradient2(limits=c(0,8), low= "#b5cce3", mid=  "white", high = "darkred", midpoint = 4,   oob=squish, name=expression(paste("NMF score (10"^"2", ")", sep = ""))) + 
  theme( axis.ticks = element_blank(), panel.background = element_blank(), panel.border=element_rect(colour = "black", size = 0.4, fill=F),   axis.line = element_blank(), axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.title = element_text(size=10), legend.text = element_text(size = 9), legend.text.align = 0.5, legend.direction = "horizontal",  plot.margin = unit(c(0,3,0.5,0.5), "cm")) +
  labs(x="\nPrograms", y="Genes") +
  scale_y_continuous(expand=c(0,0))  +                                        
  scale_x_discrete(breaks = colnames(nmf_intersect_nc_ccle)[seq(50, ncol(nmf_intersect_nc_ccle), 50)], labels = seq(50, ncol(nmf_intersect_nc_ccle), 50)) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 4.4, barheight = 0.6)) +
  geom_vline(xintercept=c(25, 53, 75, 95, 109, 132, 167, 250, 278), color = "black", size=0.2) 
                                                     
# combine plots 4
pdf("Output/module2/nmf_metaprograms_nc_ccle.pdf", height = 9.5, width = 8, onefile=FALSE)
egg::ggarrange( p4.2,p4.1,p4.3, nrow = 3, heights = c(0.5,9,7))
dev.off()                               
                                       
                                                

# ************************************************************************** 
# Identifying recurrent heterogeneous programs (RHPs) in tumors - NMF                                                              
# get gene programs (top 50 genes by NMF score)
nmf_programs_sig_tumor <- lapply(nmf_programs_genes_tumor, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
                                             
# for each cell line, select robust NMF programs (i.e. obseved using different ranks in the same cell line), remove redundancy due to multiple ranks, and apply a filter based on the similarity to programs from other cell lines. 
nmf_filter_tumor <- robust_nmf_programs(nmf_programs_sig_tumor, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10)
nmf_programs_sig_tumor <- lapply(nmf_programs_sig_tumor, function(x) x[, is.element(colnames(x), nmf_filter_tumor),drop=F])
nmf_programs_sig_tumor <- do.call(cbind, nmf_programs_sig_tumor)
                            
# calculate similarity between programs
nmf_intersect_tumor <- apply(nmf_programs_sig_tumor , 2, function(x) apply(nmf_programs_sig_tumor , 2, function(y) length(intersect(x,y)))) 
                                                                
# hierarchical clustering of the similarity matrix 
nmf_intersect_hc_tumor <- hclust(as.dist(50-nmf_intersect_tumor), method="average") 
nmf_intersect_hc_tumor <- reorder(as.dendrogram(nmf_intersect_hc_tumor), colMeans(nmf_intersect_tumor))
nmf_intersect_tumor <- nmf_intersect_tumor[order.dendrogram(nmf_intersect_hc_tumor), order.dendrogram(nmf_intersect_hc_tumor)]

# plot similarity matrix heatmap     
nmf_intersect_melt_tumor <- reshape2::melt(nmf_intersect_tumor) 

p5.1 <- ggplot(data = nmf_intersect_melt_tumor, aes(x=Var1, y=Var2, fill=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  scale_x_discrete(name="\nPrograms", breaks=unique(nmf_intersect_melt_tumor$Var1)[seq(10, length(unique(nmf_intersect_melt_tumor$Var1)), by=10)], labels= seq(10, length(unique(nmf_intersect_melt_tumor$Var1)), by=10)) + 
  scale_y_discrete(name="\nPrograms", breaks=unique(nmf_intersect_melt_tumor$Var2)[seq(10, length(unique(nmf_intersect_melt_tumor$Var2)), by=10)], labels= seq(10, length(unique(nmf_intersect_melt_tumor$Var2)), by=10))  +
  geom_vline(xintercept=c(424,504), linetype="longdash", size=0.6)+
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
                                                                           
                                                                           
# plot similarity matrix cancer type annotation
nmf_annot_type_tumor  <- data.frame("type"=sub("\\_.*", "", colnames(nmf_intersect_tumor)))

p5.2 <- ggplot(nmf_annot_type_tumor, aes(y=0, x=1:nrow(nmf_annot_type_tumor), fill=type, color=type)) +
  geom_tile() +
  scale_color_manual(values =c(brewer.pal(12, "Set3")[c(3,7)], "maroon", "slateblue2"),labels=c("Breast Cancer",  "Head and Neck Cancer", "Lung Cancer", "Melanoma"), name="") +
  scale_fill_manual(values =c(brewer.pal(12, "Set3")[c(3,7)], "maroon", "slateblue2"),labels=c("Breast Cancer",  "Head and Neck Cancer", "Lung Cancer", "Melanoma"), name="") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 12),   plot.margin = unit(c(0,4,-0.6,1), "cm"),  legend.box.margin = unit(c(4,0,1,1), "cm"), legend.position = "left", legend.text.align = 1, legend.justification = "right", legend.key.size = unit(0.4, "cm")) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x="", y="") +
  guides(fill = guide_legend(title = "", label.position="left",  override.aes = list(size = 5), ncol = 1), color=FALSE) 
                                                                           

# for each program, plot correlation between nmf cell scores and cell complexity            
complexity_tumor <- lapply(expr_tumor, function(x) apply(x, 2, function(y) length(which(y != 0))))                                                               
nmf_corr_complexity_tumor <- data.frame(matrix(ncol = 1, nrow = ncol(nmf_intersect_tumor)), row.names=colnames(nmf_intersect_tumor))
colnames(nmf_corr_complexity_tumor) <- "corr"                                             
for(i in rownames(nmf_corr_complexity_tumor)) {
    a <- gsub(".{4}$", "", i)
    b <- nmf_programs_cells_tumor[[a]][,i]
    c <- complexity_tumor[[a]]
    nmf_corr_complexity_tumor[i,1] <- cor(b,c)
}

p5.3 <- ggplot(nmf_corr_complexity_tumor, aes(y=corr, x = 1:nrow(nmf_corr_complexity_tumor))) +
  geom_smooth(span=0.05, se = FALSE, color="gray36", size= 0.8, method = "loess") + 
  geom_point(alpha=0.3, size = 1.5) + 
  theme(axis.ticks = element_blank(), axis.line = element_blank(), panel.background = element_rect(fill = "gray97"), panel.grid = element_blank(), axis.text.x = element_blank(), axis.title = element_text(size = 12), axis.text.y = element_text(size = 10), plot.margin=unit(c(1,1,-0.6,1), "cm") ) + 
  labs(y="Cell complexity\ncorrelation", x="") + 
  scale_y_continuous(expand = c(0,0), limits = c(-0.9,0.9), breaks = seq(-0.8, 0.8, 0.4)) + 
  scale_x_continuous(expand = c(0,0))                                                                           
# combine plots p5                                                                           
pdf("Output/module2/nmf_metaprograms_tumors.pdf", width = 9.5, height = 7, onefile = F)
ggarrange(p5.3,p5.2, p5.1, nrow = 3, heights = c(2,0.5,8)) 
dev.off()

# manually define metaprograms (genes observed in at least 25% of programs composing the respective metacluster) 
nmf_meta1_tumor <- sort(table(nmf_programs_sig_tumor[,colnames(nmf_intersect_tumor)[18:25]])/length(18:25), decreasing=T)
nmf_meta1_tumor_programs <- colnames(nmf_intersect_tumor)[18:25]                          
nmf_meta2_tumor <- sort(table(nmf_programs_sig_tumor[,colnames(nmf_intersect_tumor)[26:31]])/length(26:31), decreasing=T)
nmf_meta2_tumor_programs <- colnames(nmf_intersect_tumor)[26:31]
                           
nmf_meta3_tumor <- sort(table(nmf_programs_sig_tumor[,colnames(nmf_intersect_tumor)[32:37]])/length(32:37), decreasing=T)
nmf_meta3_tumor_programs <- colnames(nmf_intersect_tumor)[32:37]       
                                                                            
nmf_meta4_tumor <- sort(table(nmf_programs_sig_tumor[,colnames(nmf_intersect_tumor)[38:43]])/length(38:43), decreasing=T)
nmf_meta4_tumor_programs <- colnames(nmf_intersect_tumor)[38:43]                                                                            
nmf_meta5_tumor <- sort(table(nmf_programs_sig_tumor[,colnames(nmf_intersect_tumor)[44:52]])/length(44:52), decreasing=T)
nmf_meta5_tumor_programs <- colnames(nmf_intersect_tumor)[44:52]         
               
nmf_meta_all_tumor <- list(episen=nmf_meta1_tumor, immu=nmf_meta2_tumor, stress=nmf_meta3_tumor, hypoxia=nmf_meta4_tumor, emt=nmf_meta5_tumor)                                                                           
saveRDS(nmf_meta_all_tumor,"Output/module2/nmf_metaprograms_sig_nc_tumor.RDS")

saveRDS(lapply(nmf_meta_all_tumor, function(x) x[x>=0.25]),"Output/module2/nmf_metaprograms_sigtop25_nc_tumor.RDS")                                                         
                                      
saveRDS(list(episen=nmf_meta1_tumor_programs, immu=nmf_meta2_tumor_programs, stress=nmf_meta3_tumor_programs, hypoxia=nmf_meta4_tumor_programs, emt=nmf_meta5_tumor_programs),"Output/module2/nmf_metaprograms_programs_nc_tumor.RDS")                                                                           
