library(ggplot2)
library(Matrix)
library(gplots)
library(RColorBrewer)
library(reshape2)
library(scales)
library(ggdendro)
library(ggrepel)
library(viridis)

########################## inputing data ###############################

### reading expression data per cell line - tpm
expr <- readRDS("~/CCLE/data processing/scCCLE&pool1_tpm_by_cell_line.rds") 
expr <- lapply(expr, function(x) log2((x/10) + 1))
common_genes <- Reduce(intersect, lapply(expr, rownames))
expr <- lapply(expr, function(x) x[common_genes,])
complexity <- lapply(expr, function(x) apply(x, 2, function(y) length(which(y != 0)))) # calculates cell complexity

### reading cell scores per cell line - programs nmf
h_coef_6 <- readRDS("h_coef_nmf6.rds")
for(i in names(h_coef_6)) {
  colnames(h_coef_6[[i]]) <- paste(i, "_6.", 1:6, sep = "")
}

h_coef_7 <- readRDS("h_coef_nmf7.rds")
for(i in names(h_coef_7)) {
  colnames(h_coef_7[[i]]) <- paste(i, "_7.", 1:7, sep = "")
}

h_coef_8 <- readRDS("h_coef_nmf8.rds")
for(i in names(h_coef_8)) {
  colnames(h_coef_8[[i]]) <- paste(i, "_8.", 1:8, sep = "")
}

h_coef_9 <- readRDS("h_coef_nmf9.rds")
for(i in names(h_coef_9)) {
  colnames(h_coef_9[[i]]) <- paste(i, "_9.", 1:9, sep = "")
}

### combining cell scores
cell_score_all <-  unlist(unname(lapply(c(h_coef_6, h_coef_7, h_coef_8, h_coef_9), function(x) apply(x, 2, as.list))), recursive = F)
cell_score_all <- lapply(cell_score_all, unlist)

### reading gene scores per cell line - programs nmf
w_basis_6 <- readRDS("w_basis_nmf6.rds")
for(i in names(w_basis_6)) {
  colnames(w_basis_6[[i]]) <- paste(i, "_6.", 1:6, sep = "")
}

w_basis_7 <- readRDS("w_basis_nmf7.rds")
for(i in names(w_basis_7)) {
  colnames(w_basis_7[[i]]) <- paste(i, "_7.", 1:7, sep = "")
}

w_basis_8 <- readRDS("w_basis_nmf8.rds")
for(i in names(w_basis_8)) {
  colnames(w_basis_8[[i]]) <- paste(i, "_8.", 1:8, sep = "")
}

w_basis_9 <- readRDS("w_basis_nmf9.rds")
for(i in names(w_basis_9)) {
  colnames(w_basis_9[[i]]) <- paste(i, "_9.", 1:9, sep = "")
}

### combining all gene scores
gene_score_all <-  unlist(unname(lapply(c(w_basis_6, w_basis_7, w_basis_8, w_basis_9), function(x) apply(x, 2, as.list))), recursive = F)
gene_score_all <- lapply(gene_score_all, unlist)

### combining nmf results (w matrix) by cell line
nmf_all <- Map(cbind, w_basis_6, w_basis_7, w_basis_8, w_basis_9)
names(nmf_all) <- names(w_basis_6)

### selecting the top 50 genes in each program
nmf_all <- lapply(nmf_all, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))

### selecting recurrent programs for each cell line
nmf_all_intersect <- lapply(nmf_all, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x,y))))) ## for each cell line, calculates the intersection between programs
nmf_all_intersect <- lapply(nmf_all_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = T)[2])) ## getting the maximum intersection value for each program

nmf_all <- lapply(names(nmf_all), function(x) nmf_all[[x]][,nmf_all_intersect[[x]]>=35]) ## selecting programs with a maximum intersection of at least 35
names(nmf_all) <- names(w_basis_6)

### final selection of programs
nmf_all_combined <- do.call(cbind, nmf_all)
nmf_all_combined_intersect <- apply(nmf_all_combined , 2, function(x) apply(nmf_all_combined , 2, function(y) length(intersect(x,y)))) ## calculating intersection between all programs


prefinal_filter <- NULL # removes redundancy 

for(i in names(nmf_all)) {
  a <- nmf_all_combined_intersect[grep(i, colnames(nmf_all_combined_intersect), invert = T),grep(i, colnames(nmf_all_combined_intersect))]
  b <- sort(apply(a, 2, max), decreasing = T) # for each cell line, ranks programs based on their maximum intersection with programs of other cell lines
  if(length(b) != 0) {
    c <- names(b[1]) 
    for(y in 2:length(b)) {
      if(max(nmf_all_combined_intersect[c,names(b[y])]) <= 10) c <- c(c,names(b[y])) # selects programs interactivly from top, down. Only selects programs that have a intersection smaller than 25 with previously selected programs
    }
    prefinal_filter <- c(prefinal_filter, c)
  }
}

 # save 
  saveRDS(prefinal_filter, "redundancy_filter_nmf.RDS")

final_filter <- NULL # removes redudandancy and only selects those programs with a minimum similarity with another program

for(i in names(nmf_all)) {
  a <- nmf_all_combined_intersect[grep(i, colnames(nmf_all_combined_intersect), invert = T),grep(i, colnames(nmf_all_combined_intersect))]
  b <- sort(apply(a, 2, max), decreasing = T) # for each cell line, ranks programs based on their maximum intersection with programs of other cell lines
  b <- b[b>9] # selects programs with a maximum intersection of at least 10
  if(length(b) != 0) {
    c <- names(b[1]) 
    for(y in 2:length(b)) {
      if(max(nmf_all_combined_intersect[c,names(b[y])]) <= 10) c <- c(c,names(b[y])) # selects programs interactivly from top, down. Only selects programs that have a intersection smaller than 10 with previously selected programs
    }
    final_filter <- c(final_filter, c)
  }
}

### clustering programs
intersect_programs_filt <- nmf_all_combined_intersect[final_filter,final_filter]
hc_intersect <- hclust(as.dist(50-intersect_programs_filt), method="average") 
hc_intersect <- reorder(as.dendrogram(hc_intersect), colMeans(intersect_programs_filt))
intersect_programs_hc <- intersect_programs_filt[order.dendrogram(hc_intersect), order.dendrogram(hc_intersect)]

### plotting metaprograms
gabi_palette <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

plot_hc <- reshape2::melt(intersect_programs_hc) 

p1 <- ggplot(data = plot_hc, aes(x=Var1, y=Var2, fill=100*value/(100-value))) + 
  geom_raster() + 
  scale_fill_gradient2(limits=c(2,25), low=gabi_palette[1:111],  mid =gabi_palette[112:222], high = gabi_palette[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  scale_x_discrete(name="\nPrograms", breaks=colnames(intersect_programs_hc)[seq(100, ncol(intersect_programs_hc), by=100)], labels= seq(100, ncol(intersect_programs_hc), by=100)) + 
  scale_y_discrete(name="", breaks=colnames(intersect_programs_hc)[seq(100, ncol(intersect_programs_hc), by=100)], labels= seq(100, ncol(intersect_programs_hc), by=100))  +
  geom_vline(xintercept=c(424,504), linetype="longdash", size=0.6)+
  panel_border(colour = "gray40", size = 0.4, linetype = 1, remove = FALSE) +
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))


### plotting correlation of cell scores (for each program) vs complexity
corr_complexity <- data.frame(matrix(ncol = 2, nrow = ncol(intersect_programs_hc)))
colnames(corr_complexity) <- c("programs", "corr")
corr_complexity[,1] <- colnames(intersect_programs_hc)

for(i in 1:ncol(intersect_programs_hc)) {
  a <- cell_score_all[[colnames(intersect_programs_hc)[i]]]
  b <- complexity[[gsub(".{4}$", "", colnames(intersect_programs_hc)[i])]]
  corr_complexity[i,2] <- cor(a,b)
}

p2 <- ggplot(corr_complexity, aes(y=corr, x = 1:nrow(corr_complexity))) +
  geom_smooth(span=0.05, se = FALSE, color="gray36", size= 0.8, method = "loess") + 
  geom_point(alpha=0.3, size = 1.5) + 
  theme(axis.ticks = element_blank(), axis.line = element_blank(), panel.background = element_rect(fill = "gray97"), panel.grid = element_blank(), axis.text.x = element_blank(), axis.title = element_text(size = 12), axis.text.y = element_text(size = 10), plot.margin=unit(c(1,1,-0.6,1), "cm") ) + 
  labs(y="Cell complexity\ncorrelation", x="") + 
  scale_y_continuous(expand = c(0,0), limits = c(-0.9,0.9), breaks = seq(-0.8, 0.8, 0.4)) + 
  scale_x_continuous(expand = c(0,0))

## combining plots
egg::ggarrange(p2,p1, nrow=2, ncol = 1, heights = c(1,5))

pdf("nmf_metaprograms_complexity.pdf", height = 6, width = 7, onefile=FALSE)

  # save
  saveRDS(colnames(intersect_programs_hc)[c(424:504, 509:514)], "tech_artifact.RDS")
  
######################## removing metaprograms correlated with complexity ####################
### clustering
intersect_programs_filt_compl <- intersect_programs_hc[-c(424:504, 509:514),-c(424:504, 509:514)]
hc_intersect_compl <- hclust(as.dist(50-intersect_programs_filt_compl), method="average") 
hc_intersect_compl <- reorder(as.dendrogram(hc_intersect_compl), colMeans(intersect_programs_filt_compl))
intersect_programs_hc_compl <- intersect_programs_filt_compl[order.dendrogram(hc_intersect_compl), order.dendrogram(hc_intersect_compl)]


### plotting metaprograms
plot_hc_compl <- reshape2::melt(intersect_programs_hc_compl) 

p3 <- ggplot(data = plot_hc_compl, aes(x=Var1, y=Var2, fill=100*value/(100-value))) + 
  geom_raster() + 
  scale_fill_gradient2(limits=c(2,25), low=gabi_palette[1:111],  mid =gabi_palette[112:222], high = gabi_palette[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_text(size = 13), axis.title = element_text(size = 15), legend.title = element_text(size=12), legend.text = element_text(size = 11), legend.text.align = 0.5, legend.direction = "horizontal", legend.position = c(1.05,-0.1), plot.margin = unit(c(0.3,1,3,1), "cm")) + 
  scale_x_discrete(name="\nPrograms", breaks=colnames(intersect_programs_hc_compl)[seq(100, ncol(intersect_programs_hc_compl), by=100)], labels= seq(100, ncol(intersect_programs_hc_compl), by=100)) + 
  scale_y_discrete(name="", breaks=colnames(intersect_programs_hc_compl)[seq(100, ncol(intersect_programs_hc_compl), by=100)], labels= seq(100, ncol(intersect_programs_hc_compl), by=100)) +
  guides(fill = guide_colourbar(barheight = 0.6, barwidth = 4.8, title.position = "top", title.hjust = 0.5)) +
  panel_border(colour = "gray40", size = 0.4, linetype = 1, remove = FALSE) 

### plotting discrete programs
continous_vs_discr <- readRDS(file.choose())

continous_vs_discr_plot <- data.frame("discr"=is.element(colnames(intersect_programs_hc_compl), names(continous_vs_discr)))

p4 <- ggplot(continous_vs_discr_plot, aes(x=1:nrow(continous_vs_discr_plot), y=0))+
  geom_tile(aes(fill=discr, colour=discr)) + 
  scale_fill_manual(values = c("white", "black"), breaks="TRUE", labels="Discrete programs", name="") +
  scale_colour_manual(values = c("white", "black"), breaks="TRUE", labels="Discrete programs", name="") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), axis.line = element_blank(), axis.text = element_blank(), plot.margin = unit(c(1,5,-0.8,1), "cm"), legend.position = c(1.03,1)) + 
  scale_x_continuous(expand = c(0,0)) +
  labs(x="", y="") +
  panel_border(colour="black", remove = F)

## combining plots
egg::ggarrange(p4,p3, nrow=2, ncol = 1, heights = c(0.4,6))

pdf("nmf_metaprograms_discr_programs_bar.pdf", height = 6.4, width = 7, onefile=FALSE)


################## plotting number of programs per cell line ######################
programs_per_cellline_plot_prefinal <- nmf_all_combined_intersect[prefinal_filter,prefinal_filter]
programs_per_cellline_plot_prefinal <- data.frame(table(sapply(names(w_basis_6), function(x) length(grep(x, colnames(programs_per_cellline_plot_prefinal))))))
programs_per_cellline_plot_prefinal$label <- rep("All", 6)

programs_per_cellline_plot_final <- intersect_programs_filt_compl
programs_per_cellline_plot_final <- data.frame(table(sapply(names(w_basis_6), function(x) length(grep(x, colnames(programs_per_cellline_plot_final))))))
programs_per_cellline_plot_final$label <- rep("Recurrent", 9)

programs_per_cellline_plot <- rbind(programs_per_cellline_plot_prefinal, programs_per_cellline_plot_final)

ggplot(programs_per_cellline_plot, aes(x=as.numeric(as.character(Var1)), y=Freq*100/198)) +
  geom_bar(stat = "identity", fill="gray80") +
  labs(y="Percentage of cell lines", x="Number of continous and discrete programs") +
  facet_grid(cols = vars(label)) +
  theme(strip.background = element_blank()) +
  scale_x_continuous(breaks=0:9)

pdf("nmf_programs_per_cellline.pdf", width = 7.5, height = 3.5)


######################## subgrouping metaprograms ###################
############ non-cell cycle metaprograms
### clustering 
intersect_programs_filt_nc  <- intersect_programs_hc_compl[-c(433:800), -c(433:800)]
hc_intersect_nc <- hclust(as.dist(50-intersect_programs_filt_nc), method="average") 
hc_intersect_nc <- reorder(as.dendrogram(hc_intersect_nc), colMeans(intersect_programs_filt_nc))
intersect_programs_hc_nc <-  intersect_programs_filt_nc[order.dendrogram(hc_intersect_nc), order.dendrogram(hc_intersect_nc)]
  
  # save
  saveRDS(nmf_all_combined[,colnames(intersect_programs_hc_nc)], "nc_programs_hc_ordered.RDS")

### plotting metaprograms
plot_hc_nc <- melt(intersect_programs_hc_nc)

p5 <- ggplot(data = plot_hc_nc, aes(x=Var1, y=Var2, fill=100*value/(100-value))) + 
  geom_raster() +
  scale_fill_gradient2(limits=c(2,25), low=gabi_palette[1:111],  mid =gabi_palette[112:222], high = gabi_palette[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +  
  theme( axis.ticks = element_blank(), panel.background = element_blank(), panel.border=element_rect(colour = "gray40", size = 0.4, fill=F),  axis.line = element_blank(), axis.text.x = element_blank(), axis.text.y=element_text(size=10), axis.title = element_blank(), legend.title = element_text(size=10), legend.text = element_text(size = 9), legend.text.align = 0.5, legend.position = c(-0.4, 0), legend.direction = "horizontal", plot.margin = unit(c(-0.5,-5,-2,0), "cm")) + 
  scale_y_discrete(breaks = colnames(intersect_programs_hc_nc)[seq(50, ncol(intersect_programs_hc_nc), 50)], labels = seq(50, ncol(intersect_programs_hc_nc), 50)) +
  guides(fill = guide_colourbar(barheight = 0.6, barwidth = 4.4, title.position = "top", title.hjust = 0.5)) 
  

### plotting lineage annotation
ccle_annotation <- readRDS("~/CCLE/data processing/scCCLE&pool1_cancer_type_annotation.RDS")

annot_cancer_type_nc <- data.frame("cell_lines" =gsub(".{4}$", "",rownames(intersect_programs_hc_nc)), stringsAsFactors = F)
annot_cancer_type_nc$type <- ccle_annotation$cancer_type[match(annot_cancer_type_nc$cell_lines, rownames(ccle_annotation))]
annot_cancer_type_nc$type[is.element(annot_cancer_type_nc$type, c("Fibroblast", "Gallbladder Cancer", "Neuroblastoma", "Prostate Cancer", "Bone Cancer", "Sarcoma", "Bile Duct Cancer", "Thyroid Cancer"))] <- "Others"
annot_cancer_type_nc$type <- factor(annot_cancer_type_nc$type, levels = unique(sort(annot_cancer_type_nc$type)))

p6 <- ggplot(annot_cancer_type_nc, aes(y=0, x=1:nrow(annot_cancer_type_nc), fill=type, color=type)) +
  geom_tile() +
  scale_color_manual(breaks=c(levels(annot_cancer_type_nc$type)[levels(annot_cancer_type_nc$type)!="Others"], "Others"), values = c(brewer.pal(12, "Set3")[c(1:6,8,7,10:11)], "maroon","gray93", "yellow2", "goldenrod1", "slateblue2"), name="") +
  scale_fill_manual(breaks=c(levels(annot_cancer_type_nc$type)[levels(annot_cancer_type_nc$type)!="Others"], "Others"), values = c(brewer.pal(12, "Set3")[c(1:6,8,7,10:11)], "maroon","gray93", "yellow2", "goldenrod1", "slateblue2"), name="") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 10),  axis.title = element_text(size=8),  plot.margin = unit(c(-1,-5,-0.6,8), "cm"),  legend.position = c(-0.6, -6.5), legend.text.align = 1, legend.key.size = unit(0.4, "cm"), legend.key=element_blank()) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x="", y="") +
  guides(fill = guide_legend(title = "", label.position="left",  override.aes = list(size = 5), ncol = 1), color=FALSE) 

### getting metaprograms
meta1 <- sort(table(nmf_all_combined[,colnames(intersect_programs_hc_nc)[10:25]])/length(10:25), decreasing=T)
write.table(meta1, "metatable_mel_mitf.txt", sep = "\t", quote = F, row.names = F)
meta1 <- names(meta1[meta1 >= 0.25])
write.table(meta1, "meta_mel_mitf.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(colnames(intersect_programs_hc_nc)[10:25], "meta_mel_mitf_programs.txt", sep = "\t", quote = F, row.names = F, col.names = F)
sort(table(gsub(".{4}$", "", gsub("^[^_]*_", "", colnames(intersect_programs_hc_nc)[10:25])))/length(10:25))
                                                 

meta2 <- sort(table(nmf_all_combined[,colnames(intersect_programs_hc_nc)[46:53]])/length(46:53), decreasing=T) 
write.table(meta2, "metatable_mel_emt.txt", sep = "\t", quote = F)
meta2 <- names(meta2[meta2 >= 0.25])
write.table(meta2, "meta_mel_emt.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(colnames(intersect_programs_hc_nc)[46:53], "meta_mel_emt_programs.txt", sep = "\t", quote = F, row.names = F, col.names = F)
sort(table(gsub(".{4}$", "", gsub("^[^_]*_", "", colnames(intersect_programs_hc_nc)[46:53])))/length(46:53))

                                                            
meta3 <- sort(table(nmf_all_combined[,colnames(intersect_programs_hc_nc)[57:75]])/length(57:75), decreasing=T)
write.table(meta3, "metatable_emt.txt", sep = "\t", quote = F, row.names = F)
meta3 <- names(meta3[meta3 >= 0.25])
write.table(meta3, "meta_emt.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(colnames(intersect_programs_hc_nc)[57:75], "metatable_emt_programs.txt", sep = "\t", quote = F, row.names = F, col.names = F)
sort(table(gsub(".{4}$", "", gsub("^[^_]*_", "", colnames(intersect_programs_hc_nc)[57:75])))/length(57:75))                                                          

meta4 <- sort(table(nmf_all_combined[,colnames(intersect_programs_hc_nc)[81:95]])/length(81:95), decreasing=T)
write.table(meta4, "metatable_ifn.txt", sep = "\t", quote = F, row.names = F)
meta4 <- names(meta4[meta4 >= 0.25])
write.table(meta4, "meta_ifn.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(colnames(intersect_programs_hc_nc)[81:95], "meta_ifn_programs.txt", sep = "\t", quote = F, row.names = F, col.names = F)
sort(table(gsub(".{4}$", "", gsub("^[^_]*_", "", colnames(intersect_programs_hc_nc)[81:95])))/length(81:95))                                                             

meta5 <- sort(table(nmf_all_combined[,colnames(intersect_programs_hc_nc)[98:109]])/length(98:109), decreasing=T)
write.table(meta5, "metatable_emt_snc.txt", sep = "\t", quote = F, row.names = F)
meta5 <- names(meta5[meta5 >= 0.25])
write.table(meta5, "meta_emt_snc.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(colnames(intersect_programs_hc_nc)[98:109], "meta_emt_snc_programs.txt", sep = "\t", quote = F, row.names = F, col.names = F)
sort(table(gsub(".{4}$", "", gsub("^[^_]*_", "", colnames(intersect_programs_hc_nc)[98:109])))/length(98:109))
                                                       

meta6 <- sort(table(nmf_all_combined[,colnames(intersect_programs_hc_nc)[116:132]])/length(116:132), decreasing = T)
write.table(meta6, "metatable_snc_classic.txt", sep = "\t", quote = F, row.names = F)
meta6 <- names(meta6[meta6 >= 0.25])
write.table(meta6, "meta_snc_classic.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(colnames(intersect_programs_hc_nc)[116:132], "meta_snc_classic_programs.txt", sep = "\t", quote = F, row.names = F, col.names = F)
sort(table(gsub(".{4}$", "", gsub("^[^_]*_", "", colnames(intersect_programs_hc_nc)[116:132])))/length(116:132))
                                                           

meta7 <- sort(table(nmf_all_combined[,colnames(intersect_programs_hc_nc)[140:167]])/length(140:167), decreasing=T)
write.table(meta7, "metatable_snc_epi.txt", sep = "\t", quote = F, row.names = F)
meta7 <- names(meta7[meta7 >= 0.25])
write.table(meta7, "meta_snc_epi.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(colnames(intersect_programs_hc_nc)[140:167], "meta_snc_epi_programs.txt", sep = "\t", quote = F, row.names = F, col.names = F)
sort(table(gsub(".{4}$", "", gsub("^[^_]*_", "", colnames(intersect_programs_hc_nc)[140:167])))/length(140:167))
                                               

meta8 <- sort(table(nmf_all_combined[,colnames(intersect_programs_hc_nc)[197:250]])/length(197:250), decreasing=T)
write.table(meta8, "metatable_stress.txt", sep = "\t", quote = F, row.names = F)
meta8 <- names(meta8[meta8 >= 0.25])
write.table(meta8, "meta_stress.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(colnames(intersect_programs_hc_nc)[197:250], "meta_stress_programs.txt", sep = "\t", quote = F, row.names = F, col.names = F)
sort(table(gsub(".{4}$", "", gsub("^[^_]*_", "", colnames(intersect_programs_hc_nc)[197:250])))/length(197:250))

meta9 <- sort(table(nmf_all_combined[,colnames(intersect_programs_hc_nc)[261:278]])/length(261:278), decreasing=T)
write.table(meta9, "metatable_endo.txt", sep = "\t", quote = F, row.names = F)
meta9 <- names(meta9[meta9 >= 0.25])
write.table(meta9, "meta_endo.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(colnames(intersect_programs_hc_nc)[261:278], "meta_endo_programs.txt", sep = "\t", quote = F, row.names = F, col.names = F)
sort(table(gsub(".{4}$", "", gsub("^[^_]*_", "", colnames(intersect_programs_hc_nc)[261:278])))/length(261:278))

meta10 <- sort(table(nmf_all_combined[,colnames(intersect_programs_hc_nc)[278:432]])/length(278:432), decreasing=T)
write.table(meta10, "metatable_prot.txt", sep = "\t", quote = F, row.names = F)
meta10 <- names(meta10[meta10 >= 0.25])
write.table(meta10, "meta_prot.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(colnames(intersect_programs_hc_nc)[278:432], "meta_prot_programs.txt", sep = "\t", quote = F, row.names = F, col.names = F)
sort(table(gsub(".{4}$", "", gsub("^[^_]*_", "", colnames(intersect_programs_hc_nc)[278:432])))/length(278:432))


### plotting NMF gene scores
nc_genes <- c(rev(meta1), rev(meta2), rev(meta3), rev(meta4), rev(meta5), rev(meta6), rev(meta7), rev(meta8), rev(meta9), rev(meta10))
nc_genes_matrix <- sapply(gene_score_all[colnames(intersect_programs_hc_nc)], function(x) x[nc_genes])
rownames(nc_genes_matrix) <- NULL
nc_genes_matrix[is.na(nc_genes_matrix)] <- 0

nc_genes_matrix_plot <- melt(nc_genes_matrix)

p8 <- ggplot(data = nc_genes_matrix_plot, aes(x=Var2, y=Var1, fill=value*100)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(0,8), low= "#b5cce3", mid=  "white", high = "darkred", midpoint = 4,   oob=squish, name=expression(paste("NMF score (10"^"2", ")", sep = ""))) +
   scale_color_gradient2(limits=c(0,8), low= "#b5cce3", mid=  "white", high = "darkred", midpoint = 4,   oob=squish, name=expression(paste("NMF score (10"^"2", ")", sep = ""))) + 
  theme( axis.ticks = element_blank(), panel.background = element_blank(), panel.border=element_rect(colour = "black", size = 0.4, fill=F),   axis.line = element_blank(), axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.title = element_text(size=10), legend.text = element_text(size = 9), legend.text.align = 0.5, legend.position = c(-0.4, 0.2), legend.direction = "horizontal",  plot.margin = unit(c(0,0,1,0), "cm")) +
  labs(x="\nPrograms", y="Genes") +
  scale_x_discrete(breaks = colnames(intersect_programs_hc_nc)[seq(50, ncol(intersect_programs_hc_nc), 50)], labels = seq(50, ncol(intersect_programs_hc_nc), 50)) +
  scale_y_continuous(expand = c(0,0), breaks= seq(50, ncol(intersect_programs_hc_nc), 50)) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 4.4, barheight = 0.6)) +
  geom_vline(xintercept=c(25, 53, 75, 95, 109, 132, 167, 250, 278), color = "black", size=0.2) 
                                                     

### plotting gene annotation
go_bp <- data.frame(t(read.table("~/CCLE/datasets/genesets/GO.c5.bp.v7.1.symbols.txt", sep = "\t", header = F, comment.char = "", quote = "",row.names = 1)), stringsAsFactors = F)
snc_sig <- scan("~/CCLE/datasets/genesets/kerat_snc_hernandez-segura2017.txt", sep="\t", what="character")
hallmarks <- data.frame(t(read.table("~/CCLE/datasets/genesets/hallmarks_MSigDB.v6.2.txt", sep = "\t", header = F, comment.char = "", quote = "",row.names = 1)), stringsAsFactors = F)

extra <- scan("~/CCLE/datasets/genesets/GO_EXTRACELLULAR_SPACE.txt", sep="\t", what="character")                          
                          
gene_uni <- names(sort(rowMeans(sapply(expr, rowMeans)), decreasing=T))[1:7000]
nc_meta_list <- list(meta1=meta1, meta2=meta2, meta3=meta3, meta4=meta4, meta5=meta5, meta6=meta6, meta7=meta7, meta8=meta8, meta9=meta9, meta10=meta10)  
nc_meta_list <- lapply(nc_meta_list, function(x) x[is.element(x,gene_uni)])        
geneset_list <- list( GO_PIGMENT_CELL_DIFFERENTIATION= go_bp$GO_PIGMENT_CELL_DIFFERENTIATION, GO_WOUND_HEALING=   go_bp$GO_WOUND_HEALING, GO_RESPONSE_TO_TYPE_I_INTERFERON=  go_bp$GO_RESPONSE_TO_TYPE_I_INTERFERON, snc_sig =  snc_sig, extra =extra, GO_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS = go_bp$GO_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS, GO_PROTEIN_FOLDING= go_bp$GO_PROTEIN_FOLDING,  GO_PROTEASOMAL_PROTEIN_CATABOLIC_PROCESS= go_bp$GO_PROTEASOMAL_PROTEIN_CATABOLIC_PROCESS)
geneset_list <- lapply(geneset_list, function(x) x[is.element(x,gene_uni)])       
                       
hyper_test <-  sapply(nc_meta_list, function(x) sapply(geneset_list, function(y) phyper(q=length(intersect(x,y))-1, m=length(x), n=length(gene_uni)-length(x), k=length(y), lower.tail=F)))                                       
hyper_test <- t(apply(hyper_test, 2, function(x) p.adjust(x, "fdr", n=length(hyper_test))))                                           

nc_genes_annotation <- data.frame(programs=rep(c("meta1","meta2","meta3", "meta4", "meta5", "meta6", "meta7", "meta8", "meta9", "meta10"), sapply(list(meta1, meta2, meta3, meta4, meta5, meta6, meta7, meta8, meta9, meta10), length)))   
nc_genes_annotation <- hyper_test[match(nc_genes_annotation$programs, rownames(hyper_test)),]  
rownames(nc_genes_annotation) <- NULL

grey_palette <- brewer.pal( 9, "Greys")  
                      
nc_genes_annotation <- melt(as.matrix(nc_genes_annotation))    
                      
p9 <- ggplot(nc_genes_annotation, aes(y=Var1, x=""))+
  geom_tile(aes(fill=-log10(value), color= -log10(value))) + 
  scale_fill_gradient2(limits=c(1, 8), midpoint = 4.5, low= c( "white",grey_palette[1:3]), mid= grey_palette[4:6], high = grey_palette[7:9] ,   oob=squish, name="P adjust\n(fdr)")+
  scale_color_gradient2(limits=c(1, 8), midpoint = 4.5, low= c("white",grey_palette[1:3]), mid= grey_palette[4:6], high = grey_palette[7:9] ,   oob=squish, name="P adjust\n(fdr)") +
  theme(axis.ticks = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_blank(), panel.border = element_rect(color="black", fill=F)) + 
  labs(x="", y="") +
  facet_grid(cols = vars(Var2)) +
  theme(panel.spacing = unit(0.8, "lines"), strip.text=element_blank()) +
   scale_x_discrete(expand=c(0,0)) +                                      
   scale_y_continuous(expand=c(0,0))                  
                                                                       
nc_genes_annotation_labels <- data.frame(genes=c(rev(meta1), rev(meta2), rev(meta3), rev(meta4), rev(meta5), rev(meta6), rev(meta7), rev(meta8), rev(meta9), rev(meta10)))                     
nc_genes_annotation_labels$labels <- NA                   
nc_genes_annotation_labels[match(c("MITF", "PMEL", "CYR61", "CTGF","VIM", "AXL", "ISG15", "IFI6", "LAMB3", "JUP", "TP53TG1", "TP53I3", "AQP3", "CLDN4", "ATF3", "DDIT3", "HSPA5", "OS9", "PSMA3", "PSMA4"), nc_genes_annotation_labels$genes), "labels"] <- c("MITF", "PMEL", "CYR61", "CTGF","VIM", "AXL", "ISG15", "IFI6" ,"LAMB3", "JUP", "TP53TG1", "TP53I3", "AQP3", "CLDN4", "ATF3", "DDIT3", "HSPA5", "OS9", "PSMA3", "PSMA4")                      

pdf("metaprograms_nc_complete_legend.pdf", height = 5, width = 2)      
ggplot(nc_genes_annotation_labels, aes(y=1:nrow(nc_genes_annotation_labels), x=0, label=labels)) +
   geom_text_repel(force=2,nudge_x=-0.03, direction = "y") +
   theme(panel.background = element_blank(), axis.text=element_blank(), axis.title=element_blank(), axis.ticks = element_blank()) +
   geom_vline(xintercept=0) +
   scale_y_continuous(expand=c(0,0))                       
dev.off()    
                          
                          
### combining plots
pdf("metaprograms_nc_complete.pdf", height = 8.2, width = 11.4, onefile=FALSE)
egg::ggarrange( p6, ggplot()+theme(panel.background = element_blank()), ggplot(), p5, ggplot()+theme(panel.background = element_blank()), ggplot(), p8, p9, ggplot(),nrow = 3, ncol = 3, widths = c(5,4,0), heights = c(3,50, 40))
dev.off()



################# cell cycle metaprograms
###clustering
intersect_programs_filt_cc  <- intersect_programs_hc_compl[c(433:800),c(433:800)]
hc_intersect_cc <- hclust(as.dist(50-intersect_programs_filt_cc), method="average") 
hc_intersect_cc <- reorder(as.dendrogram(hc_intersect_cc), colMeans(intersect_programs_filt_cc))
intersect_programs_hc_cc <-  intersect_programs_filt_cc[order.dendrogram(hc_intersect_cc), order.dendrogram(hc_intersect_cc)]

  # save
  saveRDS(nmf_all_combined[,colnames(intersect_programs_hc_cc)], "cc_programs_hc_ordered.RDS")

### plotting metaprograms
plot_hc_cc <- melt(intersect_programs_hc_cc) 

p11 <- ggplot(data = plot_hc_cc, aes(x=Var1, y=Var2, fill=value)) + 
  geom_raster() +
  scale_fill_gradient2(limits=c(2,20), low=gabi_palette[1:111],  mid =gabi_palette[112:222], high = gabi_palette[223:333], midpoint = 11, oob=squish, name="Similarity") +  
  theme( axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.title = element_text(size=10), legend.text = element_text(size = 9), legend.text.align = 0.5, legend.position = c(-0.3, 0.05), legend.direction = "horizontal", plot.margin = unit(c(-0.5,-0.6,-2,0), "cm")) + 
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5,  barwidth = 4, barheight = 0.7)) +
  panel_border(colour = "gray40", size = 0.4, linetype = 1, remove = FALSE) 


### plotting lineage annotation
annot_cc <- data.frame("cell_lines" = rownames(intersect_programs_hc_cc), "lineage" =  gsub("_", " ", gsub(".{4}$", "", gsub("^[^_]*_", "", rownames(intersect_programs_hc_cc)))), stringsAsFactors = F)
annot_cc$lineage[is.element(annot_cc$lineage, c("AUTONOMIC GANGLIA", "BONE", "PROSTATE", "PLEURA", "SOFT TISSUE", "THYROID", "BILIARY TRACT", "KIDNEY", "STOMACH", "URINARY TRACT", "LIVER", "OESOPHAGUS"))] <- "OTHERS"
annot_cc$lineage <- factor(annot_cc$lineage, levels= c("BREAST", "CENTRAL NERVOUS SYSTEM", "ENDOMETRIUM", "LARGE INTESTINE",  "LUNG","OVARY" ,  "PANCREAS", "SKIN","UPPER AERODIGESTIVE TRACT", "OTHERS" ))

colors_brewer <- colorRampPalette(brewer.pal(12, "Set3"))

p12 <- ggplot(annot_cc, aes(x=factor(cell_lines, levels = cell_lines), y = 1)) + 
  geom_tile(aes(fill=lineage, colour=lineage)) + 
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 8),  legend.box.margin = unit(c(11,0,5,1), "cm"), plot.margin = unit(c(-1,0,-1,0), "cm"),  legend.position = "left", legend.text.align = 1, legend.justification = "right", legend.key.size = unit(0.4, "cm")) + 
  labs(x="", y="") + 
  guides(fill = guide_legend(title = "", label.position="left",  override.aes = list(linetype = 0, size = 5), ncol = 1), color=FALSE)  + 
  scale_color_manual(values =c(colors_set3[c(7, 2:4, 6,5,8,1, 10)], "gray95")) +
  scale_fill_manual(values =c(colors_set3[c(7, 2:4, 6,5,8,1, 10)], "gray95")) 

### getting metaprograms
meta11 <- table(nmf_all_combined[,colnames(intersect_programs_hc_cc)[25:159]])/length(25:159)
write.table(sort(meta11, decreasing = T), "metatable_g1.txt", sep = "\t", quote = F, row.names = F)
meta11 <- names(meta11[meta11 >= 0.25])
write.table(meta11, "meta_g1.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(colnames(intersect_programs_hc_cc)[25:159], "meta_g1_programs.txt", sep = "\t", quote = F, row.names = F, col.names = F)
sort(table(gsub(".{4}$", "", gsub("^[^_]*_", "", colnames(intersect_programs_hc_cc)[25:159])))/length(25:159))

meta12 <- table(nmf_all_combined[,colnames(intersect_programs_hc_cc)[165:368]])/length(165:368)
write.table(sort(meta12, decreasing = T), "metatable_g2.txt", sep = "\t", quote = F, row.names = F)
meta12 <- names(meta12[meta12 >= 0.25])
write.table(meta12, "meta_g2.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(colnames(intersect_programs_hc_cc)[165:368], "meta_g2_programs.txt", sep = "\t", quote = F, row.names = F, col.names = F)
sort(table(gsub(".{4}$", "", gsub("^[^_]*_", "", colnames(intersect_programs_hc_cc)[165:368])))/length(165:368))

### plotting correlation with  cnv clones
cnv_corr_cc <- data.frame("program" = colnames(intersect_programs_hc_cc), "cnv"=NA)
cnv_corr_cc$cnv[is.element(cnv_corr_cc$program, names(conti_gmm))] <- "cnv"
cnv_corr_cc$cnv[!is.element(cnv_corr_cc$program, names(conti_gmm))] <- "no_cnv"

p13 <- ggplot(cnv_corr_cc, aes(x=as.factor(1:nrow(cnv_corr_cc)), y = 1)) + 
  geom_tile(aes(fill=cnv, colour=cnv)) + 
  scale_fill_manual(values = c("steelblue4", "white"), breaks="cnv", labels="  CNV clone-derived") +
  scale_colour_manual(values = c("steelblue4", "white")) +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 10), legend.title = element_blank(), plot.margin = unit(c(1,-5,-5,0), "cm")) + 
  labs(x="", y="") +
  guides(colour=F) +
  panel_border(colour="gray60", remove = F)


### plotting NMF gene scores
cc_genes <- c(meta11, meta12)
cc_genes_matrix <- sapply(gene_score_all[colnames(intersect_programs_hc_cc)], function(x) x[cc_genes])
rownames(cc_genes_matrix) <- cc_genes
cc_genes_matrix[is.na(cc_genes_matrix)] <- 0

plot_cc_genes <- melt(cc_genes_matrix)

p14 <- ggplot(data = plot_cc_genes, aes(x=Var2, y=Var1, fill=value*100)) + 
  geom_raster(interpolate = T) + 
  scale_fill_gradient2(limits=c(0,6), low=c("slategray3"), mid= c("white"), high = c("darkred"), midpoint = 3,   oob=squish, name=expression(paste("NMF score (10"^"2", ")", sep = ""))) +
  theme( axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10), legend.title = element_text(size=10), legend.text = element_text(size = 9), legend.text.align = 0.5, legend.position = c(-0.3, 0.2), legend.direction = "horizontal",  plot.margin = unit(c(0,-0.6,1,0), "cm")) +
  labs(x="\nPrograms", y="") +
  scale_x_discrete(breaks = colnames(intersect_programs_hc_cc)[seq(20, ncol(intersect_programs_hc_cc), 50)], labels = seq(20, ncol(intersect_programs_hc_cc), 50)) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 4, barheight = 0.7)) +
  geom_vline(xintercept=170, color = "black", size=0.2) +
  panel_border(colour = "black", size = 0.4, linetype = 1, remove = FALSE) 


### plotting gene annotation
cc_genes_annotation <- data.frame( "G1/S transition"= is.element(cc_genes, go_bp$GO_CELL_CYCLE_G1_S_PHASE_TRANSITION), "DNA replication"= is.element(cc_genes, go_bp$GO_DNA_REPLICATION), "G2/M transition"=  is.element(cc_genes, go_bp$GO_CELL_CYCLE_G2_M_PHASE_TRANSITION), "Chromosome segregation"=  is.element(cc_genes, go_bp$GO_CHROMOSOME_SEGREGATION))

cc_genes_annotation <- melt(as.matrix(cc_genes_annotation))

p15 <- ggplot(cc_genes_annotation, aes(y=Var1, x=0))+
  geom_tile(aes(fill=value, colour=value)) + 
  scale_fill_manual(values = c("white", "black")) +
  scale_colour_manual(values = c("white", "black")) +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), axis.line = element_blank(), axis.text = element_blank()) + 
  guides(fill=F, color=F) +
  labs(x="", y="") +
  panel_border(colour="black", remove = F) +
  facet_grid(cols = vars(Var2)) +
  theme(panel.spacing = unit(0.5, "lines"), strip.text=element_blank(), plot.margin = unit(c(0,0,1,-0.3), "cm")) 

### combining plots
egg::ggarrange( p13,ggplot(), ggplot(),  ggplot(), p12, ggplot(), ggplot(), ggplot(), p11, ggplot(), ggplot(),  ggplot(), p14, p15, ggplot(),  ggplot(), nrow = 4, ncol = 4, widths = c(5,1.5,0,0), heights = c(3,3,50,15))
pdf("cellcycle&lineage.pdf", height = 7.3, width = 9.5, onefile=FALSE)














### plotting correlation with  cnv clones
conti_gmm <- readRDS(file.choose())
cnv_corr_nc <- data.frame("program" = colnames(intersect_programs_hc_nc), "cnv"=NA)
cnv_corr_nc$cnv[is.element(cnv_corr_nc$program, names(conti_gmm))] <- "cnv"
cnv_corr_nc$cnv[!is.element(cnv_corr_nc$program, names(conti_gmm))] <- "no_cnv"

p8 <- ggplot(cnv_corr_nc, aes(x=as.factor(1:nrow(cnv_corr_nc)), y = 1)) + 
  geom_tile(aes(fill=cnv, colour=cnv)) + 
  scale_fill_manual(values = c("steelblue4", "white"), breaks="cnv", labels="  CNV clone-derived") +
  scale_colour_manual(values = c("steelblue4", "white")) +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 10), legend.title = element_blank(), plot.margin = unit(c(1,-5,-5,0), "cm")) + 
  labs(x="", y="") +
  guides(colour=F) +
  panel_border(colour="gray60", remove = F)
