############################################################################################################################################################################################### DISCRETE PROGRAMS OF VARIABILITY  ########################################################

### reading TPM per cell lines
expr_tpm <- readRDS("/home/labs/tirosh/kinker/analyses/CCLE/data processing/scCCLE&pool1_tpm_by_cell_line.rds") 

##################### minpt 5 eps 1.8 
### reading discrete programs (mipt 5 eps 1.8)
discr_programs <-  readRDS("/home/labs/tirosh/kinker/analyses/CCLE/discrete programs/discr_clusters_minpt5_eps1.8.RDS") 

### creating cluster signatures
discr_programs_sig <- list()
for(i in names(discr_programs)) {
  a <- discr_programs[[i]][["clusters_sig"]]
  a <- lapply(a, function(x) x[x[,"log2.FC."] >= 1 & x[,"ttest_p"] < 0.001,])
  a <- lapply(a, function(x) rownames(x)[order(x[,"log2.FC."], decreasing = T)][1:50])
  a <- lapply(a, function(x) x[!is.na(x)])
  discr_programs_sig[[i]] <- a
}
discr_programs_sig <- unlist(discr_programs_sig, recursive = F)

### filtering out huge clusters
discr_programs_cells <- unlist(lapply(discr_programs, function(x) x[["clusters_cells"]]), recursive = F)

max_size_cutoff <- c()
for(i in names(discr_programs_cells)) {
  max_size_cutoff <- c(max_size_cutoff, length(discr_programs_cells[[i]])/ncol(expr_tpm[[sub("\\..*", "", i)]]) < 0.90)
}

discr_programs_sig <- discr_programs_sig[max_size_cutoff]
discr_programs_sig <- discr_programs_sig[sapply(discr_programs_sig, length) > 0]

### removing 2 redundant clusters
discr_programs_sig$RERFLCAI_LUNG.1 <- NULL
discr_programs_sig$RERFLCAI_LUNG.2 <- NULL

### calculating the intersection between programs (jaccard index)
programs_intersect <- sapply(discr_programs_sig, function(x) sapply(discr_programs_sig,function(y) length(intersect(y,x))/length(union(x,y))))*100

hc_intersect <- hclust(as.dist(100-programs_intersect), method="average") 
hc_intersect <- reorder(as.dendrogram(hc_intersect), colMeans(programs_intersect))
intersect_programs_hc <- programs_intersect[order.dendrogram(hc_intersect), order.dendrogram(hc_intersect)]

### plotting metaprograms
gabi_palette <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

# heatmap
plot_hc <- reshape2::melt(intersect_programs_hc) 
p1 <- ggplot(data = plot_hc, aes(x=Var1, y=Var2, fill=value, color=value)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(2,25), low=gabi_palette[1:111],  mid =gabi_palette[112:222], high = gabi_palette[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  scale_color_gradient2(limits=c(2, 25), low=gabi_palette[1:111],  mid =gabi_palette[112:222], high = gabi_palette[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  theme( axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(),  axis.text = element_text(size = 13), axis.title = element_text(size = 15), legend.title = element_text(size=12), legend.text = element_text(size = 11), legend.text.align = 0.5, legend.direction = "horizontal", legend.position = c(1.05,-0.1), plot.margin = unit(c(2,6,2,2), "cm")) + 
  scale_x_discrete(name="Discrete programs", breaks=colnames(intersect_programs_hc)[seq(5, ncol(intersect_programs_hc), by=5)], labels= seq(5, ncol(intersect_programs_hc), by=5)) + 
  scale_y_discrete(name="", breaks=colnames(intersect_programs_hc)[seq(5, ncol(intersect_programs_hc), by=5)], labels= seq(5, ncol(intersect_programs_hc), by=5)) +
  panel_border(colour = "black", size = 0.4, linetype = 1, remove = FALSE)  +
  guides(fill = guide_colourbar(barheight = 0.8, barwidth = 4.8, title.position = "top", title.hjust = 0.5))


# cancer type annotation
ccle_annotation <- readRDS(file.choose())
ccle_annotation$cancer_type[is.element(ccle_annotation$cancer_type, c("Fibroblast", "Gallbladder Cancer", "Neuroblastoma", "Prostate Cancer", "Bone Cancer", "Sarcoma", "Bile Duct Cancer", "Thyroid Cancer"))] <- "Others"

color_match <- data.frame("type"=sort(unique(ccle_annotation$cancer_type)), "color"= c(brewer.pal(12, "Set3")[c(1:6,8,7,10:11)], "maroon","gray93", "yellow2", "goldenrod1", "slateblue2"),stringsAsFactors = F)

annot_cancer_type <- data.frame("cell_lines" =sub("\\..*", "", rownames(intersect_programs_hc)), stringsAsFactors = F)
annot_cancer_type$type <- ccle_annotation$cancer_type[match(annot_cancer_type$cell_lines, rownames(ccle_annotation))]

p2 <- ggplot(annot_cancer_type, aes(y=0, x=1:nrow(annot_cancer_type), fill=type, color=type)) +
  geom_tile() +
  scale_color_manual(breaks=c(sort(unique(annot_cancer_type$type))[sort(unique(annot_cancer_type$type))!="Others"], "Others"), values = color_match$color[match(sort(unique(annot_cancer_type$type)), color_match$type)] , name="") +
  scale_fill_manual(breaks=c(sort(unique(annot_cancer_type$type))[sort(unique(annot_cancer_type$type))!="Others"], "Others"), values = color_match$color[match(sort(unique(annot_cancer_type$type)), color_match$type)] , name="") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 11),  axis.title = element_text(size=8),  plot.margin = unit(c(2,4,-2.5,1), "cm"), legend.position = c(1,-2), legend.text.align = 0, legend.key.size = unit(0.4, "cm")) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x="", y="") +
  guides(fill = guide_legend(title = "",  override.aes = list(size = 5), ncol=1), color=FALSE) 

egg::ggarrange(p2,p1, nrow=2, heights = c(0.5,9))

pdf("discrete_programs_eps1.8.pdf", height = 6.5, width = 8, onefile = F)


cluster_1 <- table(unlist(discr_programs_sig[colnames(intersect_programs_hc)[4:9]]))
cluster_1 <- sort(cluster_1, decreasing = T)[1:30]

cluster_2 <- table(unlist(discr_programs_sig[colnames(intersect_programs_hc)[12:18]]))
cluster_2 <- sort(cluster_2, decreasing = T)[1:30]

##################### minpt 5 eps 1.5
### reading discrete programs (mipt 5 eps 1.5)
discr_programs <- discr_programs <-  readRDS("/home/labs/tirosh/kinker/analyses/CCLE/discrete programs/discr_clusters_minpt5_eps1.5.RDS") 


### creating cluster signatures
discr_programs_sig <- list()
for(i in names(discr_programs)) {
  a <- discr_programs[[i]][["clusters_sig"]]
  a <- lapply(a, function(x) x[x[,"log2.FC."] >= 1 & x[,"ttest_p"] < 0.001,])
  a <- lapply(a, function(x) rownames(x)[order(x[,"log2.FC."], decreasing = T)][1:50])
  a <- lapply(a, function(x) x[!is.na(x)])
  discr_programs_sig[[i]] <- a
}
discr_programs_sig <- unlist(discr_programs_sig, recursive = F)

### filtering out huge clusters
discr_programs_cells <- unlist(lapply(discr_programs, function(x) x[["clusters_cells"]]), recursive = F)

max_size_cutoff <- c()
for(i in names(discr_programs_cells)) {
  max_size_cutoff <- c(max_size_cutoff, length(discr_programs_cells[[i]])/ncol(expr_tpm[[sub("\\..*", "", i)]]) < 0.90)
}

discr_programs_sig <- discr_programs_sig[max_size_cutoff]
discr_programs_sig <- discr_programs_sig[sapply(discr_programs_sig, length) > 0]

### removing 2 redundant clusters
discr_programs_sig$RERFLCAI_LUNG.1 <- NULL
discr_programs_sig$RERFLCAI_LUNG.2 <- NULL

### calculating the intersection between programs (jaccard index)
programs_intersect <- sapply(discr_programs_sig, function(x) sapply(discr_programs_sig,function(y) length(intersect(y,x))/length(union(x,y))))*100

hc_intersect <- hclust(as.dist(100-programs_intersect), method="average") 
hc_intersect <- reorder(as.dendrogram(hc_intersect), colMeans(programs_intersect))
intersect_programs_hc <- programs_intersect[order.dendrogram(hc_intersect), order.dendrogram(hc_intersect)]

### plotting metaprograms
gabi_palette <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

# heatmap
plot_hc <- reshape2::melt(intersect_programs_hc) 
p3 <- ggplot(data = plot_hc, aes(x=Var1, y=Var2, fill=value, color=value)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(2,25), low=gabi_palette[1:111],  mid =gabi_palette[112:222], high = gabi_palette[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  scale_color_gradient2(limits=c(2, 25), low=gabi_palette[1:111],  mid =gabi_palette[112:222], high = gabi_palette[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  theme( axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(),  axis.text = element_text(size = 13), axis.title = element_text(size = 15), legend.title = element_text(size=12), legend.text = element_text(size = 11), legend.text.align = 0.5, legend.direction = "horizontal", legend.position = c(1.05,-0.1), plot.margin = unit(c(2,4,2,2), "cm")) + 
  scale_x_discrete(name="Discrete programs", breaks=colnames(intersect_programs_hc)[seq(10, ncol(intersect_programs_hc), by=10)], labels= seq(10, ncol(intersect_programs_hc), by=10)) + 
  scale_y_discrete(name="", breaks=colnames(intersect_programs_hc)[seq(10, ncol(intersect_programs_hc), by=10)], labels= seq(10, ncol(intersect_programs_hc), by=10)) +
  panel_border(colour = "black", size = 0.4, linetype = 1, remove = FALSE)  +
  guides(fill = guide_colourbar(barheight = 0.8, barwidth = 4.8, title.position = "top", title.hjust = 0.5))

# cancer type annotation
ccle_annotation <- readRDS(file.choose())
ccle_annotation$cancer_type[is.element(ccle_annotation$cancer_type, c("Fibroblast", "Gallbladder Cancer", "Neuroblastoma", "Prostate Cancer", "Bone Cancer", "Sarcoma", "Bile Duct Cancer", "Thyroid Cancer"))] <- "Others"

color_match <- data.frame("type"=sort(unique(ccle_annotation$cancer_type)), "color"= c(brewer.pal(12, "Set3")[c(1:6,8,7,10:11)], "maroon","gray93", "yellow2", "goldenrod1", "slateblue2"),stringsAsFactors = F)

annot_cancer_type <- data.frame("cell_lines" =sub("\\..*", "", rownames(intersect_programs_hc)), stringsAsFactors = F)
annot_cancer_type$type <- ccle_annotation$cancer_type[match(annot_cancer_type$cell_lines, rownames(ccle_annotation))]

p4 <- ggplot(annot_cancer_type, aes(y=0, x=1:nrow(annot_cancer_type), fill=type, color=type)) +
  geom_tile() +
  scale_color_manual(breaks=c(sort(unique(annot_cancer_type$type))[sort(unique(annot_cancer_type$type))!="Others"], "Others"), values = color_match$color[match(sort(unique(annot_cancer_type$type)), color_match$type)] , name="") +
  scale_fill_manual(breaks=c(sort(unique(annot_cancer_type$type))[sort(unique(annot_cancer_type$type))!="Others"], "Others"), values = color_match$color[match(sort(unique(annot_cancer_type$type)), color_match$type)] , name="") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 11),  axis.title = element_text(size=8),  plot.margin = unit(c(3.5,4,-2.5,1), "cm"), legend.position = c(0,4), legend.text.align = 0, legend.key.size = unit(0.4, "cm")) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x="", y="") +
  guides(fill = guide_legend(title = "",  override.aes = list(size = 5), ncol=3), color=FALSE) 

egg::ggarrange(p4,p3, nrow=2, heights = c(0.5,9))

pdf("discrete_programs_eps1.5.pdf", height = 7, width = 7, onefile = F)

cluster_1 <- table(unlist(discr_programs_sig[colnames(intersect_programs_hc)[36:44]]))
cluster_1 <- sort(cluster_1, decreasing = T)[1:30]

cluster_2 <- table(unlist(discr_programs_sig[colnames(intersect_programs_hc)[49:55]]))
cluster_2 <- sort(cluster_2, decreasing = T)[1:30]

cluster_3 <- table(unlist(discr_programs_sig[colnames(intersect_programs_hc)[55:63]]))
cluster_3 <- sort(cluster_3, decreasing = T)[1:30]


##################### minpt 5 eps 1.2 
### reading discrete programs (mipt 5 eps 1.2)
discr_programs <- readRDS("/home/labs/tirosh/kinker/analyses/CCLE/discrete programs/discr_clusters_minpt5_eps1.2.RDS") 

### creating cluster signatures
discr_programs_sig <- list()
for(i in names(discr_programs)) {
  a <- discr_programs[[i]][["clusters_sig"]]
  a <- lapply(a, function(x) x[x[,"log2.FC."] >= 1 & x[,"ttest_p"] < 0.001,])
  a <- lapply(a, function(x) rownames(x)[order(x[,"log2.FC."], decreasing = T)][1:50])
  a <- lapply(a, function(x) x[!is.na(x)])
  discr_programs_sig[[i]] <- a
}
discr_programs_sig <- unlist(discr_programs_sig, recursive = F)

### filtering out huge clusters
discr_programs_cells <- unlist(lapply(discr_programs, function(x) x[["clusters_cells"]]), recursive = F)

max_size_cutoff <- c()
for(i in names(discr_programs_cells)) {
  max_size_cutoff <- c(max_size_cutoff, length(discr_programs_cells[[i]])/ncol(expr_tpm[[sub("\\..*", "", i)]]) < 0.90)
}

discr_programs_sig <- discr_programs_sig[max_size_cutoff]
discr_programs_sig <- discr_programs_sig[sapply(discr_programs_sig, length) > 0]


### calculating the intersection between programs (jaccard index)
programs_intersect <- sapply(discr_programs_sig, function(x) sapply(discr_programs_sig,function(y) length(intersect(y,x))/length(union(x,y))))*100

hc_intersect <- hclust(as.dist(100-programs_intersect), method="average") 
hc_intersect <- reorder(as.dendrogram(hc_intersect), colMeans(programs_intersect))
intersect_programs_hc <- programs_intersect[order.dendrogram(hc_intersect), order.dendrogram(hc_intersect)]

### plotting metaprograms
gabi_palette <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

# heatmap
plot_hc <- reshape2::melt(intersect_programs_hc) 
p5 <- ggplot(data = plot_hc, aes(x=Var1, y=Var2, fill=value, color=value)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(2,25), low=gabi_palette[1:111],  mid =gabi_palette[112:222], high = gabi_palette[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  scale_color_gradient2(limits=c(2, 25), low=gabi_palette[1:111],  mid =gabi_palette[112:222], high = gabi_palette[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  theme( axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(),  axis.text = element_text(size = 13), axis.title = element_text(size = 15), legend.title = element_text(size=12), legend.text = element_text(size = 11), legend.text.align = 0.5, legend.direction = "horizontal", legend.position = c(1.05,-0.1), plot.margin = unit(c(2,4,2,2), "cm")) + 
  scale_x_discrete(name="Discrete programs", breaks=colnames(intersect_programs_hc)[seq(20, ncol(intersect_programs_hc), by=20)], labels= seq(20, ncol(intersect_programs_hc), by=20)) + 
  scale_y_discrete(name="", breaks=colnames(intersect_programs_hc)[seq(20, ncol(intersect_programs_hc), by=20)], labels= seq(20, ncol(intersect_programs_hc), by=20)) +
  panel_border(colour = "black", size = 0.4, linetype = 1, remove = FALSE)  +
  guides(fill = guide_colourbar(barheight = 0.8, barwidth = 4.8, title.position = "top", title.hjust = 0.5))

# cancer type annotation
ccle_annotation <- readRDS(file.choose())
ccle_annotation$cancer_type[is.element(ccle_annotation$cancer_type, c("Fibroblast", "Gallbladder Cancer", "Neuroblastoma", "Prostate Cancer", "Bone Cancer", "Sarcoma", "Bile Duct Cancer", "Thyroid Cancer"))] <- "Others"

color_match <- data.frame("type"=sort(unique(ccle_annotation$cancer_type)), "color"= c(brewer.pal(12, "Set3")[c(1:6,8,7,10:11)], "maroon","gray93", "yellow2", "goldenrod1", "slateblue2"),stringsAsFactors = F)

annot_cancer_type <- data.frame("cell_lines" =sub("\\..*", "", rownames(intersect_programs_hc)), stringsAsFactors = F)
annot_cancer_type$type <- ccle_annotation$cancer_type[match(annot_cancer_type$cell_lines, rownames(ccle_annotation))]

p6 <- ggplot(annot_cancer_type, aes(y=0, x=1:nrow(annot_cancer_type), fill=type, color=type)) +
  geom_tile() +
  scale_color_manual(breaks=c(sort(unique(annot_cancer_type$type))[sort(unique(annot_cancer_type$type))!="Others"], "Others"), values = color_match$color[match(sort(unique(annot_cancer_type$type)), color_match$type)] , name="") +
  scale_fill_manual(breaks=c(sort(unique(annot_cancer_type$type))[sort(unique(annot_cancer_type$type))!="Others"], "Others"), values = color_match$color[match(sort(unique(annot_cancer_type$type)), color_match$type)] , name="") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 11),  axis.title = element_text(size=8),  plot.margin = unit(c(3.5,4,-2.5,1), "cm"), legend.position = c(0,4), legend.text.align = 0, legend.key.size = unit(0.4, "cm")) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x="", y="") +
  guides(fill = guide_legend(title = "",  override.aes = list(size = 5), ncol=3), color=FALSE) 

egg::ggarrange(p6,p5, nrow=2, heights = c(0.5,9))

pdf("discrete_programs_eps1.2.pdf", height = 7, width = 7, onefile = F)

# color legend
ggplot(annot_cancer_type, aes(y=0, x=1:nrow(annot_cancer_type), fill=type, color=type)) +
  geom_tile() +
  scale_color_manual(breaks=c(sort(unique(annot_cancer_type$type))[sort(unique(annot_cancer_type$type))!="Others"], "Others"), values = color_match$color[match(sort(unique(annot_cancer_type$type)), color_match$type)] , name="") +
  scale_fill_manual(breaks=c(sort(unique(annot_cancer_type$type))[sort(unique(annot_cancer_type$type))!="Others"], "Others"), values = color_match$color[match(sort(unique(annot_cancer_type$type)), color_match$type)] , name="") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 10),  axis.title = element_text(size=8), legend.text.align = 0, legend.key.size = unit(0.4, "cm")) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x="", y="") +
  guides(fill = guide_legend(title = "",  override.aes = list(size = 5), nrow=3), color=FALSE) 

pdf("discrete_programs_heatmap_legend.pdf", height = 2, width = 12, onefile = F)

cluster_1 <- table(unlist(discr_programs_sig[colnames(intersect_programs_hc)[83:95]]))
cluster_1 <- sort(cluster_1, decreasing = T)[1:30]

cluster_2 <- table(unlist(discr_programs_sig[colnames(intersect_programs_hc)[106:119]]))
cluster_2 <- sort(cluster_2, decreasing = T)[1:30]

cluster_3 <- table(unlist(discr_programs_sig[colnames(intersect_programs_hc)[123:141]]))
cluster_3 <- sort(cluster_3, decreasing = T)[1:30]



########################### comparing number of clusters ############################
discr_eps1.2 <- readRDS(file.choose())
discr_eps1.5 <- readRDS(file.choose())
discr_eps1.8 <- readRDS(file.choose())


three_more <- c(length(which(sapply(discr_eps1.2, function(x) length(x$clusters_cells)>2))), length(which(sapply(discr_eps1.5, function(x) length(x$clusters_cells)>2))), length(which(sapply(discr_eps1.8, function(x) length(x$clusters_cells)>2))))

one_large <- c(length(which(sapply(discr_eps1.2, function(x)  (min(sapply(x$clusters_cells, length))/max(sapply(x$clusters_cells, length))) > 0.2 ))), length(which(sapply(discr_eps1.5, function(x)  (min(sapply(x$clusters_cells, length))/max(sapply(x$clusters_cells, length))) > 0.2 ))), length(which(sapply(discr_eps1.8, function(x)  (min(sapply(x$clusters_cells, length))/max(sapply(x$clusters_cells, length))) > 0.2 ))))

one_small <- c(length(discr_eps1.2), length(discr_eps1.5), length(discr_eps1.8)) - (three_more + one_large)

none <- c(198-49, 198-30, 198-22)

pie_chart <- rbind(three_more, one_large, one_small, none)
pie_chart <- 100*pie_chart/colSums(pie_chart)
colnames(pie_chart) <- c("eps1.2", "eps1.5", "eps1.8")
  

ggplot(data.frame(pie_chart), aes(x="", y=eps1.2, fill=rownames(pie_chart))) +
  geom_bar(stat="identity", alpha=0.8) +
  coord_polar("y") +
  scale_fill_brewer(palette = "Pastel2", name="", breaks= c("none", "one_small", "one_large", "three_more"), labels=c("None", "1 small", "1 large", "3 or more")) +
  labs(x="", y="") +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), legend.text = element_text(size=12))

pdf("pie_discr_eps1.2.pdf", width = 5, height = 3 )

ggplot(data.frame(pie_chart), aes(x="", y=eps1.5, fill=rownames(pie_chart))) +
  geom_bar(stat="identity", alpha=0.8) +
  coord_polar("y") +
  scale_fill_brewer(palette = "Pastel2", name="", breaks= c("none", "one_small", "one_large", "three_more"), labels=c("None", "1 small", "1 large", "3 or more")) +
  labs(x="", y="") +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), legend.text = element_text(size=12))

pdf("pie_discr_eps1.5.pdf", width = 5, height = 3 )

ggplot(data.frame(pie_chart), aes(x="", y=eps1.8, fill=rownames(pie_chart))) +
  geom_bar(stat="identity", alpha=0.8) +
  coord_polar("y") +
  scale_fill_brewer(palette = "Pastel2", name="", breaks= c("none", "one_small", "one_large", "three_more"), labels=c("None", "1 small", "1 large", "3 or more")) +
  labs(x="", y="") +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), legend.text = element_text(size=12))

pdf("pie_discr_eps1.8.pdf", width = 5, height = 3 )


length(which(sapply(discr_eps1.2, function(x) length(x$clusters_cells)==2 & (min(sapply(x$clusters_cells, length))/max(sapply(x$clusters_cells, length))) > 0.8 )))



length(which(sapply(discr_eps1.2, function(x) length(x$clusters_cells)==2 & (min(sapply(x$clusters_cells, length))/max(sapply(x$clusters_cells, length))) > 0.8 )))

length(which(sapply(discr_eps1.8, function(x)  (min(sapply(x$clusters_cells, length))/max(sapply(x$clusters_cells, length))) > 0.2 )))

length(which(sapply(discr_eps1.2, function(x) length(x$clusters_cells)==2)))
dd

length(which(sapply(discr_eps1.8, function(x) length(x$clusters_cells)==2 & (min(sapply(x$clusters_cells, length))/max(sapply(x$clusters_cells, length))) > 0.8 )))














p1 <- ggplot(data = plot_hc, aes(x=Var1, y=Var2, fill=value, color=value)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(2,20), low=gabi_palette[1:111],  mid =gabi_palette[112:222], high = gabi_palette[223:333], midpoint = 11, oob=squish, name="Similarity")+
  scale_color_gradient2(limits=c(2,20), low=gabi_palette[1:111],  mid =gabi_palette[112:222], high = gabi_palette[223:333], midpoint = 11, oob=squish, name="Similarity")+
  theme( axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(),  axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification ="bottom") + 
  scale_x_discrete(name="\nPrograms", breaks=colnames(intersect_programs_hc)[seq(3, ncol(intersect_programs_hc), by=5)], labels= seq(3, ncol(intersect_programs_hc), by=5)) + 
  scale_y_discrete(name="", breaks=colnames(intersect_programs_hc)[seq(3, ncol(intersect_programs_hc), by=5)], labels= seq(3, ncol(intersect_programs_hc), by=5)) +
  panel_border(colour = "black", size = 0.4, linetype = 1, remove = FALSE)  
  



### plotting cell line annotation
annot <- data.frame("cell_lines" = sub("\\..*", "", colnames(intersect_programs_hc)), stringsAsFactors = F)
annot$lineage <- gsub("^[^_]*_", "", annot$cell_lines)
annot$lineage <- gsub("_", " ",annot$lineage)
annot$cell_lines <- sub('_', ' - ', annot$cell_lines)
annot$cell_lines <- gsub("_", " ", annot$cell_lines)
annot$cell_lines <- factor(annot$cell_lines, levels = unique(annot$cell_lines[order(annot$lineage, annot$cell_lines)]))


p2 <- ggplot(annot, aes(x=as.factor(1:nrow(annot)), y = 1)) + 
  geom_tile(aes(fill=lineage, colour=lineage)) + 
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 8),  legend.box.margin = unit(c(10,1,5,0), "cm"), plot.margin = unit(c(0,0,-0.7,0), "cm"),  legend.position = "left", legend.text.align = 1, legend.justification = "right", legend.key.size = unit(0.4, "cm")) + 
  labs(x="", y="") + 
  guides(fill = guide_legend(title = "", label.position="left",  override.aes = list(linetype = 0, size = 0.4), ncol = 1), color=FALSE)  + 
  scale_fill_brewer(palette = "Set3") +
  scale_color_brewer(palette = "Set3") 


### plotting correlation with  cnv clones
cnv_corr <- data.frame("program" = colnames(intersect_programs_hc), "cnv"=NA)
cnv_corr$cnv[is.element(cnv_corr$program, names(discr_gmm))] <- "cnv"
cnv_corr$cnv[!is.element(cnv_corr$program, names(discr_gmm))] <- "no_cnv"

p3 <- ggplot(cnv_corr, aes(x=as.factor(1:nrow(cnv_corr)), y = 1)) + 
  geom_tile(aes(fill=cnv, colour=cnv)) + 
  scale_fill_manual(values = c("darkblue", "gray97"), breaks="cnv", labels="  CNV clone-derived") +
  scale_colour_manual(values = c("darkblue", "gray97")) +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 10), legend.title = element_blank(), plot.margin = unit(c(1,0,-0.7,0), "cm")) + 
  labs(x="", y="") +
  guides(colour=F)


### combining
egg::ggarrange(p2, p1, nrow = 2, heights = c(0.7,8))

pdf("discrete_programs.pdf", height = 5.7, width = 9, onefile=FALSE)













### reading gmm test - cnv clones
gmm_final <- readRDS(file.choose())

### reading inferred cnv 
cnv_smooth_infer <- readRDS(file.choose())
cnv_smooth_infer_arm <- readRDS(file.choose())
for(i in names(cnv_smooth_infer_arm)) { # correcting cell names
  colnames(cnv_smooth_infer_arm[[i]]) <- gsub("[.]", "-" ,colnames(cnv_smooth_infer_arm[[i]]))
  colnames(cnv_smooth_infer_arm[[i]]) <- gsub("X", "" ,colnames(cnv_smooth_infer_arm[[i]]))
}

### reading final discrete programs and their signature genes
discr_programs <- readRDS(file.choose()) # cluster cell composition
discr_programs_sig <- readRDS(file.choose()) # cluster signature

### association between discrete programs and cnv clones
discr_gmm <- list()

for(i in names(discr_programs)) {
  a <- rownames(cnv_smooth_infer[[sub("\\..*", "", i)]]) # gets all cells from the selected cell line
  b <- as.character(is.element(a,discr_programs[[i]])) # checks which cells are part of the discrete cluster
  c <- gmm_final[[sub("\\..*", "", i)]] # gets gmm test for the selected cell line for each chromosome arm
  d <- c()
  for(j in names(c)){ 
    e <- data.frame("clone" = c[[j]][,"class"], "discr" = b) # combines gmm test results and discrete program classification
    d[j] <- fisher.test(table(e))$p.value # fisher test - classification vs discrete program classification
  }
  discr_gmm[[i]] <- d
}

discr_gmm <- discr_gmm[sapply(discr_gmm, function(x) length(which(x<0.05)) != 0)]












### plotting correlation with cell complexity
complexity <- lapply(expr, function(x) apply(x,2, function(y) length(which(y!=0))))

corr_complexity <- data.frame("program" = colnames(intersect_programs_hc), "corr"=NA, row.names = colnames(intersect_programs_hc))

for(i in colnames(intersect_programs_hc)) {
  a <- complexity[[sub("\\..*", "", i)]]
  b <- discr_programs[[i]]
  corr_complexity[i,"corr"] <- log2(mean(a[is.element(names(a), b)])/mean(a[!is.element(names(a), b)]))
}


corr_complexity <- data.frame("program" = colnames(intersect_programs_hc), "corr"=NA)
for(i in 1:ncol(intersect_programs_hc)) {
  a <- colnames(intersect_programs_hc)[i]
  b <- complexity[[sub("\\..*", "", a)]]
  c <- colMeans(expr[[sub("\\..*", "", a)]][discr_programs_sig[,a][!is.na(discr_programs_sig[,a])],])
  corr_complexity[i,"corr"] <- cor(b,c)
}

p2 <- ggplot(corr_complexity, aes(y=corr, x = 1:nrow(corr_complexity))) +
  geom_smooth(span=0.2, se = FALSE, color="gray36", size= 0.8, method = "loess") + 
  geom_point(alpha=0.3, size = 1.5) + 
  theme(axis.ticks = element_blank(), axis.line = element_blank(), panel.background = element_rect(fill = "gray97"), panel.grid = element_blank(), axis.text.x = element_blank(), axis.title = element_text(size = 12), axis.text.y = element_text(size = 10), plot.margin=unit(c(1,1,-0.6,1), "cm") ) + 
  labs(y="Cell complexity\ncorrelation", x="") + 
  scale_y_continuous(expand = c(0,0), limits = c(-0.6,0.6), breaks = seq(-0.5, 0.5, 0.1)) + 
  scale_x_continuous(expand = c(0,0))




### plotting lineage abundance
p3 <- ggplot(data.frame()) + 
  geom_blank() + 
  scale_x_continuous(limits = c(0, ncol(intersect_programs_hc)), expand = c(0,0)) +
  scale_y_continuous(limits = c(-3, 25), expand = c(0,0), breaks = seq(0, 20, 5)) +
  geom_segment(aes(x = 16, y = 8, xend = 25, yend = 8), size=1) +
  labs(x="", y="Lineages") +
  theme(axis.ticks = element_blank(), axis.line = element_blank(), panel.background = element_rect(fill = "gray97"), panel.grid = element_blank(), axis.text.x = element_blank(), axis.title = element_text(size = 12), axis.text.y = element_text(size = 10), plot.margin = unit(c(1,1,-0.8,1), "cm")) 

egg::ggarrange(p3, p2, p1, nrow = 3, heights = c(10,2,50))


### ploting cnv signal and cnv clones
cnv_ave_signal <- sapply(cnv_smooth_infer, function(x) mean(x^2))
cnv_ave_signal_plot <- data.frame("cnv_ave_signal"=cnv_ave_signal)
cnv_ave_signal_plot$clone_status <- ""
cnv_ave_signal_plot$clone_status[sapply(cnv_clones, function(x) length(which(x<0.05)) != 0)] <- "CNV clone"
cnv_ave_signal_plot <-  cnv_ave_signal_plot[order(cnv_ave_signal_plot$cnv_ave_signal),]

ggplot(cnv_ave_signal_plot, aes(x=1:nrow(cnv_ave_signal_plot), y=cnv_ave_signal, color=clone_status)) +
  geom_point() +
  scale_color_manual(values=c("gray80", "coral2"), breaks="CNV clone") +
  theme(legend.title = element_blank())




permutation_discr <- list()

for(i in "SCC9_UPPER_AERODIGESTIVE_TRACT.1.1") {
  a <- discr_programs[[i]]
  b <- cnv_smooth_infer[[sub("\\..*", "", i)]]
  c <- is.element(rownames(b), a)
  d <- c()
  for(j in 1:10) {
    e <- b[sample(rownames(b)),]
    d[j] <- dist(rbind(colMeans(e[c,]),colMeans(e[!c,]))) # euclidean distance: sqrt(sum((e[1,]-e[2,])^2))
  }
  percentile <- ecdf(d)
  f <- dist(rbind(colMeans(b[c,]),colMeans(b[!c,])))
  permutation_discr[[i]] <- list("disrtribution"= d, "p"= 1-percentile(f), "z" = (f - mean(d))/sd(d))
}


