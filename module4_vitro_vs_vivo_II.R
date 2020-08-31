# ----------------------------------------------------------------------------------------------------
# Module 4. Evaluating expression heterogeneity in HNSCC and melanoma by analyzing human tumors and cell lines as a single dataset
# ----------------------------------------------------------------------------------------------------

# **************************************************************************
# Basic setup

# load necessary R packages and functions
library(reshape2)
library(ggplot2)
library(scales)
                                               
# read scRNA-seq data from cell lines and tumors
expr_ccle <- readRDS("CCLE_heterogeneity_Rfiles/CCLE_scRNAseq_CPM.RDS") # CCLE cell lines
expr_tumor <- readRDS("CCLE_heterogeneity_Rfiles/tumors_scRNAseq_logTPM.RDS") # human tumors

# calculate average gene expression
rowmeans_ccle <- lapply(expr_ccle, rowMeans)                 
rowmeans_tumor <- lapply(expr_tumor, function(x) rowMeans(10*((2^x)-1)))  
                          
# process data                            
expr_ccle<- lapply(expr_ccle, function(x) {log2((x/10) + 1)})                        
expr_ccle <- lapply(expr_ccle, function(x) {x-rowMeans(x)})
 
expr_tumor <- lapply(expr_tumor, function(x) {x-rowMeans(x)})

# select relevant cell lines and tumors   
nmf_metaprograms_programs_nc_ccle<- readRDS("Expected_results/module2/nmf_metaprograms_programs_nc_ccle.RDS")   # programs in each cell line metaprogram
nmf_metaprograms_programs_nc_tumor<- readRDS("Expected_results/module2/nmf_metaprograms_programs_nc_tumor.RDS") # programs in each tumor metaprogram                             
                          
episen_emt_ccle <- intersect(gsub(".{4}$", "",nmf_metaprograms_programs_nc_ccle$episen), gsub(".{4}$", "",nmf_metaprograms_programs_nc_ccle$emtII))
episen_emt_ccle <- unique(episen_emt_ccle[grep("UPPER", episen_emt_ccle)])  # HNCSS cell lines harboring both EpiSen and EMTII                                                                               
episen_emt_tumor <- unique(intersect(gsub(".{4}$", "",nmf_metaprograms_programs_nc_tumor$episen), gsub(".{4}$", "",nmf_metaprograms_programs_nc_tumor$emt))) # HNCSS tumors harboring both EpiSen and EMTII    

skinpig_emt_ccle <- intersect(gsub(".{4}$", "",nmf_metaprograms_programs_nc_ccle$skinpig), gsub(".{4}$", "",nmf_metaprograms_programs_nc_ccle$emtI)) # melanoma cell lines harboring both skinpig and EMTI      
                                  
skinpig_emt_tumor <- names(expr_tumor)[grep("mel", names(expr_tumor))] # melanoma tumors harboring both skinpig and EMTI  

# read metaprogram signatures from cell lines
meta_sig_ccle <- readRDS("Expected_results/module2/nmf_metaprograms_sig_nc_ccle.RDS")
meta_sig_ccle <- meta_sig_ccle[c("skinpig", "emtI", "episen", "emtII")]
names(meta_sig_ccle) <- paste0(names(meta_sig_ccle), "_vitro")                         
meta_sig_ccle <- lapply(meta_sig_ccle, function(x) names(x[1:100]))           
                                                                                    
# read metaprogram signatures from tumors                                             
meta_sig_tumor  <- unlist(list(read.table("CCLE_heterogeneity_Rfiles/metaprograms_tumors_literature.txt", sep = "\t", header = T,stringsAsFactors = F)), recursive=F)
meta_sig_tumor <- meta_sig_tumor[c("melanoma.MITF.program", "melanoma.AXL.program", "HNSCC.PEMT","HNSCC.Epidif.1" )]  
meta_sig_tumor <- lapply(meta_sig_tumor, function(x) x[x!=""])                                             
  
# get genes shared between programs in tumors and cell lines                
meta_intersect <- lapply(meta_sig_ccle, function(x) lapply(meta_sig_tumor, function(y) intersect(x,y)))

# **************************************************************************
# Combine datasets and run PCA - HNSCC
                                                         
# get expression data from selected cell lines
expr_emt_episen_vitro <- expr_ccle[episen_emt_ccle]
expr_emt_episen_vivo <- expr_tumor[episen_emt_tumor]                                         
                                                           
# select common genes in the cell line and tumor datasets             
expr_emt_episen_vitro <- sapply(expr_emt_episen_vitro, function(x) x[Reduce(intersect, c(lapply( expr_emt_episen_vitro, rownames), lapply(expr_emt_episen_vivo, rownames))),])
expr_emt_episen_vivo <- sapply(expr_emt_episen_vivo, function(x) x[Reduce(intersect, c(lapply(expr_emt_episen_vitro, rownames), lapply(expr_emt_episen_vivo, rownames))),]) 
 
# unlist datasets                               
expr_emt_episen_vitro <- do.call(cbind,expr_emt_episen_vitro )        
expr_emt_episen_vivo <-  do.call(cbind, expr_emt_episen_vivo)                  
                                                         
# select top expressed genes among cell lines and tumors      
rowmeans_ccle_hnscc <-rowMeans(sapply(rowmeans_ccle[episen_emt_ccle], function(x) x[rownames(expr_emt_episen_vitro)]))                                           
rowmeans_tumor_hnscc <- rowMeans(sapply(rowmeans_tumor[episen_emt_tumor], function(x) x[rownames(expr_emt_episen_vitro)]))
                
expr_emt_episen_vitro_filt <- expr_emt_episen_vitro[order(rowMeans(cbind(rowmeans_ccle_hnscc, rowmeans_tumor_hnscc)), decreasing=T)[1:4500] ,]
expr_emt_episen_vivo_filt <- expr_emt_episen_vivo[order(rowMeans(cbind(rowmeans_ccle_hnscc, rowmeans_tumor_hnscc)), decreasing=T)[1:4500]  ,]
 
# combine the two datasets                                        
expr_emt_episen_final <- cbind(expr_emt_episen_vitro_filt,expr_emt_episen_vivo_filt)    

# run pca
expr_emt_episen_pca <- prcomp(t(expr_emt_episen_final))                       
                               

# **************************************************************************
# Combine datasets and run PCA - melanoma
                                                         
# get expression data from selected cell lines
expr_skinpig_emt_vitro <- expr_ccle[skinpig_emt_ccle]
expr_skinpig_emt_vivo <- expr_tumor[skinpig_emt_tumor]                
                                        
# select common genes in the cell line and tumor datasets             
expr_skinpig_emt_vitro <- sapply(expr_skinpig_emt_vitro, function(x) x[Reduce(intersect, c(lapply( expr_skinpig_emt_vitro, rownames), lapply(expr_skinpig_emt_vivo, rownames))),])
expr_skinpig_emt_vivo <- sapply(expr_skinpig_emt_vivo, function(x) x[Reduce(intersect, c(lapply(expr_skinpig_emt_vitro, rownames), lapply(expr_skinpig_emt_vivo, rownames))),])    
                              
# unlist datasets                               
expr_skinpig_emt_vitro <- do.call(cbind, expr_skinpig_emt_vitro)        
expr_skinpig_emt_vivo <- do.call(cbind, expr_skinpig_emt_vivo)                               
                                                         
# select top expressed genes among cell lines and tumors      
rowmeans_ccle_melanoma <-rowMeans(sapply(rowmeans_ccle[skinpig_emt_ccle], function(x) x[rownames(expr_skinpig_emt_vitro)]))                                        
rowmeans_tumor_melanoma<- rowMeans(sapply(rowmeans_tumor[skinpig_emt_tumor], function(x) x[rownames(expr_skinpig_emt_vitro)]))
 
expr_skinpig_emt_vitro_filt <- expr_skinpig_emt_vitro[order(rowMeans(cbind(rowmeans_ccle_melanoma, rowmeans_tumor_melanoma)), decreasing=T)[1:4500] ,]
expr_skinpig_emt_vivo_filt <- expr_skinpig_emt_vivo[order(rowMeans(cbind(rowmeans_ccle_melanoma, rowmeans_tumor_melanoma)), decreasing=T)[1:4500] ,]

# combine the two datasets                                                       
expr_skinpig_emt_final <- cbind(expr_skinpig_emt_vitro_filt, expr_skinpig_emt_vivo_filt) 
                                          
# run PCA                                          
expr_skinpig_emt_pca <- prcomp(t(expr_skinpig_emt_final))                                          
                                                                                                                   
# **************************************************************************
# Plot PCA results - HNSCC
                                                                                            
# store pca coordinates in  dataframe                        
expr_emt_episen_plot <- data.frame(expr_emt_episen_pca$x[,1:5])
 
# calculate program scores 
expr_emt_episen_plot$emtII <- colMeans(expr_emt_episen_final[is.element(rownames(expr_emt_episen_final), meta_intersect$emtII_vitro$HNSCC.PEMT),])
                                          
expr_emt_episen_plot$episen=colMeans(expr_emt_episen_final[is.element(rownames(expr_emt_episen_final), meta_intersect$episen_vitro$`HNSCC.Epidif.1`),])    

# add  metadata     
expr_emt_episen_plot$type <- rep(c("in vitro", "in vivo"), c(ncol(expr_emt_episen_vitro), ncol(expr_emt_episen_vivo)))                                          
expr_emt_episen_plot$names <- c(rep(episen_emt_ccle, sapply(expr_ccle[episen_emt_ccle],ncol)),rep(episen_emt_tumor, sapply(expr_tumor[episen_emt_tumor],ncol)) )
expr_emt_episen_plot$cell_line <- NA
expr_emt_episen_plot$cell_line[expr_emt_episen_plot$type=="in vitro"] <- expr_emt_episen_plot$names[expr_emt_episen_plot$type=="in vitro"] 
expr_emt_episen_plot$cell_line <-   gsub("_UPPER_AERODIGESTIVE_TRACT", "", expr_emt_episen_plot$cell_line)   
expr_emt_episen_plot$tumor <- NA
expr_emt_episen_plot$tumor[expr_emt_episen_plot$type=="in vivo"] <- expr_emt_episen_plot$names[expr_emt_episen_plot$type=="in vivo"]      
                                                                                            
# plot pca coordinates     
pdf("Output/module4/HNSCC_pca_program.pdf", width = 3.9, height = 3.8)  
ggplot(expr_emt_episen_plot, aes(x=PC2, y=PC4+PC5, color=emtII-episen)) +
            geom_point( size=0.5) +
            scale_color_gradient2(limits=c(-3, 3), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(EMT II - Episen)") +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))
dev.off()
                                                                                                                                         
pdf("Output/module4/HNSCC_pca_type.pdf", width = 3.95, height = 3.7)  
ggplot(expr_emt_episen_plot, aes(x=PC2, y= PC4+PC5, color=type)) +
            geom_point( size=0.5) +
             scale_color_manual(values=c("goldenrod", "mediumpurple"), name="") +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
            guides(color=guide_legend(override.aes = list(size=4)))  
dev.off()
                                                                                                                                         
pdf("Output/module4/HNSCC_pca_cellline.pdf", width = 3.9, height = 3.8)  
ggplot(expr_emt_episen_plot, aes(x=PC2, y= PC4 + PC5, color=cell_line)) +
            geom_point(data=expr_emt_episen_plot[is.na(expr_emt_episen_plot$cell_line),], aes(x=PC2, y= PC4+PC5), size=0.5, color="gray90") + 
             geom_point( size=0.5) +
             scale_color_brewer(palette = "Set2", name="", breaks=unique(expr_emt_episen_plot$cell_line)[!is.na(unique(expr_emt_episen_plot$cell_line))]) +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
            guides(color=guide_legend(override.aes = list(size=3), nrow=2) )           
dev.off()
                                                                                                       
pdf("Output/module4/HNSCC_pca_tumor.pdf", width = 3.9, height = 3.8) 
ggplot(expr_emt_episen_plot, aes(x=PC2, y= PC4 + PC5, color=tumor)) +
            geom_point( size=0.5) +
             scale_color_brewer(palette = "Set2", name="", na.value="gray90", breaks=unique(expr_emt_episen_plot$tumor)[!is.na(unique(expr_emt_episen_plot$tumor))], labels=c("HNSCC 18","HNSCC 26", "HNSCC 25", "HNSCC 5", "HNSCC 17" )) +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
            guides(color=guide_legend(override.aes = list(size=3), nrow=2))                                              
dev.off()                                   
                                                                                             # plot top genes with positive pca loadings                                                 
pca_loadings_pos_hnscc <- rev(melt(apply(expr_emt_episen_pca$rotation[,1:5], 2, function(x) names(sort(x, decreasing=T))[1:10]))$value)                                               
pca_loadings_pos_hnscc_plot <- as.matrix(expr_emt_episen_pca$rotation[as.character(pca_loadings_pos_hnscc),1:5] )
rownames(pca_loadings_pos_hnscc_plot) <- NULL
pca_loadings_pos_hnscc_plot <- melt(as.matrix(pca_loadings_pos_hnscc_plot))
                                           
pdf("Output/module4/HNSCC_pca_pos.pdf", width = 4.5, height = 7.5)                                     
ggplot(data = pca_loadings_pos_hnscc_plot, aes(x=Var2, y=factor(Var1), fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(-0.12, 0.12), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="PC loadings")+
  scale_color_gradient2(limits=c(-0.12, 0.12), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="PC loadings") +
  scale_y_discrete(labels=pca_loadings_pos_hnscc)  +     
  geom_hline(yintercept = seq(10.5,40.5,10)) + 
  geom_vline(xintercept = c(1.5,2.5,3.5, 4.5))       +                              
  scale_x_discrete(expand=c(0,0)) +                                   
  labs(y="", x="") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), axis.line = element_blank(), legend.title= element_text(size=11), legend.text = element_text(size=10), axis.text = element_text(size=11), legend.text.align = 1, axis.text.y = element_text(face= ifelse(is.element(pca_loadings_pos_hnscc, c(meta_intersect$emtII_vitro$HNSCC.PEMT,meta_intersect$episen_vitro$`HNSCC.Epidif.1`)), "bold","plain")))          
 dev.off()
                                         
# plot top genes with negative pca loadings                            
pca_loadings_neg_hnscc <- rev(melt(apply(expr_emt_episen_pca$rotation[,1:5], 2, function(x) names(sort(x, decreasing=F))[1:10]))$value)                                               
pca_loadings_neg_hnscc_plot <- as.matrix(expr_emt_episen_pca$rotation[as.character(pca_loadings_neg_hnscc),1:5] )
rownames(pca_loadings_neg_hnscc_plot) <- NULL
pca_loadings_neg_hnscc_plot <- melt(as.matrix(pca_loadings_neg_hnscc_plot))

pdf("Output/module4/HNSCC_pca_neg.pdf", width = 4.5, height = 7.5)                                     
ggplot(data = pca_loadings_neg_hnscc_plot, aes(x=Var2, y=factor(Var1), fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(-0.12, 0.12), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="PC loadings")+
  scale_color_gradient2(limits=c(-0.12, 0.12), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="PC loadings") +
  scale_y_discrete(labels=pca_loadings_neg_hnscc)  +     
  geom_hline(yintercept = seq(10.5,40.5,10)) + 
  geom_vline(xintercept = c(1.5,2.5,3.5, 4.5))       +                              
  scale_x_discrete(expand=c(0,0)) +                                   
  labs(y="", x="") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), axis.line = element_blank(), legend.title= element_text(size=11), legend.text = element_text(size=10), axis.text = element_text(size=11), legend.text.align = 1, axis.text.y = element_text(face= ifelse(is.element(pca_loadings_neg_hnscc, c(meta_intersect$emtII_vitro$HNSCC.PEMT,meta_intersect$episen_vitro$`HNSCC.Epidif.1`)), "bold","plain")))  
dev.off()

# **************************************************************************
# Plot PCA results - melanoma
                                                                                            
# store pca coordinates in a dataframe
expr_skinpig_emt_plot <- data.frame(expr_skinpig_emt_pca$x[,1:5] )     
                                         
# calculate program scores 
expr_skinpig_emt_plot$emtI <- colMeans(expr_skinpig_emt_final[is.element(rownames(expr_skinpig_emt_final), meta_intersect$emtI_vitro$melanoma.AXL.program),])                                        
expr_skinpig_emt_plot$skinpig <- colMeans(expr_skinpig_emt_final[is.element(rownames(expr_skinpig_emt_final), meta_intersect$skinpig_vitro$melanoma.MITF.program),])    

# add additional metadata         
expr_skinpig_emt_plot$type <- rep(c("in vitro", "in vivo"), c(ncol(expr_skinpig_emt_vitro), ncol(expr_skinpig_emt_vivo)))                                  
expr_skinpig_emt_plot$names <- c(rep(skinpig_emt_ccle, sapply(expr_ccle[skinpig_emt_ccle],ncol)),rep(skinpig_emt_tumor, sapply(expr_tumor[skinpig_emt_tumor],ncol)) )
expr_skinpig_emt_plot$cell_line <- NA
expr_skinpig_emt_plot$cell_line[expr_skinpig_emt_plot$type=="in vitro"] <- expr_skinpig_emt_plot$names[expr_skinpig_emt_plot$type=="in vitro"] 
expr_skinpig_emt_plot$cell_line <-   gsub("_SKIN", "", expr_skinpig_emt_plot$cell_line)   
expr_skinpig_emt_plot$tumor <- NA
expr_skinpig_emt_plot$tumor[expr_skinpig_emt_plot$type=="in vivo"] <- expr_skinpig_emt_plot$names[expr_skinpig_emt_plot$type=="in vivo"]                                                     
# plot pca coordinates 
                                         
pdf("Output/module4/melanoma_pca_program.pdf", width = 3.9, height = 3.8)  
ggplot(expr_skinpig_emt_plot, aes(x=PC2+PC3, y=PC4, color=emtI-skinpig)) +
            geom_point( size=0.5) +
            scale_color_gradient2(limits=c(-3, 3), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(EMT I - Skinpig)") +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))
dev.off()
                                                                                                                                         
pdf("Output/module4/melanoma_pca_type.pdf", width = 3.95, height = 3.7)  
ggplot(expr_skinpig_emt_plot, aes(x=PC2+PC3, y=PC4, color=type)) +
            geom_point( size=0.5) +
             scale_color_manual(values=c("goldenrod", "mediumpurple"), name="") +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
            guides(color=guide_legend(override.aes = list(size=4)))  
dev.off()
                                                                                                                                         
pdf("Output/module4/melanoma_pca_cellline.pdf", width = 3.9, height = 3.8)  
ggplot(expr_skinpig_emt_plot, aes(x=PC2+PC3, y= PC4, color=cell_line)) +
            geom_point(data=expr_skinpig_emt_plot[is.na(expr_skinpig_emt_plot$cell_line),], aes(x=PC2+PC3, y= PC4), size=0.5, color="gray90") + 
             geom_point( size=0.5) +
             scale_color_brewer(palette = "Set2", name="", breaks=unique(expr_skinpig_emt_plot$cell_line)[!is.na(unique(expr_skinpig_emt_plot$cell_line))]) +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=10.4), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
            guides(color=guide_legend(override.aes = list(size=3), nrow=2) )    
dev.off()
                                                                                                       
pdf("Output/module4/melanoma_pca_tumor.pdf", width = 3.9, height = 3.8) 
ggplot(expr_skinpig_emt_plot, aes(x=PC2+PC3, y= PC4, color=tumor)) +
            geom_point( size=0.5) +
             scale_color_brewer(palette = "Set2", name="", na.value="gray90", breaks=unique(expr_skinpig_emt_plot$tumor)[!is.na(unique(expr_skinpig_emt_plot$tumor))], labels=c("Mel 59","Mel 71", "Mel 78", "Mel 19", "Mel 80", "Mel 81", "Mel 88", "Mel 89" )) +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
            guides(color=guide_legend(override.aes = list(size=3), nrow=2))                                        
dev.off()                                   
                                                                                             # plot top genes with positive pca loadings                                                 
pca_loadings_pos_mel <- rev(melt(apply(expr_skinpig_emt_pca$rotation[,1:5], 2, function(x) names(sort(x, decreasing=T))[1:10]))$value)              
pca_loadings_pos_mel_plot <- as.matrix(expr_mitf_emt_pca$rotation[as.character(pca_loadings_pos_mel),1:5] )
rownames(pca_loadings_pos_mel_plot) <- NULL
pca_loadings_pos_mel_plot <- melt(as.matrix(pca_loadings_pos_mel_plot))

                                           
pdf("Output/module4/melanoma_pca_pos.pdf", width = 4.5, height = 7.5)                                     
ggplot(data = pca_loadings_pos_mel_plot, aes(x=Var2, y=factor(Var1), fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(-0.12, 0.12), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="PC loadings")+
  scale_color_gradient2(limits=c(-0.12, 0.12), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="PC loadings") +
  scale_y_discrete(labels=pca_loadings_pos_mel)  +     
  geom_hline(yintercept = seq(10.5,40.5,10)) + 
  geom_vline(xintercept = c(1.5,2.5,3.5, 4.5))       +                              
  scale_x_discrete(expand=c(0,0)) +                                   
  labs(y="", x="") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), axis.line = element_blank(), legend.title= element_text(size=11), legend.text = element_text(size=10), axis.text = element_text(size=11), legend.text.align = 1, axis.text.y = element_text(face= ifelse(is.element(pca_loadings_pos_mel, c(meta_intersect$emtI_vitro$melanoma.AXL.program,meta_intersect$skinpig_vitro$melanoma.MITF.program)), "bold","plain")))  
 dev.off()
                                         
# plot top genes with negative pca loadings                            
pca_loadings_neg_mel <- rev(melt(apply(expr_skinpig_emt_pca$rotation[,1:5], 2, function(x) names(sort(x, decreasing=F))[1:10]))$value)                                              
pca_loadings_neg_mel_plot <- as.matrix(expr_skinpig_emt_pca$rotation[as.character(pca_loadings_neg_mel),1:5] )
rownames(pca_loadings_neg_mel_plot) <- NULL
pca_loadings_neg_mel_plot <- melt(as.matrix(pca_loadings_neg_mel_plot))

pdf("Output/module4/melanoma_pca_neg.pdf", width = 4.5, height = 7.5)                       
ggplot(data = pca_loadings_neg_mel_plot, aes(x=Var2, y=factor(Var1), fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(-0.12, 0.12), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="PC loadings")+
  scale_color_gradient2(limits=c(-0.12, 0.12), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="PC loadings") +
  scale_y_discrete(labels=pca_loadings_neg_mel)  +     
  geom_hline(yintercept = seq(10.5,40.5,10)) + 
  geom_vline(xintercept = c(1.5,2.5,3.5, 4.5))       +                              
  scale_x_discrete(expand=c(0,0)) +                                   
  labs(y="", x="") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), axis.line = element_blank(), legend.title= element_text(size=11), legend.text = element_text(size=10), axis.text = element_text(size=11), legend.text.align = 1, axis.text.y = element_text(face= ifelse(is.element(pca_loadings_neg_mel, c(meta_intersect$emtI_vitro$melanoma.AXL.program,meta_intersect$skinpig_vitro$melanoma.MITF.program)), "bold","plain")))            

dev.off()
                                                                                              