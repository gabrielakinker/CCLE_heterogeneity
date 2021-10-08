# ----------------------------------------------------------------------------------------------------
# Module 3. Comparing RHPs found in cell lines and in human tumor samples.
# ----------------------------------------------------------------------------------------------------

# **************************************************************************
# Basic setup

# load necessary R packages and functions
library(reshape2)
library(ggplot2)
library(scales)
source("custom_magma.R")
source("control_geneset.R")


# read scRNA-seq data (CPM) from cell lines 
expr_ccle <- readRDS("CCLE_heterogeneity_Rfiles/CCLE_scRNAseq_CPM.RDS") # CCLE cell lines
common_genes_ccle  <- Reduce(intersect, lapply(expr_ccle, rownames))
expr_ccle <-  lapply(expr_ccle, function(x) x[common_genes_ccle,])              
                     
# calculate average gene expression 
ave_expr_ccle <- rowMeans(do.call(cbind, expr_ccle))
gene_uni_ccle <- names(sort(ave_expr_ccle, decreasing = T)[1:7000])

# process data                     
expr_ccle <- lapply(expr_ccle, function(x) log2((x/10) + 1))
ave_logexpr_ccle <- rowMeans(do.call(cbind, expr_ccle))          
expr_ccle_cen <- lapply(expr_ccle, function(x) x-ave_logexpr_ccle)                        
  
# read nmf programs from cell lines 
nmf_programs_sig_ccle <- readRDS("Expected_results/module2/nmf_programs_sig_nc_ccle.RDS")   
nmf_meta_programs_ccle <- readRDS("Expected_results/module2/nmf_metaprograms_programs_nc_ccle.RDS")   
                        
# read  metaprograms from cell lines and tumors    
nmf_meta_sig_ccle <- readRDS("Expected_results/module2/nmf_metaprograms_sigtop25_nc_ccle.RDS") # nmf - cell lines
                     
nmf_meta_sig_tumor <-  readRDS("Expected_results/module2/nmf_metaprograms_sigtop25_nc_tumor.RDS") # nmf- tumors 
meta_sig_tumor_lit  <- unlist(list(read.table("CCLE_heterogeneity_Rfiles/metaprograms_tumors_literature.txt", sep = "\t", header = T,stringsAsFactors = F)), recursive=F)  # literature - tumors                     
meta_sig_tumor_all <- c(meta_sig_tumor_lit, nmf_meta_sig_tumor)
meta_sig_tumor_all <- lapply(meta_sig_tumor_all, function(x) x[x!=""])
                             

# ************************************************************************** 
# Compare metaprograms signatures from cell lines and tumors
# (jaccard index and hypergeometric test)                             
                             
# hypergeometric test and jaccard index
nmf_meta_sig_ccle_filt <- lapply(nmf_meta_sig_ccle, function(x){ x[is.element(x, gene_uni_ccle)]})   
meta_sig_tumor_all_filt <- lapply(meta_sig_tumor_all, function(x){ x[is.element(x, gene_uni_ccle)]})                             
                          
meta_jaccard_vivovitro <- sapply(nmf_meta_sig_ccle_filt, function(x) sapply(meta_sig_tumor_all, function(y) length(intersect(x,y))/length(union(x,y) )))             
meta_phyper_vivovitro  <- sapply(nmf_meta_sig_ccle_filt, function(x) sapply(meta_sig_tumor_all, function(y) phyper(q=length(intersect(x,y)), m= length(x), 7000-length(x), k=length(y), lower.tail = F)))   
                                                             
meta_phyper_vivovitro  <- apply(meta_phyper_vivovitro, 2, function(x) p.adjust(x, n=length(meta_phyper_vivovitro), "fdr"))                                     
                      
# remove tumor metaprograms with low similarity to cell line metaprograms
meta_jaccard_vivovitro <- meta_jaccard_vivovitro[apply(meta_jaccard_vivovitro, 1, function(x) length(which(x > 0.05)))!=0,]                      
meta_jaccard_vivovitro <- meta_jaccard_vivovitro[rev(c( 'melanoma.MITF.program', 'HNSCC.Epidif.1', 'episen', 'melanoma.AXL.program', 'HNSCC.PEMT', 'emt', 'GBM.MES1', 'ovarian.ifn', 'HNSCC.Stress', 'stress','GBM.MES2')),   c('skinpig','episen','p53snc', 'emtI','emtIII' , 'emtII',  'ifn', 'stress','protmat', 'protdreg')]

meta_phyper_vivovitro <- meta_phyper_vivovitro[rownames(meta_jaccard_vivovitro), colnames(meta_jaccard_vivovitro)]
 
# plot heatmap                                     
meta_jaccard_vivovitro_melt <- reshape2::melt(meta_jaccard_vivovitro)                      
meta_phyper_vivovitro_melt  <-  reshape2::melt(meta_phyper_vivovitro)                                                                                                                           
red_palette <- brewer.pal(9, "Reds")   
                                     
pdf("Output/module3/meta_vitroVSvivo_jaccard.pdf", width = 6, height = 5)     
ggplot(meta_jaccard_vivovitro_melt, aes(x=Var2, y=Var1, fill=value, color=value)) +
    geom_tile() +
    labs(x="In vitro", y="In vivo") +                                 
     scale_fill_gradient2(limits=c(0.01, 0.16), midpoint = 0.085, low= c( "white",red_palette[1:3]), mid= red_palette[4:6], high = red_palette[7:9] ,   oob=squish, name="Jaccard\nindex")+
  scale_color_gradient2(limits=c(0.01, 0.16), midpoint = 0.085, low= c("white",red_palette[1:3]), mid= red_palette[4:6], high = red_palette[7:9] ,   oob=squish, name="Jaccard\nindex") +
    theme(axis.text.x=element_text(angle=45, hjust=1), panel.background = element_blank(), panel.border = element_rect(fill=F), axis.text=element_text(size=11), axis.title=element_text(size=12))+
 scale_x_discrete(expand=c(0,0), labels = c("Skin pig.", "Epi. sen.", "p53-dep. sen.",  "EMT I", "EMT III" , "EMT II",  "IFN resp.", "Stress", "Prot. mat.", "Prot. deg."   )) +
 scale_y_discrete(expand=c(0,0), labels = rev(c("Melanoma MITF", "HSNCC epi. dif. 1", "Meta epi. sen.", "Melanoma AXL", "HNSCC pEMT", "Meta pEMT", "GBS mes. 1", "Ovarian IFN", "HNSCC pEMT", "Meta stress", "GBM mes. 2")))                                                                          
dev.off()    
    
pdf("Output/module3/meta_vitroVSvivo_phyper.pdf", width = 6, height = 5)                                     
ggplot(meta_phyper_vivovitro_melt, aes(x=Var2, y=Var1, fill=-log10(value), -log10(color=value))) +
    geom_tile() +
    labs(x="In vitro", y="In vivo") +                                 
     scale_fill_gradient2(limits=c(1, 12), midpoint = 6.5, low= c( "white",red_palette[1:3]), mid= red_palette[4:6], high = red_palette[7:9] ,   oob=squish, name="P adjust\n(fdr)")+
  scale_color_gradient2(limits=c(1, 12), midpoint = 6.5, low= c("white",red_palette[1:3]), mid= red_palette[4:6], high = red_palette[7:9] ,   oob=squish, name="P adjust\n(fdr)") +
    theme(axis.text.x=element_text(angle=45, hjust=1), panel.background = element_blank(), panel.border = element_rect(fill=F), axis.text=element_text(size=11), axis.title=element_text(size=12))+                           
 scale_x_discrete(expand=c(0,0), labels = c("Skin pig.", "Epi. sen.", "p53-dep. sen.",  "EMT I", "EMT III" , "EMT II",  "IFN resp.", "Stress", "Prot. mat.", "Prot. deg."   )) +
 scale_y_discrete(expand=c(0,0), labels = rev(c("Melanoma MITF", "HSNCC epi. dif. 1", "Meta epi. sen.", "Melanoma AXL", "HNSCC pEMT", "Meta pEMT", "GBS mes. 1", "Ovarian IFN", "HNSCC stress", "Meta stress", "GBM mes. 2")))                                    
dev.off()    
                                                       
# ************************************************************************** 
# Compare individual programs (not metaprograms) from cell lines to tumor metaprograms 
# (mean jaccard index and mean correlation of program scores)                                          
  
# calculate the jaccard index between nmf programs from cell lines and metaprograms from tumors 
indprog_jaccard_vivovitro <- sapply(meta_sig_tumor_all, function(x) apply(nmf_programs_sig_ccle[,unlist(nmf_meta_programs_ccle)], 2, function(y) length(intersect(x,y))/length(union(x,y))))

# calculate program scores - cell line programs
nmf_meta_sig_ccle_filtII <- apply(nmf_programs_sig_ccle[,unlist(nmf_meta_programs_ccle)], 2, function(x) x[is.element(x, common_genes_ccle)])
nmf_meta_sig_ccle_filtII_ctr <- lapply(nmf_meta_sig_ccle_filtII, function(x) control_geneset(ave_tpm = ave_expr_ccle, program = x, bins = 42))

vitro_scores_ccle <- list()

for(i in names(nmf_meta_sig_ccle_filtII)) {
  a <- gsub(".{4}$", "", i)
  b <- expr_ccle[[a]]
  vitro_scores_ccle[[i]] <- colMeans(b[ nmf_meta_sig_ccle_filtII[[i]],]) - colMeans(b[nmf_meta_sig_ccle_filtII_ctr[[i]],])
}

# calculate program scores - tumor metaprograms
meta_sig_tumor_all_filtII <- lapply(meta_sig_tumor_all, function(x) x[is.element(x, common_genes_ccle)])
meta_sig_tumor_all_filtII_ctr <- lapply(meta_sig_tumor_all_filtII, function(x) control_geneset(ave_tpm = ave_expr_ccle, program = x, bins = 42))

vivo_scores_ccle <- list()
for(i in names(nmf_meta_sig_ccle_filtII)) {
  a <- gsub(".{4}$", "", i)
  b <- expr_ccle[[a]]
  vivo_scores_ccle[[i]] <- sapply(meta_sig_tumor_all_filtII, function(x) colMeans(b[x,])) - sapply(meta_sig_tumor_all_filtII_ctr, function(x) colMeans(b[x,]))
}

# calculate correlation between program scores - vitro vs vivo
indprog_corr_vivovitro <- sapply(unlist(nmf_meta_programs_ccle), function(x) cor(vitro_scores_ccle[[x]], vivo_scores_ccle[[x]]))
rownames(indprog_corr_vivovitro) <- names(meta_sig_tumor_all)                       
colnames(indprog_corr_vivovitro) <- unlist(nmf_meta_programs_ccle)
                               
# calculating random correlation and jaccard index
jaccard_perm <- list()
corr_perm <- list()

for(i in 1:100) {
  a <- lapply(nmf_meta_sig_ccle_filtII, function(x) control_geneset(ave_tpm = ave_expr_ccle, program = x, bins = 42, size=1, seed = i))
  a_crt <- lapply(a, function(x) control_geneset(ave_tpm = ave_expr_ccle, program = x, bins = 42))     
  b <- list()
  for(j in names(a)) {
    c <- gsub(".{4}$", "", j)
    d <- expr_ccle[[c]]
    b[[j]] <- colMeans(d[a[[j]],]) - colMeans(d[a_crt[[j]],])
  }
  jaccard_perm[[i]] <- sapply(meta_sig_tumor_all, function(x) sapply(a, function(y) length(intersect(x,y))/length(union(x,y))))
  corr_perm[[i]] <- sapply(names(a), function(x) cor(b[[x]], vivo_scores_ccle[[x]]))                                   
}


# aggregate results by metaprogram
indprog_corr_vivovitro <- data.frame(aggregate(t(indprog_corr_vivovitro), list(melt(nmf_meta_programs_ccle)$L1), mean), row.names=1) 
indprog_jaccard_vivovitro <- data.frame(aggregate(indprog_jaccard_vivovitro, list(melt(nmf_meta_programs_ccle)$L1), mean), row.names=1) 

# plot                                                
jaccard_corr_plot <- data.frame(melt(as.matrix(indprog_corr_vivovitro)),  melt(as.matrix(indprog_jaccard_vivovitro)))                 
                                                        
pdf("Output/module3/indprog_vitroVSvivo_jaccard&corr.pdf", width = 9, height = 5, onefile = F)                                                                                                
ggplot(jaccard_corr_plot, aes(x=value, y=value.1)) +
  geom_point(size=4, shape=21, fill="gray90") +
  facet_wrap(facets = vars(Var1), nrow = 2) +
  geom_hline(yintercept =  quantile(unlist(jaccard_perm), 0.999), linetype="dashed") +
  geom_vline(xintercept =  quantile(unlist(corr_perm), 0.999), linetype="dashed") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), strip.text = element_text(size=13), strip.background = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=13)) +
  labs(y="Mean similarity\n(Jaccard Index)", x="Mean correlation") +
  scale_x_continuous(breaks=seq(-0.5, 0.5, 0.5))
dev.off()
                                                            
                                     
