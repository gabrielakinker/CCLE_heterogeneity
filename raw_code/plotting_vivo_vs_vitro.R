library(RColorBrewer)
library(reshape2)
library(scales)
library(ggplot2)
library(ggrepel)
source("~/CCLE/nmf/function_get_control_geneset.R")
        
                 
#########################################################################################################################
################################ plotting in vivo vs in vitro programs - heatmap p hyper test per program ###################################################

### reading and processing TPM data                 
expr <- readRDS("~/CCLE/data processing/scCCLE&pool1_tpm_by_cell_line.rds") 
expr <- lapply(expr, function(x) log2((x/10) + 1))
expr <- lapply(expr, function(x) x[apply(x, 1, function(y) length(which(y > 3.5)) > ncol(x)*0.02),])
               

### reading ccle nmf programs (already ordered by hc)
nmf_ccle <- readRDS("~/CCLE/nmf/nc_programs_hc_ordered.RDS")

### reading metaprograms
# literature
meta_vivo_literature <- read.table("~/CCLE/datasets/genesets/in_vivo_signatures.txt", sep = "\t", header = T,stringsAsFactors = F) 
            
# nmf
meta_vivo_nmf <-  read.table("~/CCLE/nmf/in vivo/tumors/in_vivo_signatures_nmf.txt", sep = "\t", header = T,stringsAsFactors = F)

# combine
vivo_signatures <- unlist(unname(list(meta_vivo_literature,meta_vivo_nmf)), recursive = F )
vivo_signatures <- lapply(vivo_signatures, function(x) x[x!=""])

### calculating p hyper test - similarity beween in vitro vs. in vivo programs
p_hyper <- data.frame(matrix(ncol=length(vivo_signatures), nrow=ncol(nmf_ccle)), row.names=colnames(nmf_ccle))            
colnames(p_hyper) <- names(vivo_signatures)
                          
for(i in colnames(nmf_ccle)) {
    a <- nmf_ccle[,i]
    univ <- rownames(expr[[gsub(".{4}$", "", i)]])
    for(j in names(vivo_signatures)) { 
    b <- vivo_signatures[[j]][is.element(vivo_signatures[[j]], univ)]
    p_hyper[i,j] <- phyper(length(intersect(a,b)), length(a), length(univ) - length(a), length(b), lower.tail = F)
    }
}  
                          
                          
p_hyper <- apply(p_hyper, 2, function(x) p.adjust(x, "fdr"))                          
p_hyper <- -log10(p_hyper)       
p_hyper <- p_hyper[,order.dendrogram(reorder(as.dendrogram(hclust(as.dist(1-cor(p_hyper)))), colMeans(p_hyper))) ]                                 
                 
                                                                         
red_palette <- brewer.pal(9, "Greys")                                                                          
 
pdf("nc_vivo_vs_vitro_phyper.pdf", width = 5.5, height = 7.1)                           
ggplot(data = melt(p_hyper), aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(1, 12), midpoint = 6.5, low= c( "white",red_palette[1:3]), mid= red_palette[4:6], high = red_palette[7:9] ,   oob=squish, name="P adjust\n(fdr)")+
  scale_color_gradient2(limits=c(1, 12), midpoint = 6.5, low= c("white",red_palette[1:3]), mid= red_palette[4:6], high = red_palette[7:9] ,   oob=squish, name="P adjust\n(fdr)") +
  labs(y="", x="") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), axis.line = element_blank(), axis.ticks = element_blank(), legend.text = element_text(size=12),  axis.text = element_blank(), legend.text.align = 0.5)    +
  geom_hline(yintercept = c(25,53, 75, 95, 109, 132, 167, 250, 278), size=0.3)         
dev.off()                     
                          
                             
#########################################################################################################################
################################ plotting in vivo vs in vitro programs - heatmap p hyper test per metaprogram (RHP) ###################################################
 
                 
### reading and processing TPM data                 
expr <- readRDS("~/CCLE/data processing/scCCLE&pool1_tpm_by_cell_line.rds")                 
common_genes <- Reduce(intersect, lapply(expr, rownames))
expr <- lapply(expr, function(x) x[common_genes,])
               
### calculating average TPM values
ave_expr <- rowMeans(do.call(cbind, expr))
                 
### reading  programs       
mel_mitf <- scan("~/CCLE/nmf/metaprograms/metaprograms - top genes/meta_mel_mitf.txt", "character")   
mel_emt <- scan("~/CCLE/nmf/metaprograms/metaprograms - top genes/meta_mel_emt.txt", "character")     
emt <- scan("~/CCLE/nmf/metaprograms/metaprograms - top genes/meta_emt.txt", "character")
emt_snc <- scan("~/CCLE/nmf/metaprograms/metaprograms - top genes/meta_emt_snc.txt", "character")   
emt <- scan("~/CCLE/nmf/metaprograms/metaprograms - top genes/meta_emt.txt", "character")
emt_snc <- scan("~/CCLE/nmf/metaprograms/metaprograms - top genes/meta_emt_snc.txt", "character")
ifn <- scan("~/CCLE/nmf/metaprograms/metaprograms - top genes/meta_ifn.txt", "character")
snc_classic <- scan("~/CCLE/nmf/metaprograms/metaprograms - top genes/meta_snc_classic.txt", "character")
snc_epi <- scan("~/CCLE/nmf/metaprograms/metaprograms - top genes/meta_snc_epi.txt", "character")
stress <- scan("~/CCLE/nmf/metaprograms/metaprograms - top genes/meta_stress.txt", "character")
endo <- scan("~/CCLE/nmf/metaprograms/metaprograms - top genes/meta_endo.txt", "character")
prot <- scan("~/CCLE/nmf/metaprograms/metaprograms - top genes/meta_prot.txt", "character")

### list of metaprograms
metaprograms_vitro <- list("mel_mitf"= mel_mitf,"mel_emt"= mel_emt,"emt"= emt, "emt_snc"= emt_snc, "ifn"=ifn, "snc_classic"=snc_classic, "snc_epi"=snc_epi, "stress"=stress, "endo"=endo, "prot"=prot)
                 
                 
### reading metaprograms
# literature
meta_vivo_literature <- read.table("~/CCLE/datasets/genesets/in_vivo_signatures.txt", sep = "\t", header = T,stringsAsFactors = F) 
            
# nmf
meta_vivo_nmf <-  read.table("~/CCLE/nmf/in vivo/tumors/in_vivo_signatures_nmf.txt", sep = "\t", header = T,stringsAsFactors = F)

# combine
metaprograms_vivo <- unlist(unname(list(meta_vivo_literature,meta_vivo_nmf)), recursive = F )
metaprograms_vivo <- lapply(metaprograms_vivo, function(x) x[x!=""])
                 
### p hyper test
gene_uni <- names(sort(ave_expr, decreasing = T)[1:7000])
 
metaprograms_vitro <- lapply(metaprograms_vitro, function(x){ x[is.element(x, gene_uni)]})   
metaprograms_vivo <- lapply(metaprograms_vivo, function(x){ x[is.element(x, gene_uni)]})                             
                          
jaccard_index <- sapply(metaprograms_vitro, function(x) sapply(metaprograms_vivo, function(y) length(intersect(x,y))/length(union(x,y) )))                          

phyper_test <- sapply(metaprograms_vitro, function(x) sapply(metaprograms_vivo, function(y) phyper(q=length(intersect(x,y)), m= length(x), 7000-length(x), k=length(y), lower.tail = F)))   
                                                             
phyper_test  <- apply(phyper_test, 2, function(x) p.adjust(x, n=length(phyper_test), "fdr"))                                     
                      
### plotting    
jaccard_index <- jaccard_index[apply(jaccard_index, 1, function(x) length(which(x > 0.05)))!=0,]                      
jaccard_index <- jaccard_index[rev(c( 'melanoma.MITF.program', 'HNSCC.Epidif.1', 'snc_epi', 'melanoma.AXL.program', 'HNSCC.PEMT', 'pemt', 'GBM.MES1', 'ovarian.ifn', 'HNSCC.Stress', 'stress','GBM.MES2')),   c('mel_mitf','snc_epi','snc_classic', 'mel_emt','emt_snc' , 'emt',  'ifn', 'stress','endo', 'prot')]

phyper_test <- phyper_test[rownames(jaccard_index), colnames(jaccard_index)]
                      
jaccard_index_melt <- reshape2::melt(jaccard_index)                      
phyper_test_melt  <-  reshape2::melt(phyper_test)                                                              
                                                             
red_palette <- brewer.pal(9, "Reds")   
 
pdf("RHP_vivo_phyper.pdf", width = 6, height = 5)                                     
ggplot(phyper_test_melt, aes(x=Var2, y=Var1, fill=-log10(value), -log10(color=value))) +
    geom_tile() +
    labs(x="In vitro", y="In vivo") +                                 
     scale_fill_gradient2(limits=c(1, 12), midpoint = 6.5, low= c( "white",red_palette[1:3]), mid= red_palette[4:6], high = red_palette[7:9] ,   oob=squish, name="P adjust\n(fdr)")+
  scale_color_gradient2(limits=c(1, 12), midpoint = 6.5, low= c("white",red_palette[1:3]), mid= red_palette[4:6], high = red_palette[7:9] ,   oob=squish, name="P adjust\n(fdr)") +
    theme(axis.text.x=element_text(angle=45, hjust=1), panel.background = element_blank(), panel.border = element_rect(fill=F), axis.text=element_text(size=11), axis.title=element_text(size=12))+                           
 scale_x_discrete(expand=c(0,0), labels = c("Skin pig.", "Epi. sen.", "p53-dep. sen.",  "EMT I", "EMT III" , "EMT II",  "IFN resp.", "Stress", "Prot.mat.", "Prot. deg."   )) +
 scale_y_discrete(expand=c(0,0), labels = rev(c("Melanoma MITF", "HSNCC epi. dif. 1", "Meta epi. sen.", "Melanoma AXL", "HNSCC pEMT", "Meta pEMT", "GBS mes. 1", "Ovarian IFN", "HNSCC stress", "Meta stress", "GBM mes. 2")))                                    
dev.off()    
                                     
pdf("RHP_vivo_jaccard.pdf", width = 6, height = 5)     
ggplot(jaccard_index_melt, aes(x=Var2, y=Var1, fill=value, color=value)) +
    geom_tile() +
    labs(x="In vitro", y="In vivo") +                                 
     scale_fill_gradient2(limits=c(0.01, 0.16), midpoint = 0.085, low= c( "white",red_palette[1:3]), mid= red_palette[4:6], high = red_palette[7:9] ,   oob=squish, name="Jaccard\nindex")+
  scale_color_gradient2(limits=c(0.01, 0.16), midpoint = 0.085, low= c("white",red_palette[1:3]), mid= red_palette[4:6], high = red_palette[7:9] ,   oob=squish, name="Jaccard\nindex") +
    theme(axis.text.x=element_text(angle=45, hjust=1), panel.background = element_blank(), panel.border = element_rect(fill=F), axis.text=element_text(size=11), axis.title=element_text(size=12))+
 scale_x_discrete(expand=c(0,0), labels = c("Skin pig.", "Epi. sen.", "p53-dep. sen.",  "EMT I", "EMT III" , "EMT II",  "IFN resp.", "Stress", "Prot.mat.", "Prot. deg."   )) +
 scale_y_discrete(expand=c(0,0), labels = rev(c("Melanoma MITF", "HSNCC epi. dif. 1", "Meta epi. sen.", "Melanoma AXL", "HNSCC pEMT", "Meta pEMT", "GBS mes. 1", "Ovarian IFN", "HNSCC pEMT", "Meta stress", "GBM mes. 2")))                                                                          
dev.off()    
                 
#########################################################################################################################
################################ plotting in vivo vs in vitro programs -  mean correlation/similarity ###################################################
                
                 
### reading TPM data                     
expr <- readRDS("~/CCLE/data processing/scCCLE&pool1_tpm_by_cell_line.rds") 
common_genes <- Reduce(intersect, lapply(expr, rownames))
expr <- lapply(expr, function(x) x[common_genes,])

### calculating average TPM values
ave_tpm <- rowMeans(sapply(expr, rowMeans))

### processing TPM data
expr <- lapply(expr, function(x) log2((x/10) + 1))
ave_log_tpm <- rowMeans(do.call(cbind, expr))
expr <- lapply(expr, function(x) x-ave_log_tpm) 

### reading ccle nmf programs (already ordered by hc)
nmf_ccle <- readRDS("~/CCLE/nmf/nc_programs_hc_ordered.RDS")

### reading metaprograms
# literature
meta_vivo_literature <- read.table("~/CCLE/datasets/genesets/in_vivo_signatures.txt", sep = "\t", header = T,stringsAsFactors = F) 
            
# nmf
meta_vivo_nmf <-  read.table("~/CCLE/nmf/in vivo/tumors/in_vivo_signatures_nmf.txt", sep = "\t", header = T,stringsAsFactors = F)

# combine
vivo_signatures <- unlist(unname(list(meta_vivo_literature,meta_vivo_nmf)), recursive = F )
vivo_signatures <- lapply(vivo_signatures, function(x) x[x!=""])

### calculating the similarity between in vitro and in vivo
intersect_vivo <- sapply(vivo_signatures, function(x) apply(nmf_ccle, 2, function(y) length(intersect(x,y))/length(union(x,y))))


### calculating program scores - in vitro
nmf_ccle_filt <- apply(nmf_ccle, 2, function(x) x[is.element(x, common_genes)])
nmf_ccle_filt_ctr <- lapply(nmf_ccle_filt, function(x) get_control_geneset(ave_tpm = ave_tpm, program = x, bins = 42))

nmf_ccle_scores <- list()

for(i in names(nmf_ccle_filt)) {
  a <- gsub(".{4}$", "", i)
  b <- expr[[a]]
  nmf_ccle_scores[[i]] <- colMeans(b[ nmf_ccle_filt[[i]],]) - colMeans(b[nmf_ccle_filt_ctr[[i]],])
}

### calculating program scores - in vivo
vivo_signatures_filt <- lapply(vivo_signatures, function(x) x[is.element(x, common_genes)])
vivo_signatures_filt_crt <- lapply(vivo_signatures_filt, function(x) get_control_geneset(ave_tpm = ave_tpm, program = x, bins = 42))

in_vivo_scores <- list()
for(i in names(nmf_ccle_filt)) {
  a <- gsub(".{4}$", "", i)
  b <- expr[[a]]
  in_vivo_scores[[i]] <- sapply(vivo_signatures_filt, function(x) colMeans(b[x,])) - sapply(vivo_signatures_filt_crt, function(x) colMeans(b[x,]))
}

### calculating correlation vitro vs vivo
corr_vitro_vivo <- sapply(colnames(nmf_ccle), function(x) cor(nmf_ccle_scores[[x]], in_vivo_scores[[x]]))
rownames(corr_vitro_vivo) <- names(vivo_signatures)                       

### calculating random maximum correlation and maximum similarity
sim_perm <- list()
cor_perm <- list()

for(i in 1:100) {
  a <- future.apply::future_lapply(nmf_ccle_filt, function(x) get_control_geneset(ave_tpm = ave_tpm, program = x, bins = 42, size=1, seed = i))
  a_crt <- future.apply::future_lapply(a, function(x) get_control_geneset(ave_tpm = ave_tpm, program = x, bins = 42))     
  b <- list()
  for(j in names(a)) {
    c <- gsub(".{4}$", "", j)
    d <- expr[[c]]
    b[[j]] <- colMeans(d[a[[j]],]) - colMeans(d[a_crt[[j]],])
  }
  sim_perm[[i]] <- future.apply::future_sapply(vivo_signatures, function(x) sapply(a, function(y) length(intersect(x,y))/length(union(x,y))))
  cor_perm[[i]] <- future.apply::future_sapply(names(a), function(x) cor(b[[x]], in_vivo_scores[[x]]))                                   
}



### plotting mean correlations/similarities - by in vitro program

### aggregate corr and sim per metaprogrma
metaprograms <- rep(NA, ncol(nmf_ccle))  
metaprograms[10:25] <- "Skin\npigmenation"
metaprograms[46:53] <- "EMT I"
metaprograms[ 57:75] <- "EMT II"
metaprograms[ 81:95] <- "IFN\nresponse"
metaprograms[ 98:109] <- "EMT III"
metaprograms[ 116:132] <- "p53-dep.\nsenescence"
metaprograms[ 140:167] <- "Epi.\nsenescence"
metaprograms[197:250] <- "Stress"
metaprograms[261:278] <- "Protein\nmaturation"
metaprograms[279:432] <- "Proteasomal\ndegradation"                                                
   
cor_aggre <- data.frame(aggregate(t(corr_vitro_vivo), list(metaprograms), mean), row.names=1) 
sim_aggre <- data.frame(aggregate(intersect_vivo, list(metaprograms), mean), row.names=1) 

### plot                                                
sim_cor_plot <- data.frame(melt(as.matrix(cor_aggre)),  melt(as.matrix(sim_aggre)))                 
                                               
          
pdf("nc_vitro_vivo_sim_cor.pdf", width = 9, height = 5, onefile = F)                                                                                                 
ggplot(sim_cor_plot, aes(x=value, y=value.1)) +
  geom_point(size=4, shape=21, fill="gray90") +
  facet_wrap(facets = vars(Var1), nrow = 2) +
  geom_hline(yintercept =  quantile(unlist(sim_perm), 0.999), linetype="dashed") +
  geom_vline(xintercept =  quantile(unlist(cor_perm), 0.999), linetype="dashed") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), strip.text = element_text(size=13), strip.background = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=13)) +
  labs(x="Mean similarity\n(Jaccard Index)", y="Mean correlation") +
  scale_x_continuous(breaks=seq(-0.5, 0.5, 0.5))
dev.off()
                                                
                                                
    
#########################################################################################################################
################################ plotting in vivo vs in vitro programs - aggregating cell lines and tumors ###################################################
                                               
### reading TPM expression data cell lines
expr_ccle <- readRDS("~/CCLE/data processing/scCCLE&pool1_tpm_by_cell_line.rds")                       
rowmeans_ccle <- lapply(expr_ccle, rowMeans)                                                
expr_ccle<- lapply(expr_ccle, function(x) {log2((x/10) + 1)})                                           
expr_ccle <- lapply(expr_ccle, function(x) {x-rowMeans(x)})
                                                                                                                                                                                         
### reading log TPM expression data - tumors
expr_tumors <- readRDS("~/CCLE/datasets/scRNA-seq/scRNAseq_tumors_logtpm.RDS")
expr_tumors <- unlist(expr_tumors, recursive = F)
expr_tumors <- expr_tumors[sapply(expr_tumors, function(x) nrow(x) >= 50)]
expr_tumors <- lapply(expr_tumors, t)
rowmeans_tumors <- lapply(expr_tumors, function(x) rowMeans(10*((2^x)-1)))  
expr_tumors <- lapply(expr_tumors, function(x) {x-rowMeans(x)})


### getting cell lines and tumors selected programs                                

# episen and emt                                  
episen_emt_programs_ccle <- intersect(gsub(".{4}$", "",scan("~/CCLE/nmf/metaprograms/meta_snc_epi_programs.txt", "character")), gsub(".{4}$", "",scan("~/CCLE/nmf/metaprograms/meta_emt_programs.txt", "character")))
episen_emt_programs_ccle <- unique(episen_emt_programs_ccle[grep("UPPER", episen_emt_programs_ccle)])                                              
                                  
episen_emt_programs_tumor <- unique(intersect(gsub(".{4}$", "",scan("~/CCLE/nmf/in vivo/tumors/meta_epi_snc_vivo_programs.txt", "character")), gsub(".{4}$", "",scan("~/CCLE/nmf/in vivo/tumors/meta_emt_vivo_programs.txt", "character"))))

# mel and emt      
mitf_emt_programs_ccle <- intersect(gsub(".{4}$", "",scan("~/CCLE/nmf/metaprograms/meta_mel_mitf_programs.txt", "character")), gsub(".{4}$", "",scan("~/CCLE/nmf/metaprograms/meta_mel_emt_programs.txt", "character")))
                                  
mitf_emt_programs_tumor <- names(expr_tumors)[grep("mel", names(expr_tumors))]                                  
### reading programs                                                                   

# vitro
emt <- read.table("~/CCLE/nmf/metaprograms/metatable_emt.txt", sep="\t", header = T, stringsAsFactors = F)[1:100,1]
snc_epi <- read.table("~/CCLE/nmf/metaprograms/metatable_snc_epi.txt", sep="\t", header = T, stringsAsFactors = F)[1:100,1]      
mel_emt <- read.table("~/CCLE/nmf/metaprograms/metatable_mel_emt.txt", sep="\t", header = T, stringsAsFactors = F)[1:100,1]
mel_mitf <- read.table("~/CCLE/nmf/metaprograms/metatable_mel_mitf.txt", sep="\t", header = T, stringsAsFactors = F)[1:100,1]                                        
                                  
vitro_signatures <- c(list("emt_vitro"= emt, "snc_epi_vitro"=snc_epi, "mel_emt_vitro"= mel_emt, "mel_mitf_vitro"=mel_mitf))                             

# vivo                                  
meta_vivo_literature <- read.table("~/CCLE/datasets/genesets/in_vivo_signatures.txt", sep = "\t", header = T,stringsAsFactors = F)  # literature
            
meta_vivo_nmf <-  read.table("~/CCLE/nmf/in vivo/tumors/in_vivo_signatures_nmf.txt", sep = "\t", header = T,stringsAsFactors = F) # nmf

vivo_signatures <- unlist(unname(list(meta_vivo_literature,meta_vivo_nmf)), recursive = F )
vivo_signatures <- lapply(vivo_signatures, function(x) x[x!=""])                                  
  
# intersection betweeb vitro and vivo                          
intersect <- lapply(vitro_signatures, function(x) lapply(vivo_signatures, function(y) intersect(x,y)))

### pca
                                                         
# emt and episen     
expr_emt_episen_vitro <- expr_ccle[episen_emt_programs_ccle]
expr_emt_episen_vivo <- expr_tumors[episen_emt_programs_tumor]                                                                                                   
expr_emt_episen_vitro <- sapply(expr_emt_episen_vitro, function(x) x[Reduce(intersect, c(lapply( expr_emt_episen_vitro, rownames), lapply(expr_emt_episen_vivo, rownames))),])
expr_emt_episen_vivo <- sapply(expr_emt_episen_vivo, function(x) x[Reduce(intersect, c(lapply(expr_emt_episen_vitro, rownames), lapply(expr_emt_episen_vivo, rownames))),]) 
                               
cell_number_emt_episen <- min(c(sapply(expr_emt_episen_vitro, ncol),sapply(expr_emt_episen_vivo, ncol)))
                               
expr_emt_episen_vitro <- do.call(cbind,expr_emt_episen_vitro )        
expr_emt_episen_vivo <-  do.call(cbind, expr_emt_episen_vivo)                  
                                                         
                               
expr_emt_episen_vitro_filt <- expr_emt_episen_vitro[order(rowMeans(cbind(rowMeans(sapply(rowmeans_ccle[episen_emt_programs_ccle], function(x) x[rownames(expr_emt_episen_vitro)])), rowMeans(sapply(rowmeans_tumors[episen_emt_programs_tumor], function(x) x[rownames(expr_emt_episen_vitro)])))), decreasing=T)[1:4500] ,]
expr_emt_episen_vivo_filt <- expr_emt_episen_vivo[order(rowMeans(cbind(rowMeans(sapply(rowmeans_ccle[episen_emt_programs_ccle], function(x) x[rownames(expr_emt_episen_vivo)])), rowMeans(sapply(rowmeans_tumors[episen_emt_programs_tumor], function(x) x[rownames(expr_emt_episen_vivo)])))), decreasing=T)[1:4500] ,]
                                                                                          
expr_emt_episen_final <- cbind(expr_emt_episen_vitro_filt,expr_emt_episen_vivo_filt)    
expr_emt_episen_pca <- prcomp(t(expr_emt_episen_final))                       
                               
                                                                                                                                               
# emt and mel    

expr_mitf_emt_vitro <- expr_ccle[mitf_emt_programs_ccle]
expr_mitf_emt_vivo <- expr_tumors[mitf_emt_programs_tumor]                                                                                                   
expr_mitf_emt_vitro <- sapply(expr_mitf_emt_vitro, function(x) x[Reduce(intersect, c(lapply( expr_mitf_emt_vitro, rownames), lapply(expr_mitf_emt_vivo, rownames))),])
expr_mitf_emt_vivo <- sapply(expr_mitf_emt_vivo, function(x) x[Reduce(intersect, c(lapply(
expr_mitf_emt_vitro, rownames), lapply(expr_mitf_emt_vivo, rownames))),])    
                             
expr_mitf_emt_vitro <- do.call(cbind, expr_mitf_emt_vitro)        
expr_mitf_emt_vivo <- do.call(cbind, expr_mitf_emt_vivo)                               
                                                         
expr_mitf_emt_vitro_filt <- expr_mitf_emt_vitro[order(rowMeans(cbind(rowMeans(sapply(rowmeans_ccle[mitf_emt_programs_ccle], function(x) x[rownames(expr_mitf_emt_vitro)])), rowMeans(sapply(rowmeans_tumors[mitf_emt_programs_tumor], function(x) x[rownames(expr_mitf_emt_vitro)])))), decreasing=T)[1:4500] ,]
expr_mitf_emt_vivo_filt <- expr_mitf_emt_vivo[order(rowMeans(cbind(rowMeans(sapply(rowmeans_ccle[mitf_emt_programs_ccle], function(x) x[rownames(expr_mitf_emt_vivo)])), rowMeans(sapply(rowmeans_tumors[mitf_emt_programs_tumor], function(x) x[rownames(expr_mitf_emt_vivo)])))), decreasing=T)[1:4500] ,]
        
expr_mitf_emt_final <- cbind(expr_mitf_emt_vitro_filt, expr_mitf_emt_vivo_filt)                      
expr_mitf_emt_pca <- prcomp(t(expr_mitf_emt_final))                       
                                                                                                                                                                              
### plot PCA
                                                                                                   
# episen and emt                                  
expr_emt_episen_plot <- data.frame(expr_emt_episen_pca$x[,1:5], type=rep(c("in vitro", "in vivo"), c(ncol(expr_emt_episen_vitro), ncol(expr_emt_episen_vivo))), emt=colMeans(expr_emt_episen_final[is.element(rownames(expr_emt_episen_final), intersect$emt_vitro$HNSCC.PEMT),]), snc=colMeans(expr_emt_episen_final[is.element(rownames(expr_emt_episen_final), intersect$snc_epi_vitro$`HNSCC.Epidif.1`),]))    

expr_emt_episen_plot$names <- c(rep(episen_emt_programs_ccle, sapply(expr_ccle[episen_emt_programs_ccle],ncol)),rep(episen_emt_programs_tumor, sapply(expr_tumors[episen_emt_programs_tumor],ncol)) )
expr_emt_episen_plot$cell_line <- NA
expr_emt_episen_plot$cell_line[expr_emt_episen_plot$type=="in vitro"] <- expr_emt_episen_plot$names[expr_emt_episen_plot$type=="in vitro"] 
expr_emt_episen_plot$cell_line <-   gsub("_UPPER_AERODIGESTIVE_TRACT", "", expr_emt_episen_plot$cell_line)   
expr_emt_episen_plot$tumor <- NA
expr_emt_episen_plot$tumor[expr_emt_episen_plot$type=="in vivo"] <- expr_emt_episen_plot$names[expr_emt_episen_plot$type=="in vivo"]      

                                                                                                                                         
pdf("emt_episen_pca_program.pdf", width = 3.9, height = 3.8)  
ggplot(expr_emt_episen_plot, aes(x=PC2, y=PC4+PC5, color=emt-snc)) +
            geom_point( size=0.5) +
            scale_color_gradient2(limits=c(-3, 3), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(EMT II - Episen)") +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))
dev.off()
                                                                                                                                         
pdf("emt_episen_pca_type.pdf", width = 3.95, height = 3.7)  
ggplot(expr_emt_episen_plot, aes(x=PC2, y= PC4+PC5, color=type)) +
            geom_point( size=0.5) +
             scale_color_manual(values=c("goldenrod", "mediumpurple"), name="") +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
            guides(color=guide_legend(override.aes = list(size=4)))  
dev.off()
                                                                                                                                         
pdf("emt_episen_pca_cellline.pdf", width = 3.9, height = 3.8)  
ggplot(expr_emt_episen_plot, aes(x=PC2, y= PC4 + PC5, color=cell_line)) +
            geom_point(data=expr_emt_episen_plot[is.na(expr_emt_episen_plot$cell_line),], aes(x=PC2, y= PC4+PC5), size=0.5, color="gray90") + 
             geom_point( size=0.5) +
             scale_color_brewer(palette = "Set2", name="", breaks=unique(expr_emt_episen_plot$cell_line)[!is.na(unique(expr_emt_episen_plot$cell_line))]) +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
            guides(color=guide_legend(override.aes = list(size=3), nrow=2) )           
dev.off()
                                                                                                       
pdf("emt_episen_pca_tumor.pdf", width = 3.9, height = 3.8) 
ggplot(expr_emt_episen_plot, aes(x=PC2, y= PC4 + PC5, color=tumor)) +
            geom_point( size=0.5) +
             scale_color_brewer(palette = "Set2", name="", na.value="gray90", breaks=unique(expr_emt_episen_plot$tumor)[!is.na(unique(expr_emt_episen_plot$tumor))], labels=c("HNSCC 18","HNSCC 26", "HNSCC 25", "HNSCC 5", "HNSCC 17" )) +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
            guides(color=guide_legend(override.aes = list(size=3), nrow=2))                                              
dev.off()                                   
                                                                                                                                         
  
pdf("emt_episen_pca_programII.pdf", width = 3.9, height = 3.8)  
ggplot(expr_emt_episen_plot, aes(x=PC2, y=PC4, color=emt-snc)) +
            geom_point( size=0.5) +
            scale_color_gradient2(limits=c(-3, 3), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(EMT II - Episen)") +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))
dev.off()
                                                                                                                                         
pdf("emt_episen_pca_typeII.pdf", width = 3.95, height = 3.7)  
ggplot(expr_emt_episen_plot, aes(x=PC2, y= PC4, color=type)) +
            geom_point( size=0.5) +
            scale_color_manual(values=c("goldenrod", "mediumpurple"), name="") +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
            guides(color=guide_legend(override.aes = list(size=4)))  
dev.off()                                                                                                                                         
                                                                                        pdf("emt_episen_pca_programIII.pdf", width = 3.9, height = 3.8)  
ggplot(expr_emt_episen_plot, aes(x=PC2, y=PC5, color=emt-snc)) +
            geom_point( size=0.5) +
            scale_color_gradient2(limits=c(-3, 3), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(EMT II - Episen)") +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))
dev.off()
                                                                                                                                         
pdf("emt_episen_pca_typeIII.pdf", width = 3.95, height = 3.7)  
ggplot(expr_emt_episen_plot, aes(x=PC2, y= PC5, color=type)) +
            geom_point( size=0.5) +
             scale_color_manual(values=c("goldenrod", "mediumpurple"), name="") +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
            guides(color=guide_legend(override.aes = list(size=4)))  
dev.off()                                                                                                                                               
                                                                                           
pca_loadings_pos_hnscc <- rev(melt(apply(expr_emt_episen_pca$rotation[,1:5], 2, function(x) names(sort(x, decreasing=T))[1:10]))$value)              
                                 
pca_loadings_pos_hnscc_plot <- as.matrix(expr_emt_episen_pca$rotation[as.character(pca_loadings_pos_hnscc),1:5] )
rownames(pca_loadings_pos_hnscc_plot) <- NULL
pca_loadings_pos_hnscc_plot <- melt(as.matrix(pca_loadings_pos_hnscc_plot))

pdf("hnscc_pca_pos.pdf", width = 4.5, height = 7.5)                                     
ggplot(data = pca_loadings_pos_hnscc_plot, aes(x=Var2, y=factor(Var1), fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(-0.12, 0.12), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="PC loadings")+
  scale_color_gradient2(limits=c(-0.12, 0.12), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="PC loadings") +
  scale_y_discrete(labels=pca_loadings_pos_hnscc)  +     
  geom_hline(yintercept = seq(10.5,40.5,10)) + 
  geom_vline(xintercept = c(1.5,2.5,3.5, 4.5))       +                              
  scale_x_discrete(expand=c(0,0)) +                                   
  labs(y="", x="") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), axis.line = element_blank(), legend.title= element_text(size=11), legend.text = element_text(size=10), axis.text = element_text(size=11), legend.text.align = 1, axis.text.y = element_text(face= ifelse(is.element(pca_loadings_pos_hnscc, c(intersect$emt_vitro$HNSCC.PEMT,intersect$snc_epi_vitro$`HNSCC.Epidif.1`)), "bold","plain")))          
 dev.off()
                                         
                 
pca_loadings_neg_hnscc <- rev(melt(apply(expr_emt_episen_pca$rotation[,1:5], 2, function(x) names(sort(x, decreasing=F))[1:10]))$value)              
                                 
pca_loadings_neg_hnscc_plot <- as.matrix(expr_emt_episen_pca$rotation[as.character(pca_loadings_neg_hnscc),1:5] )
rownames(pca_loadings_neg_hnscc_plot) <- NULL
pca_loadings_neg_hnscc_plot <- melt(as.matrix(pca_loadings_neg_hnscc_plot))

pdf("hnscc_pca_neg.pdf", width = 4.5, height = 7.5)                                     
ggplot(data = pca_loadings_neg_hnscc_plot, aes(x=Var2, y=factor(Var1), fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(-0.12, 0.12), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="PC loadings")+
  scale_color_gradient2(limits=c(-0.12, 0.12), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="PC loadings") +
  scale_y_discrete(labels=pca_loadings_neg_hnscc)  +     
  geom_hline(yintercept = seq(10.5,40.5,10)) + 
  geom_vline(xintercept = c(1.5,2.5,3.5, 4.5))       +                              
  scale_x_discrete(expand=c(0,0)) +                                   
  labs(y="", x="") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), axis.line = element_blank(), legend.title= element_text(size=11), legend.text = element_text(size=10), axis.text = element_text(size=11), legend.text.align = 1, axis.text.y = element_text(face= ifelse(is.element(pca_loadings_neg_hnscc, c(intersect$emt_vitro$HNSCC.PEMT,intersect$snc_epi_vitro$`HNSCC.Epidif.1`)), "bold","plain")))  
dev.off()

                                                                                                
# mitf and mel                                  
expr_mitf_emt_plot <- data.frame(expr_mitf_emt_pca$x[,1:5], type=rep(c("in vitro", "in vivo"), c(ncol(expr_mitf_emt_vitro), ncol(expr_mitf_emt_vivo))), emt=colMeans(expr_mitf_emt_final[is.element(rownames(expr_mitf_emt_final), intersect$mel_emt_vitro$melanoma.AXL.program),]), mitf=colMeans(expr_mitf_emt_final[is.element(rownames(expr_mitf_emt_final), intersect$mel_mitf_vitro$melanoma.MITF.program),]))    

expr_mitf_emt_plot$names <- c(rep(mitf_emt_programs_ccle, sapply(expr_ccle[mitf_emt_programs_ccle],ncol)),rep(mitf_emt_programs_tumor, sapply(expr_tumors[mitf_emt_programs_tumor],ncol)) )
expr_mitf_emt_plot$cell_line <- NA
expr_mitf_emt_plot$cell_line[expr_mitf_emt_plot$type=="in vitro"] <- expr_mitf_emt_plot$names[expr_mitf_emt_plot$type=="in vitro"] 
expr_mitf_emt_plot$cell_line <-   gsub("_SKIN", "", expr_mitf_emt_plot$cell_line)   
expr_mitf_emt_plot$tumor <- NA
expr_mitf_emt_plot$tumor[expr_mitf_emt_plot$type=="in vivo"] <- expr_mitf_emt_plot$names[expr_mitf_emt_plot$type=="in vivo"]      
                                    
                                                                                                     
pdf("mitf_emt_pca_program.pdf", width = 3.9, height = 3.8)  
ggplot(expr_mitf_emt_plot, aes(x=PC2+PC3, y=PC4, color=emt-mitf)) +
            geom_point( size=0.5) +
            scale_color_gradient2(limits=c(-3, 3), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(EMT II - Episen)") +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))
dev.off()
                                                                                                                                         
pdf("mitf_emt_pca_type.pdf", width = 4, height = 3.7)  
ggplot(expr_mitf_emt_plot, aes(x=PC2+PC3, y=PC4, color=type)) +
            geom_point( size=0.5) +
             scale_color_manual(values=c("goldenrod", "mediumpurple"), name="") +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
            guides(color=guide_legend(override.aes = list(size=4)))  
dev.off()

                                         
                                         
pdf("mitf_emt_pca_cellline.pdf", width = 3.9, height = 3.8)  
ggplot(expr_mitf_emt_plot, aes(x=PC2+PC3, y= PC4, color=cell_line)) +
            geom_point(data=expr_mitf_emt_plot[is.na(expr_mitf_emt_plot$cell_line),], aes(x=PC2+PC3, y= PC4), size=0.5, color="gray90") + 
             geom_point( size=0.5) +
             scale_color_brewer(palette = "Set2", name="", breaks=unique(expr_mitf_emt_plot$cell_line)[!is.na(unique(expr_mitf_emt_plot$cell_line))]) +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=10.4), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
            guides(color=guide_legend(override.aes = list(size=3), nrow=2) )           
dev.off()
                                                                                                       
pdf("mitf_emt_pca_tumor.pdf", width = 3.9, height = 3.8) 
ggplot(expr_mitf_emt_plot, aes(x=PC2+PC3, y= PC4, color=tumor)) +
            geom_point( size=0.5) +
             scale_color_brewer(palette = "Set2", name="", na.value="gray90", breaks=unique(expr_mitf_emt_plot$tumor)[!is.na(unique(expr_mitf_emt_plot$tumor))], labels=c("Mel 59","Mel 71", "Mel 78", "Mel 19", "Mel 80", "Mel 81", "Mel 88", "Mel 89" )) +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
            guides(color=guide_legend(override.aes = list(size=3), nrow=2))                                              
dev.off()                                   
                                                                                                                                         
  
pdf("mitf_emt_pca_programII.pdf", width = 3.9, height = 3.8)  
ggplot(expr_mitf_emt_plot, aes(x=PC2, y=PC4, color=emt-mitf)) +
            geom_point( size=0.5) +
            scale_color_gradient2(limits=c(-3, 3), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(EMT II - Episen)") +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))
dev.off()
                                                                                                                                         
pdf("mitf_emt_pca_typeII.pdf", width = 3.95, height = 3.7)  
ggplot(expr_mitf_emt_plot, aes(x=PC2, y= PC4, color=type)) +
            geom_point( size=0.5) +
            scale_color_manual(values=c("goldenrod", "mediumpurple"), name="") +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
            guides(color=guide_legend(override.aes = list(size=4)))  
dev.off()                                                                                                                                         
                                                                                        pdf("mitf_emt_pca_programIII.pdf", width = 3.9, height = 3.8)  
ggplot(expr_mitf_emt_plot, aes(x=PC3, y=PC4, color=emt-mitf)) +
            geom_point( size=0.5) +
            scale_color_gradient2(limits=c(-3, 3), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(EMT II - Episen)") +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))
dev.off()
                                                                                                                                         
pdf("mitf_emt_pca_typeIII.pdf", width = 3.95, height = 3.7)  
ggplot(expr_mitf_emt_plot, aes(x=PC3, y= PC4, color=type)) +
            geom_point( size=0.5) +
             scale_color_manual(values=c("goldenrod", "mediumpurple"), name="") +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
            guides(color=guide_legend(override.aes = list(size=4)))  
dev.off()                                                                                                                                               
                                                                                           
pca_loadings_pos_mel <- rev(melt(apply(expr_mitf_emt_pca$rotation[,1:5], 2, function(x) names(sort(x, decreasing=T))[1:10]))$value)              
                                 
pca_loadings_pos_mel_plot <- as.matrix(expr_mitf_emt_pca$rotation[as.character(pca_loadings_pos_mel),1:5] )
rownames(pca_loadings_pos_mel_plot) <- NULL
pca_loadings_pos_mel_plot <- melt(as.matrix(pca_loadings_pos_mel_plot))

pdf("mel_pca_pos.pdf", width = 4.5, height = 7.5)                                     
ggplot(data = pca_loadings_pos_mel_plot, aes(x=Var2, y=factor(Var1), fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(-0.12, 0.12), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="PC loadings")+
  scale_color_gradient2(limits=c(-0.12, 0.12), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="PC loadings") +
  scale_y_discrete(labels=pca_loadings_pos_mel)  +     
  geom_hline(yintercept = seq(10.5,40.5,10)) + 
  geom_vline(xintercept = c(1.5,2.5,3.5, 4.5))       +                              
  scale_x_discrete(expand=c(0,0)) +                                   
  labs(y="", x="") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), axis.line = element_blank(), legend.title= element_text(size=11), legend.text = element_text(size=10), axis.text = element_text(size=11), legend.text.align = 1, axis.text.y = element_text(face= ifelse(is.element(pca_loadings_pos_mel, c(intersect$mel_emt_vitro$melanoma.AXL.program,intersect$mel_mitf_vitro$melanoma.MITF.program)), "bold","plain"))) +
  scale_x_discrete(breaks= unique(Var2)[seq(500,300,500)], labels=seq(500,300,500))                                     
 dev.off()
                                         
                 
pca_loadings_neg_mel <- rev(melt(apply(expr_mitf_emt_pca$rotation[,1:5], 2, function(x) names(sort(x, decreasing=F))[1:10]))$value)              
                                 
pca_loadings_neg_mel_plot <- as.matrix(expr_mitf_emt_pca$rotation[as.character(pca_loadings_neg_mel),1:5] )
rownames(pca_loadings_neg_mel_plot) <- NULL
pca_loadings_neg_mel_plot <- melt(as.matrix(pca_loadings_neg_mel_plot))

pdf("mel_pca_neg.pdf", width = 4.5, height = 7.5)                                     
ggplot(data = pca_loadings_neg_mel_plot, aes(x=Var2, y=factor(Var1), fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(-0.12, 0.12), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="PC loadings")+
  scale_color_gradient2(limits=c(-0.12, 0.12), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="PC loadings") +
  scale_y_discrete(labels=pca_loadings_neg_mel)  +     
  geom_hline(yintercept = seq(10.5,40.5,10)) + 
  geom_vline(xintercept = c(1.5,2.5,3.5, 4.5))       +                              
  scale_x_discrete(expand=c(0,0)) +                                   
  labs(y="", x="") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), axis.line = element_blank(), legend.title= element_text(size=11), legend.text = element_text(size=10), axis.text = element_text(size=11), legend.text.align = 1, axis.text.y = element_text(face= ifelse(is.element(pca_loadings_neg_mel, c(intersect$mel_emt_vitro$melanoma.AXL.program,intersect$mel_mitf_vitro$melanoma.MITF.program)), "bold","plain")))  
dev.off()                                         
                                         
                                         
                                         
### plotting heatmap

                                       
# episen emt
                                       
emt_episen_heatmap <- cbind(expr_emt_episen_vitro,expr_emt_episen_vivo)    
emt_episen_heatmap <- emt_episen_heatmap[c(intersect$emt_vitro$HNSCC.PEMT,intersect$snc_epi_vitro$`HNSCC.Epidif.1`),]
                                      
                                       
emt_episen_heatmap_order <- order(colMeans(emt_episen_heatmap[1:25,]) - colMeans(emt_episen_heatmap[26:54,]) )                                      
emt_episen_heatmap <- emt_episen_heatmap[hclust(as.dist(1-cor(t(emt_episen_heatmap))))$order,emt_episen_heatmap_order]                                       
                        
emt_episen_heatmap_melt <- melt(as.matrix(emt_episen_heatmap))                                               
                                               
p1 <- ggplot(data = emt_episen_heatmap_melt, aes(x=Var2, y=Var1, fill=value, color=value)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(-3, 3), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="Relative\nexpression\n(log2)")+
  scale_color_gradient2(limits=c(-3, 3), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="Relative\nexpression\n(log2)") +                           
  labs(y="", x="Cells") +
  geom_hline(yintercept=29.5) +       
   scale_x_discrete(breaks=unique(emt_episen_heatmap_melt$Var2)[seq(500,5000,500)], labels=seq(500,5000,500))       +                         
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), panel.border = element_rect(fill=F), legend.text = element_text(hjus=1), legend.title = element_text(hjust = 0))   
                                       
                                    

                                       
emt_episen_heatmap_vivo_label  <- data.frame(label= expr_emt_episen_plot$type[emt_episen_heatmap_order] )                                  
      
p2 <- ggplot(data = emt_episen_heatmap_vivo_label, aes(x=1:nrow(emt_episen_heatmap_vivo_label), y="", fill=label, color=label)) + 
  geom_tile() + 
 scale_color_manual(values=c("goldenrod", "mediumpurple"), name="") +
  scale_fill_manual(values=c("goldenrod", "mediumpurple"), name="")  +
  scale_x_continuous(expand= c(0,0))   +
  theme(panel.background = element_blank(), axis.text = element_blank(), axis.ticks=element_blank(), axis.line=element_blank(), axis.title = element_blank(), plot.margin = unit(c(1,1, 0,1), "cm"), legend.position = "top", legend.text = element_text(size=12) )                                    
                                       
                                       
emt_episen_heatmap_vivo_gene <- data.frame(gene= rep(NA, nrow(emt_episen_heatmap)))
                                                     
 emt_episen_heatmap_vivo_gene$gene[ is.element(rownames(emt_episen_heatmap), c("SLPI", "S100A9", "AQP3", "CLDN4", "VIM", "PDPN", "LAMC2", "MYH9"))] <-  rownames(emt_episen_heatmap)[is.element(rownames(emt_episen_heatmap), c("SLPI", "S100A9", "AQP3", "CLDN4", "VIM", "PDPN", "LAMC2", "MYH9"))]    
                                       
p3 <-ggplot(emt_episen_heatmap_vivo_gene, aes(y=1:nrow(emt_episen_heatmap_vivo_gene), x="", label=gene))+
  geom_text_repel(nudge_x = -0.2,  direction = "y", hjust=1) +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title=element_blank(), plot.margin = unit(c(0,-5, 0,0), "cm"), panel.background = element_blank()) +
  scale_y_continuous(expand = c(0,0), position = "right")  
                                       
pdf("hnscc_heatmap_vivo&vitro.pdf", width = 10, 3, onefile = F)                                       
egg::ggarrange(ggplot()+theme(panel.background=element_blank(),  plot.margin = unit(c(0,-5, 0,0), "cm")), p2, p3,p1, nrow=2, ncol=2, widths=c(4,10), heights=c(1.5,10))                                       
dev.off()                                       
                                       

# mitf emt
                                       
mitf_emt_heatmap <- cbind(expr_mitf_emt_vitro,expr_mitf_emt_vivo)    
mitf_emt_heatmap <- mitf_emt_heatmap[c(intersect$mel_mitf_vitro$melanoma.MITF.program,intersect$mel_emt_vitro$melanoma.AXL.program),]
                                      
                                       
mitf_emt_heatmap_order <- order( colMeans(mitf_emt_heatmap[24:38,]) - colMeans(mitf_emt_heatmap[1:23,]))                                      
mitf_emt_heatmap <- mitf_emt_heatmap[rev(hclust(as.dist(1-cor(t(mitf_emt_heatmap))))$order),mitf_emt_heatmap_order]                                       
                        
mitf_emt_heatmap_melt <- melt(as.matrix(mitf_emt_heatmap))                                               
                                               
p4 <- ggplot(data = mitf_emt_heatmap_melt, aes(x=Var2, y=Var1, fill=value, color=value)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(-3, 3), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="Relative\nexpression\n(log2)")+
  scale_color_gradient2(limits=c(-3, 3), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F"),   oob=squish, name="Relative\nexpression\n(log2)") +                           
  labs(y="", x="Cells") +
  geom_hline(yintercept=23.5) +        scale_x_discrete(breaks=unique(mitf_emt_heatmap_melt$Var2)[seq(500,5000,500)], labels=seq(500,5000,500))       +                         
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), panel.border = element_rect(fill=F), legend.text = element_text(hjus=1), legend.title = element_text(hjust = 0))   
                                       
                                    
                                       
mitf_emt_heatmap_vivo_label  <- data.frame(label= expr_mitf_emt_plot$type[mitf_emt_heatmap_order] )                                  
      
p5 <- ggplot(data = mitf_emt_heatmap_vivo_label, aes(x=1:nrow(mitf_emt_heatmap_vivo_label), y="", fill=label, color=label)) + 
  geom_tile() + 
 scale_color_manual(values=c("goldenrod", "mediumpurple"), name="") +
  scale_fill_manual(values=c("goldenrod", "mediumpurple"), name="")  +
  scale_x_continuous(expand= c(0,0))   +
  theme(panel.background = element_blank(), axis.text = element_blank(), axis.ticks=element_blank(), axis.line=element_blank(), axis.title = element_blank(), plot.margin = unit(c(1,1, 0,1), "cm"), legend.position = "top", legend.text = element_text(size=12) )                                    
                                       
                                       
mitf_emt_heatmap_vivo_gene <- data.frame(gene= rep(NA, nrow(mitf_emt_heatmap)))
                                                     
 mitf_emt_heatmap_vivo_gene$gene[ is.element(rownames(mitf_emt_heatmap), c("MITF", "MLANA", "PMEL", "IL8", "PLAUR", "NGFR", "TNFRSF12A"))] <-  rownames(mitf_emt_heatmap)[is.element(rownames(mitf_emt_heatmap), c("MITF", "MLANA", "PMEL", "IL8", "PLAUR", "NGFR", "TNFRSF12A"))]    
                                       
p6 <-ggplot(mitf_emt_heatmap_vivo_gene, aes(y=1:nrow(mitf_emt_heatmap_vivo_gene), x="", label=gene))+
  geom_text_repel(nudge_x = -0.2,  direction = "y", hjust=1) +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title=element_blank(), plot.margin = unit(c(0,-5, 0,0), "cm"), panel.background = element_blank()) +
  scale_y_continuous(expand = c(0,0), position = "right")  
                                       
pdf("mel_heatmap_vivo&vitro.pdf", width = 10, 3, onefile = F)                                       
egg::ggarrange(ggplot()+theme(panel.background=element_blank(),  plot.margin = unit(c(0,-5, 0,0), "cm")), p5, p6,p4, nrow=2, ncol=2, widths=c(4,10), heights=c(1,10))                                       
dev.off()                           
                                       
  
                                       
                                       
                                       
                                       
                                       
                                       
                                       
                                       
                                       
                                       
                                       
                                       
                                       
                                       
                                       
                                       
                                       
                                       
                                       
                                       
                            
#########################################################################################################################
################################ plotting in vivo vs in vitro programs - heatmap ###################################################

expr <- readRDS("/home/labs/tirosh/kinker/analyses/CCLE/data processing/scCCLE&pool1_tpm_by_cell_line.rds") 
common_genes <- Reduce(intersect, lapply(expr, rownames))
expr <- lapply(expr, function(x) x[common_genes,])

### calculating average TPM values
ave_tpm <- rowMeans(do.call(cbind, expr))

### processing TPM data
expr <- lapply(expr, function(x) log2((x/10) + 1))
ave_log_tpm <- rowMeans(do.call(cbind, expr))
expr <- lapply(expr, function(x) x-ave_log_tpm)

### reading ccle nmf programs (already ordered by hc)
nmf_ccle <- readRDS(file.choose())

### reading metaprograms
# literature
meta_vivo_literature <- read.table(file.choose(), sep = "\t", header = T,stringsAsFactors = F) 
meta_vivo_literature <- meta_vivo_literature[1:50,]

# nmf
meta_vivo_nmf <-  read.table(file.choose(), sep = "\t", header = T,stringsAsFactors = F)

# combine
vivo_signatures <- unlist(unname(list(meta_vivo_literature,meta_vivo_nmf)), recursive = F )
vivo_signatures <- lapply(vivo_signatures, function(x) x[x!=""])

### calculating the similarity between in vitro and in vivo
intersect_vivo <- sapply(vivo_signatures, function(x) apply(nmf_ccle, 2, function(y) length(intersect(x,y))/length(union(x,y))))

### clustering
hc_intersect_vivo <- hclust(as.dist(1-cor(intersect_vivo)), method="average") 
hc_intersect_vivo <- reorder(as.dendrogram(hc_intersect_vivo), colMeans(intersect_vivo))

### plotting intersection
intersect_vivo_plot <-  reshape2::melt(intersect_vivo[,order.dendrogram(hc_intersect_vivo)]) 

gabi_palette <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

p1 <- ggplot(data = intersect_vivo_plot, aes(x=Var2, y=Var1, fill=value*100, color=value*100)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(2,25), low=gabi_palette[1:111],  mid =gabi_palette[112:222], high = gabi_palette[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  scale_color_gradient2(limits=c(2,25), low=gabi_palette[1:111],  mid =gabi_palette[112:222], high = gabi_palette[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  theme( axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank()) + 
  guides(fill = F, color=F) +
  panel_border(colour = "black", size = 0.4, linetype = 1, remove = FALSE) +
  geom_hline(yintercept = c(25,53, 75, 95, 109, 132, 167, 250, 278), size=0.3)


### calculating program scores - in vitro
nmf_ccle_filt <- apply(nmf_ccle, 2, function(x) x[is.element(x, common_genes)])
nmf_ccle_filt_ctr <- lapply(nmf_ccle_filt, function(x) get_control_geneset(ave_tpm = ave_tpm, program = x, bins = 42))

nmf_ccle_scores <- list()

for(i in names(nmf_ccle_filt)) {
  a <- gsub(".{4}$", "", i)
  b <- expr[[a]]
  nmf_ccle_scores[[i]] <- colMeans(b[ nmf_ccle_filt[[i]],]) - colMeans(b[nmf_ccle_filt_ctr[[i]],])
}

### calculating program scores - in vivo
vivo_signatures_filt <- lapply(vivo_signatures, function(x) x[is.element(x, common_genes)])
vivo_signatures_filt_crt <- lapply(vivo_signatures_filt, function(x) get_control_geneset(ave_tpm = ave_tpm, program = x, bins = 42))

in_vivo_scores <- list()
for(i in names(nmf_ccle_filt)) {
  a <- gsub(".{4}$", "", i)
  b <- expr[[a]]
  in_vivo_scores[[i]] <- sapply(vivo_signatures_filt, function(x) colMeans(b[x,])) - sapply(vivo_signatures_filt_crt, function(x) colMeans(b[x,]))
}

### calculating correlation vitro vs vivo
corr_vitro_vivo <- sapply(colnames(nmf_ccle), function(x) cor(nmf_ccle_scores[[x]], in_vivo_scores[[x]]))

### calculating random maximum correlation and maximum similarity
max_sim_perm <- c()
max_cor_perm <- c()

for(i in 1:100) {
  a <- lapply(nmf_ccle_filt, function(x) get_control_geneset(ave_tpm = ave_tpm, program = x, bins = 42, size=1, seed = i))
  b <- list()
  for(j in names(a)) {
    c <- gsub(".{4}$", "", j)
    d <- expr[[c]]
    b[[j]] <- colMeans(d[a[[j]],])
  }
  max_sim_perm <- c(max_sim_perm, max(sapply(vivo_signatures, function(x) sapply(a, function(y) length(intersect(x,y))/length(union(x,y))))))
  max_cor_perm <- c(max_cor_perm, max(sapply(names(a), function(x) cor(b[[x]], in_vivo_scores[[x]]))))
}

### plotting correlations 
max_cor_plot <- data.frame("max_sim"=apply(corr_vitro_vivo,2,max))

p2 <- ggplot(max_cor_plot, aes(y=max_sim, x=1:nrow(max_cor_plot))) +
  geom_smooth(span=0.04, se=F, size=0.8) + 
  coord_flip() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(position = "right") +
  theme(axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(2.5,1.5,0.5,0), "cm"), axis.text.x = element_text(size=11), axis.title.x = element_text(size=12)) +
  panel_border(colour="black", size=0.4) +
  labs(y="Maximum\ncorrelation") +
  geom_hline(yintercept = max(max_cor_perm), linetype="dashed")

egg::ggarrange(p1, p2, ncol=2, widths = c(6,2.5))

pdf("nc_vitro_vs_vivo.pdf", width = 4.3, height = 7.1, onefile = F)

### ploting maximum correlations/similarities of in vivo programs
sim_corr_plot_vivo <- data.frame("max_sim"= apply(intersect_vivo, 2, max), "max_cor"=apply(corr_vitro_vivo,1,max))
sim_corr_plot_vivo$label <- NA
sim_corr_plot_vivo$label[order(sim_corr_plot_vivo$max_cor, decreasing = T)[1:10]] <- 1:10 

sim_corr_plot_vivo["HNSCC.Epidif.1","label"] <- "HNSCC Epi. dif. 1"
sim_corr_plot_vivo["HNSCC.Stress","label"] <- "HNSCC Stress"
sim_corr_plot_vivo["melanoma.MITF.program","label"] <- "Melanoma MITF"
sim_corr_plot_vivo["HNSCC.PEMT","label"] <- "HNSCC pEMT"
sim_corr_plot_vivo["melanoma.Inflammatory","label"] <- "Melanoma stress"
sim_corr_plot_vivo["melanoma.AXL.program","label"] <- "Melanoma AXL"
sim_corr_plot_vivo["snc_epi","label"] <- "Meta epi. sen."
sim_corr_plot_vivo["stress","label"] <- "Meta stress"
sim_corr_plot_vivo["pemt","label"] <- "Meta pEMT"
sim_corr_plot_vivo["GBM.MES2","label"] <- "GBM mes 2"


ggplot(sim_corr_plot_vivo, aes(x=max_sim, y=max_cor, label=label)) +
  geom_point(size=5, shape=21, fill="gray95") +
  geom_hline(yintercept = max(max_cor_perm), linetype="dashed") +
  geom_vline(xintercept = max(max_sim_perm), linetype="dashed") +
  geom_text_repel(size=3, point.padding = 0, box.padding = 0.1) +
  labs(y="Maximum correlation", x= "Maximum similarity\n(Jaccard index)") +
  theme(axis.text = element_text(size=11), axis.title = element_text(size=12), axis.line = element_blank(), plot.margin = unit(c(1,8,1,1), "cm")) +
  panel_border(colour = "black", size = 0.4, linetype = 1, remove = FALSE)+
  ylim(0.34, 1)

pdf("nc_vivo_sim_cor.pdf", width = 6.8, height = 4, onefile = F)

### plotting maximum correlations/similarities of in vitro metaprogram by in vivo program

# getting similarities/correlations for the selected in vivo programs
sim_corr_vitro <- c()

for(i in c( "melanoma.MITF.program", "melanoma.AXL.program", "HNSCC.Epidif.1", "HNSCC.PEMT", "melanoma.Inflammatory", "HNSCC.Stress")) {
  a <- data.frame("sim"= intersect_vivo[,i], "cor"=corr_vitro_vivo[which(names(vivo_signatures) == i),])
  a$facet <- i
  sim_corr_vitro <- rbind(sim_corr_vitro,a)
}

sim_corr_vitro$facet <-factor(sim_corr_vitro$facet, levels = c( "melanoma.MITF.program", "melanoma.AXL.program","melanoma.Inflammatory", "HNSCC.Epidif.1", "HNSCC.PEMT", "HNSCC.Stress"))
sim_corr_vitro$order <- rep(1:432, 6)

# annotation the in vitro programs
sim_corr_vitro$metaprogram <- NA
sim_corr_vitro$metaprogram[is.element(sim_corr_vitro$order, 10:25)] <- "Skin\npigmenation"
sim_corr_vitro$metaprogram[is.element(sim_corr_vitro$order, 46:53)] <- "EMT I"
sim_corr_vitro$metaprogram[is.element(sim_corr_vitro$order, 57:75)] <- "EMT II"
sim_corr_vitro$metaprogram[is.element(sim_corr_vitro$order, 81:95)] <- "IFN\nresponse"
sim_corr_vitro$metaprogram[is.element(sim_corr_vitro$order, 98:109)] <- "EMT III"
sim_corr_vitro$metaprogram[is.element(sim_corr_vitro$order, 116:132)] <- "p53-dep.\nsenescence"
sim_corr_vitro$metaprogram[is.element(sim_corr_vitro$order, 140:167)] <- "Epi.\nsenescence"
sim_corr_vitro$metaprogram[is.element(sim_corr_vitro$order, 197:250)] <- "Stress"
sim_corr_vitro$metaprogram[is.element(sim_corr_vitro$order, 261:278)] <- "Protein\nmaturation"
sim_corr_vitro$metaprogram[is.element(sim_corr_vitro$order, 279:432)] <- "Proteasomal\ndegradation"

sim_corr_vitro <- sim_corr_vitro[!is.na(sim_corr_vitro$metaprogram),]

# averaging correlations/similarities per metaprogram
sim_corr_plot_vitro_mean <- c()

for(i in unique(sim_corr_vitro$facet)) {
  a <- sim_corr_vitro[sim_corr_vitro$facet==i,]
  b <- cbind(aggregate(a$sim, list(a$metaprogram), mean), aggregate(a$cor, list(a$metaprogram), mean)[,2])
  b$facet <- i
  sim_corr_plot_vitro_mean <- rbind(sim_corr_plot_vitro_mean,b)
}

colnames(sim_corr_plot_vitro_mean) <- c("metaprogram", "mean_sim", "mean_cor", "facet")

# plot
sim_corr_plot_vitro_mean$facet <-factor(sim_corr_plot_vitro_mean$facet, levels = c( "melanoma.MITF.program", "melanoma.AXL.program","melanoma.Inflammatory", "HNSCC.Epidif.1", "HNSCC.PEMT", "HNSCC.Stress"))

ggplot(sim_corr_plot_vitro_mean, aes(x=mean_sim, y=mean_cor, fill=metaprogram)) +
  geom_point(size=4.5, shape=21) +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(facets = vars(facet), nrow = 2, labeller = as_labeller(c("melanoma.MITF.program"="Melanoma MITF",  "melanoma.AXL.program"= "Melanoma AXL","melanoma.Inflammatory"= "Melanoma Stress", "HNSCC.Epidif.1"="HNSCC Epi. dif. 1", "HNSCC.PEMT"= "HNSCC pEMT", "HNSCC.Stress"="HNSCC Stress"))) +
  theme(strip.background = element_blank(), legend.title = element_blank(), axis.line = element_blank()) +
  labs(x="Mean similarity\n(Jaccard Index)", y="Mean correlation") +
  panel_border(colour = "black", size = 0.4, linetype = 1, remove = FALSE) 

pdf("nc_vitro_sim_cor.pdf", width = 7.4, height = 4.4, onefile = F)


                                                                                
                                             
                                             
#########################################################################################################################
################################ plotting in vivo vs in vitro programs - scatterplots ###################################################


### reading TPM expression data
expr <- readRDS("/home/labs/tirosh/kinker/analyses/CCLE/data processing/scCCLE&pool1_tpm_by_cell_line.rds") 
common_genes <- Reduce(intersect, lapply(expr, rownames))
expr <- lapply(expr, function(x) x[common_genes,])

### calculating average TPM values
ave_tpm <- rowMeans(do.call(cbind, expr))

### processing TPM data
expr <- lapply(expr, function(x) log2((x/10) + 1))
ave_log_tpm <- rowMeans(do.call(cbind, expr))
expr <- lapply(expr, function(x) x-ave_log_tpm)

### reading in vitro metaprograms
mel_mitf <- scan("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/metaprograms - top genes/meta_mel_mitf.txt", "character")
mel_emt <- scan("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/metaprograms - top genes/meta_mel_emt.txt", "character")
emt <- scan("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/metaprograms - top genes/meta_emt.txt", "character")
snc_epi <- scan("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/metaprograms - top genes/meta_snc_epi.txt", "character")
stress <- scan("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/metaprograms - top genes/meta_stress.txt", "character")

### reading in vivo metaprograms - nmf
# literature
meta_vivo_literature <- read.table(file.choose(), sep = "\t", header = T,stringsAsFactors = F) 
meta_vivo_literature <- meta_vivo_literature[1:50,]

# nmf
meta_vivo_nmf <-  read.table(file.choose(), sep = "\t", header = T,stringsAsFactors = F)

# combine
vivo_signatures <- unlist(unname(list(meta_vivo_literature,meta_vivo_nmf)), recursive = F )
vivo_signatures <- lapply(vivo_signatures, function(x) x[x!=""])

### list of metaprograms
metaprograms <- c(list("mel_mitf_vitro" = mel_mitf, "mel_emt_vitro" = mel_emt, "emt_vitro"= emt, "snc_epi_vitro"=snc_epi, "stress_vitro"=stress), vivo_signatures)

### filtering genes
metaprograms <- lapply(metaprograms, function(x) x[is.element(x, common_genes)])

### creating control gensets
control_genests <- lapply(metaprograms, function(x) get_control_geneset(ave_tpm = ave_tpm, program = x, bins = 100))

### calculating metaprograms scores
metaprograms_scores <- lapply(names(metaprograms), function(x) lapply(expr, function(y) colMeans(y[metaprograms[[x]],]) - colMeans(y[control_genests[[x]],])))
names(metaprograms_scores) <- names(metaprograms)

### getting model cell lines for each metaprogram
mel_mitf_models <- gsub(".{4}$", "", scan("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/meta_mel_mitf_programs.txt", "character"))
mel_emt_models <-  gsub(".{4}$", "",scan("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/meta_mel_emt_programs.txt", "character"))
emt_models <-  gsub(".{4}$", "",scan("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/meta_emt_programs.txt", "character"))
snc_epi_models <- gsub(".{4}$", "",scan("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/meta_snc_epi_programs.txt", "character"))
stress_models <- gsub(".{4}$", "",scan("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/meta_stress_programs.txt", "character"))


mel_mitf_plot <- data.frame("vitro"=unlist(metaprograms_scores$mel_mitf_vitro$SKMEL30_SKIN), "vivo"= unlist(metaprograms_scores$melanoma.MITF.program$SKMEL30_SKIN))

p1 <- ggplot(mel_mitf_plot, aes(x=vitro, y=vivo)) +
  geom_point(alpha=0.4, color="gray80", size=4) +
  theme(axis.line=element_blank(), plot.title=element_text(hjust = 0)) +
  panel_border(colour = "black", 0.4) +
  labs(x="Skin pigmentation\nprogram score", y = "Melanoma MITF\nprogram score" , title = "SKMEL30 - Melanoma" )

mel_emt_plot <- data.frame("vitro"=unlist(metaprograms_scores$mel_emt_vitro$IGR1_SKIN), "vivo"= unlist(metaprograms_scores$melanoma.AXL.program$IGR1_SKIN))

p2 <- ggplot(mel_emt_plot, aes(x=vitro, y=vivo)) +
  geom_point(alpha=0.4, color="gray80", size=4) +
  theme(axis.line=element_blank(), plot.title=element_text(hjust = 0)) +
  panel_border(colour = "black", 0.4) +
  labs(x="EMT I \nprogram score", y = "Melanoma AXL\nprogram score" , title = "IGR1 - Melanoma" )

emt_plot <-  data.frame("vitro"=unlist(metaprograms_scores$emt_vitro$SCC47_UPPER_AERODIGESTIVE_TRACT), "vivo"= unlist(metaprograms_scores$HNSCC.PEMT$SCC47_UPPER_AERODIGESTIVE_TRACT))

p3 <- ggplot(emt_plot, aes(x=vitro, y=vivo)) +
  geom_point(alpha=0.4, color="gray80", size=4) +
  theme(axis.line=element_blank(), plot.title=element_text(hjust = 0)) +
  panel_border(colour = "black", 0.4) +
  labs(x="EMT II\nprogram score", y = "HNSCC pEMT\nprogram score", title = "SCC47 - Oropharyngeal\n squamous cell carcinoma"   )

epi_snc_plot <- data.frame("vitro"=unlist(metaprograms_scores$snc_epi_vitro$JHU006_UPPER_AERODIGESTIVE_TRACT), "vivo"=unlist(metaprograms_scores$HNSCC.Epidif.1$JHU006_UPPER_AERODIGESTIVE_TRACT))

p4 <- ggplot(epi_snc_plot[epi_snc_plot$vitro<3,], aes(x=vitro, y=vivo)) +
  geom_point(alpha=0.4, color="gray80", size=4) +
  theme(axis.line=element_blank(), plot.title=element_text(hjust = 0)) +
  panel_border(colour = "black", 0.4) +
  labs(x="Epithelial senescence\nprogram score", y = "HNSCC epithelial diff. I\nprogram score", title= "JHU006 - Laryngea\nsquamous cell carcinoma" )

stress_plot <- data.frame("vitro"=unlist(metaprograms_scores$stress_vitro$JHOS2_OVARY), "vivo"= unlist(metaprograms_scores$stress$JHOS2_OVARY))

p5 <- ggplot(stress_plot, aes(x=vitro, y=vivo)) +
  geom_point(alpha=0.4, color="gray80", size=4) +
  theme(axis.line=element_blank(),plot.title=element_text(hjust = 0)) +
  panel_border(colour = "black", 0.4) +
  labs(x="Stress response\nprogram score", y = "Melanoma/HNSCC stress\nprogram score", title =  "JHOS2 - Ovarian cancer") 

egg::ggarrange(p1, p2, p3, p4, p5, ncol=5)

pdf("nc_vitro_vs_vivo_scatter.pdf", width = 18, height = 4, onefile = F)


#########################################################################################################################
################################ plotting in vivo vs in vitro programs - heatmaps ###################################################


### reading TPM expression data - cell lines
expr_ccle <- readRDS("/home/labs/tirosh/kinker/analyses/CCLE/data processing/scCCLE&pool1_tpm_by_cell_line.rds") 
common_genes_ccle <- Reduce(intersect, lapply(expr_ccle, rownames))
expr_ccle <- lapply(expr_ccle, function(x) x[common_genes_ccle,])
expr_ccle <- lapply(expr_ccle, function(x) log2((x/10) + 1))

### reading log TPM expression data - tumors
expr_tumors <- readRDS("/home/labs/tirosh/kinker/datasets/scRNA-seq/scRNAseq_tumors_logtpm.RDS")
expr_tumors <- unlist(expr_tumors, recursive = F)
expr_tumors <- expr_tumors[sapply(expr_tumors, function(x) nrow(x) >= 50)]
expr_tumors <- lapply(expr_tumors, t)

### reading in vitro metaprograms
mel_mitf <- read.table("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/metatable_mel_mitf.txt", sep="\t", header = T, stringsAsFactors = F)[1:100,1]
mel_emt <- read.table("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/metatable_mel_emt.txt", sep="\t", header = T, stringsAsFactors = F)[1:100,1]
emt <- read.table("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/metatable_emt.txt", sep="\t", header = T, stringsAsFactors = F)[1:100,1]
snc_epi <- read.table("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/metatable_snc_epi.txt", sep="\t", header = T, stringsAsFactors = F)[1:100,1]
stress <- read.table("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/metatable_stress.txt", sep="\t", header = T, stringsAsFactors = F)[1:100,1]

### combine
vitro_signatures <- c(list("mel_mitf_vitro" = mel_mitf, "mel_emt_vitro" = mel_emt, "emt_vitro"= emt, "snc_epi_vitro"=snc_epi, "stress_vitro"=stress))
vitro_signatures <- lapply(vitro_signatures, function(x) x[is.element(x, common_genes_ccle)])

### reading in vivo metaprograms - nmf
# literature
meta_vivo_literature <- read.table(file.choose(), sep = "\t", header = T,stringsAsFactors = F) 
meta_vivo_literature <- meta_vivo_literature[1:100,]

# nmf
stress_nmf_vivo <- read.table("/home/labs/tirosh/kinker/analyses/CCLE/nmf/in vivo/tumors/metatable_stress_vivo.txt", sep="\t", header = T, stringsAsFactors = F)[1:100,1]

### combine
vivo_signatures <- unlist(unname(list(meta_vivo_literature,list("stress"=stress_nmf_vivo))), recursive = F )
vivo_signatures <- lapply(vivo_signatures, function(x) x[x!=""])

### intersection
intersect <- lapply(vitro_signatures, function(x) lapply(vivo_signatures, function(y) intersect(x,y)))

### plotting
# melanoma
SKMEL30_plot <- expr_ccle$SKMEL30_SKIN[c(intersect$mel_mitf_vitro$melanoma.MITF.program, intersect$mel_emt_vitro$melanoma.AXL.program),]
SKMEL30_plot <- SKMEL30_plot - rowMeans(SKMEL30_plot)
SKMEL30_plot <- SKMEL30_plot[,order(colMeans(SKMEL30_plot[intersect$mel_mitf_vitro$melanoma.MITF.program,]) - colMeans(SKMEL30_plot[intersect$mel_emt_vitro$melanoma.AXL.program,]))]
colnames(SKMEL30_plot) <- 1:ncol(SKMEL30_plot)
SKMEL30_plot <- melt(SKMEL30_plot)

IGR1_plot <- expr_ccle$IGR1_SKIN[c(intersect$mel_mitf_vitro$melanoma.MITF.program, intersect$mel_emt_vitro$melanoma.AXL.program),]
IGR1_plot <- IGR1_plot - rowMeans(IGR1_plot)
IGR1_plot <- IGR1_plot[,order(colMeans(IGR1_plot[intersect$mel_mitf_vitro$melanoma.MITF.program,]) - colMeans(IGR1_plot[intersect$mel_emt_vitro$melanoma.AXL.program,]))]
colnames(IGR1_plot) <- 1:ncol(IGR1_plot)
IGR1_plot <- melt(IGR1_plot)

meltumor79_plot <- expr_tumors$melanoma_tirosh_2016.79[c(intersect$mel_mitf_vitro$melanoma.MITF.program, intersect$mel_emt_vitro$melanoma.AXL.program),]
meltumor79_plot <- meltumor79_plot - rowMeans(meltumor79_plot)
meltumor79_plot <- meltumor79_plot[,order(colMeans(meltumor79_plot[intersect$mel_mitf_vitro$melanoma.MITF.program,])-colMeans(meltumor79_plot[intersect$mel_emt_vitro$melanoma.AXL.program,]))]
colnames(meltumor79_plot) <- 1:ncol(meltumor79_plot)
meltumor79_plot <- melt(meltumor79_plot)

meltumor80_plot <- expr_tumors$melanoma_tirosh_2016.80[c(intersect$mel_mitf_vitro$melanoma.MITF.program, intersect$mel_emt_vitro$melanoma.AXL.program),]
meltumor80_plot <- meltumor80_plot - rowMeans(meltumor80_plot)
meltumor80_plot <- meltumor80_plot[,order(colMeans(meltumor80_plot[intersect$mel_mitf_vitro$melanoma.MITF.program,])-colMeans(meltumor80_plot[intersect$mel_emt_vitro$melanoma.AXL.program,]))]
colnames(meltumor80_plot) <- 1:ncol(meltumor80_plot)
meltumor80_plot <- melt(meltumor80_plot)

melfinal_plot <- rbind(SKMEL30_plot,IGR1_plot, meltumor79_plot, meltumor80_plot)
melfinal_plot$facet <- rep(c("SKMEL30\n", "IGR1\n",  "Melanoma 79\n", "Melanoma 80\n"), c(12355, 25375,16380, 4375 ))
melfinal_plot$facet <- factor(melfinal_plot$facet, levels = c("SKMEL30\n", "IGR1\n", "Melanoma 79\n", "Melanoma 80\n"))

count <- 0
p1 <- ggplot(data = melfinal_plot, aes(y=Var1, x=as.numeric(Var2))) + 
  geom_tile(aes(fill=value, color=value)) + 
  scale_fill_gradient2(limits=c(-4, 4), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "#F7F7F7", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative\nexpression (log2)") +
  scale_color_gradient2(limits=c(-4, 4), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "#F7F7F7", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative\nexpression (log2)") +
  facet_grid(cols=vars(facet), scales = "free_x") +
  scale_x_continuous(breaks= function(x) {
    count <<- count + 1L
    switch(
      count,
      seq(100, 300, 100),
      seq(200, 600, 200),
      seq(100, 500, 100),
      seq(20, 100, 20),
  )
  },  expand = c(0,0)) +
  theme(axis.ticks.y = element_blank(),  axis.line = element_blank(), axis.text.y = element_blank(), strip.background = element_blank(), strip.text = element_text(size=14), legend.text.align = 1, legend.margin = margin(c(1,1,1,150)), panel.spacing = unit(1.2, "lines")) +
  labs(x="Cells", y="") +
  panel_border(colour="black", remove = F, size = 0.4)

labels_melanoma <- data.frame("labels"=c(intersect$mel_mitf_vitro$melanoma.MITF.program, intersect$mel_emt_vitro$melanoma.AXL.program))
labels_melanoma$labels[!is.element(labels_melanoma$labels, c("MITF", "MLANA", "PMEL", "NGFR", "FN1", "PLAUR"))] <- NA
 
p2 <-  ggplot(labels_melanoma, aes(y=1:nrow(labels_melanoma), x="", label=labels))+
  geom_text_repel(nudge_x = -0.05,  direction = "y", hjust=1) +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title=element_blank(), plot.margin = unit(c(1,-8.3, 1,-2), "cm")) +
  scale_y_continuous(expand = c(0,0), position = "right")

pdf("heatmap_vivo_vitro_melanoma.pdf", width = 17, height = 2.8, onefile = F)
egg::ggarrange(p2, p1, ncol = 2, widths = c(5,10))
dev.off()

# hnscc
JHU006_plot <- expr_ccle$JHU006_UPPER_AERODIGESTIVE_TRACT[c(intersect$snc_epi_vitro$HNSCC.Epidif.1, intersect$emt_vitro$HNSCC.PEMT),]
colnames(JHU006_plot) <- paste("V1",colnames(JHU006_plot))
JHU006_plot <- JHU006_plot - rowMeans(JHU006_plot)
JHU006_plot <- JHU006_plot[,order(colMeans(JHU006_plot[intersect$snc_epi_vitro$HNSCC.Epidif.1,]) )]
colnames(JHU006_plot) <- 1:ncol(JHU006_plot)
JHU006_plot <- melt(JHU006_plot)

SCC47_plot <- expr_ccle$SCC47_UPPER_AERODIGESTIVE_TRACT[c(intersect$snc_epi_vitro$HNSCC.Epidif.1, intersect$emt_vitro$HNSCC.PEMT),]
colnames(SCC47_plot) <- paste("V1",colnames(SCC47_plot))
SCC47_plot <- SCC47_plot - rowMeans(SCC47_plot)
SCC47_plot <- SCC47_plot[,order(colMeans(SCC47_plot[intersect$snc_epi_vitro$HNSCC.Epidif.1,]) - colMeans(SCC47_plot[ intersect$emt_vitro$HNSCC.PEMT,]))]
colnames(SCC47_plot) <- 1:ncol(SCC47_plot)
SCC47_plot <- melt(SCC47_plot)

hnscctumor17_plot <- expr_tumors$hnscc_puram_2017.17[c(intersect$snc_epi_vitro$HNSCC.Epidif.1, intersect$emt_vitro$HNSCC.PEMT),]
hnscctumor17_plot <- hnscctumor17_plot - rowMeans(hnscctumor17_plot)
hnscctumor17_plot <- hnscctumor17_plot[,order(colMeans(hnscctumor17_plot[intersect$snc_epi_vitro$HNSCC.Epidif.1,]) - colMeans(hnscctumor17_plot[ intersect$emt_vitro$HNSCC.PEMT,]))]
colnames(hnscctumor17_plot) <- 1:ncol(hnscctumor17_plot)
hnscctumor17_plot <- melt(hnscctumor17_plot)

hnscctumor25_plot <- expr_tumors$hnscc_puram_2017.25[c(intersect$snc_epi_vitro$HNSCC.Epidif.1, intersect$emt_vitro$HNSCC.PEMT),]
hnscctumor25_plot <- hnscctumor25_plot - rowMeans(hnscctumor25_plot)
hnscctumor25_plot <- hnscctumor25_plot[,order(colMeans(hnscctumor25_plot[intersect$snc_epi_vitro$HNSCC.Epidif.1,]) - colMeans(hnscctumor25_plot[ intersect$emt_vitro$HNSCC.PEMT,]))]
colnames(hnscctumor25_plot) <- 1:ncol(hnscctumor25_plot)
hnscctumor25_plot <- melt(hnscctumor25_plot)

hnsccfinal_plot <- rbind(JHU006_plot, SCC47_plot, hnscctumor17_plot, hnscctumor25_plot)
hnsccfinal_plot$facet <- rep(c("JHU006\n", "SCC47\n", "HNSCC 17\n", "HNSCC 25\n"), c(17550, 31968, 17820, 11286 ))
hnsccfinal_plot$facet <- factor(hnsccfinal_plot$facet, levels = c("JHU006\n", "SCC47\n", "HNSCC 17\n", "HNSCC 25\n"))


count <- 0
p3 <- ggplot(data = hnsccfinal_plot, aes(y=Var1, x=as.numeric(Var2))) + 
  geom_tile(aes(fill=value, color=value)) + 
  scale_fill_gradient2(limits=c(-4, 4), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "#F7F7F7", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative\nexpression (log2)") +
  scale_color_gradient2(limits=c(-4, 4), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "#F7F7F7", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative\nexpression (log2)") +
  facet_grid(cols=vars(facet), scales = "free_x") +
  scale_x_continuous(breaks= function(x) {
    count <<- count + 1L
    switch(
      count,
      seq(50, 600, 50),
      seq(100, 600, 100),
      seq(50, 600, 50),
      seq(50, 600, 50),
    )
  },  expand = c(0,0)) +
  theme(axis.ticks.y = element_blank(),  axis.line = element_blank(), axis.text.y = element_blank(), strip.background = element_blank(), strip.text = element_text(size=14), legend.text.align = 1, legend.margin = margin(c(1,1,1,150)), panel.spacing = unit(1.2, "lines")) +
  labs(x="Cells", y="") +
  panel_border(colour="black", remove = F, size = 0.4)


labels_hnscc <- data.frame("labels"=c(intersect$snc_epi_vitro$HNSCC.Epidif.1, intersect$emt_vitro$HNSCC.PEMT))
labels_hnscc$labels[!is.element(labels_hnscc$labels, c(  "S100P" ,"AQP3", "CLDN4", "VIM", "COL5A2", "PDPN"))] <- NA

p4 <- ggplot(labels_hnscc, aes(y=1:nrow(labels_hnscc), x="", label=labels))+
  geom_text_repel(nudge_x = -0.05,  direction = "y", hjust=1) +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title=element_blank(), plot.margin = unit(c(1,-8.3, 1,-2), "cm")) +
  scale_y_continuous(expand = c(0,0), position = "right")

pdf("heatmap_vivo_vitro_hnscc.pdf", width = 17, height = 2.8, onefile = F)
egg::ggarrange(p4, p3, ncol = 2, widths = c(5,10))
dev.off()













ggplot(data = JHU006_plot, aes(y=Var1, x=Var2)) + 
  geom_tile(aes(fill=value, color=value)) + 
  scale_fill_gradient2(limits=c(-4, 4), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "#F7F7F7", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative\nexpression (log2)") +
  scale_color_gradient2(limits=c(-4, 4), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "#F7F7F7", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative\nexpression (log2)") 
  
        
        
        axis.text.y  = element_blank(), axis.text.x = element_text(size=9), legend.title = element_text(size=9, vjust = 0.8), axis.title= element_text(size=10),   legend.text = element_text(size = 8), legend.text.align = 0.5, legend.justification = "bottom", plot.margin = unit(c(3.5,0,3.5,1), "cm"))+ 
  labs(x="\nChromosome", y="") +
  scale_x_discrete( labels= plot_label$label) +
  panel_border(colour="black", remove = F) +
  geom_vline(xintercept = vlines[1:22]) +
  geom_hline(yintercept = 37,  linetype="dashed")


### filtering genes
metaprograms <- lapply(metaprograms, function(x) x[is.element(x, common_genes)])

### creating control gensets
control_genests <- lapply(metaprograms, function(x) get_control_geneset(ave_tpm = ave_tpm, program = x, bins = 100))

### calculating metaprograms scores
metaprograms_scores <- lapply(names(metaprograms), function(x) lapply(expr, function(y) colMeans(y[metaprograms[[x]],]) - colMeans(y[control_genests[[x]],])))
names(metaprograms_scores) <- names(metaprograms)

### getting model cell lines for each metaprogram
mel_mitf_models <- gsub(".{4}$", "", scan("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/meta_mel_mitf_programs.txt", "character"))
mel_emt_models <-  gsub(".{4}$", "",scan("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/meta_mel_emt_programs.txt", "character"))
emt_models <-  gsub(".{4}$", "",scan("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/meta_emt_programs.txt", "character"))
snc_epi_models <- gsub(".{4}$", "",scan("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/meta_snc_epi_programs.txt", "character"))
stress_models <- gsub(".{4}$", "",scan("/home/labs/tirosh/kinker/analyses/CCLE/nmf/metaprograms/meta_stress_programs.txt", "character"))


mel_mitf_plot <- data.frame("vitro"=unlist(metaprograms_scores$mel_mitf_vitro$SKMEL30_SKIN), "vivo"= unlist(metaprograms_scores$melanoma.MITF.program$SKMEL30_SKIN))

p1 <- ggplot(mel_mitf_plot, aes(x=vitro, y=vivo)) +
  geom_point(alpha=0.4, color="gray80", size=4) +
  theme(axis.line=element_blank(), plot.title=element_text(hjust = 0)) +
  panel_border(colour = "black", 0.4) +
  labs(x="Skin pigmentation\nprogram score", y = "Melanoma MITF\nprogram score" , title = "SKMEL30 - Melanoma" )

mel_emt_plot <- data.frame("vitro"=unlist(metaprograms_scores$mel_emt_vitro$IGR1_SKIN), "vivo"= unlist(metaprograms_scores$melanoma.AXL.program$IGR1_SKIN))

p2 <- ggplot(mel_emt_plot, aes(x=vitro, y=vivo)) +
  geom_point(alpha=0.4, color="gray80", size=4) +
  theme(axis.line=element_blank(), plot.title=element_text(hjust = 0)) +
  panel_border(colour = "black", 0.4) +
  labs(x="EMT I \nprogram score", y = "Melanoma AXL\nprogram score" , title = "IGR1 - Melanoma" )

emt_plot <-  data.frame("vitro"=unlist(metaprograms_scores$emt_vitro$SCC47_UPPER_AERODIGESTIVE_TRACT), "vivo"= unlist(metaprograms_scores$HNSCC.PEMT$SCC47_UPPER_AERODIGESTIVE_TRACT))

p3 <- ggplot(emt_plot, aes(x=vitro, y=vivo)) +
  geom_point(alpha=0.4, color="gray80", size=4) +
  theme(axis.line=element_blank(), plot.title=element_text(hjust = 0)) +
  panel_border(colour = "black", 0.4) +
  labs(x="EMT II\nprogram score", y = "HNSCC pEMT\nprogram score", title = "SCC47 - Tongue squamous\ncell carcinoma"   )

epi_snc_plot <- data.frame("vitro"=unlist(metaprograms_scores$snc_epi_vitro$JHU006_UPPER_AERODIGESTIVE_TRACT), "vivo"=unlist(metaprograms_scores$HNSCC.Epidif.1$JHU006_UPPER_AERODIGESTIVE_TRACT))

p4 <- ggplot(epi_snc_plot[epi_snc_plot$vitro<3,], aes(x=vitro, y=vivo)) +
  geom_point(alpha=0.4, color="gray80", size=4) +
  theme(axis.line=element_blank(), plot.title=element_text(hjust = 0)) +
  panel_border(colour = "black", 0.4) +
  labs(x="Epithelial senescence\nprogram score", y = "HNSCC epithelial diff. I\nprogram score", title= "JHU006 - Tongue squamous\ncell carcinoma" )

stress_plot <- data.frame("vitro"=unlist(metaprograms_scores$stress_vitro$JHOS2_OVARY), "vivo"= unlist(metaprograms_scores$stress$JHOS2_OVARY))

p5 <- ggplot(stress_plot, aes(x=vitro, y=vivo)) +
  geom_point(alpha=0.4, color="gray80", size=4) +
  theme(axis.line=element_blank(),plot.title=element_text(hjust = 0)) +
  panel_border(colour = "black", 0.4) +
  labs(x="Stress response\nprogram score", y = "Melanoma/HNSCC stress\nprogram score", title =  "JHOS2 - Ovarian cancer") 

egg::ggarrange(p1, p2, p3, p4, p5, ncol=5)

pdf("nc_vitro_vs_vivo_scatter.pdf", width = 18, height = 4, onefile = F)
