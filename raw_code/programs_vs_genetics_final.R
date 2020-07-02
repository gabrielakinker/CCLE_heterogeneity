############################################################################################################################################################################################### PROGRAMS OF VARIABILITY VS GENETICS ########################################################

### reading expression data
expr <- readRDS("~/CCLE/data processing/scCCLE&pool1_tpm_by_cell_line.rds") 

### reading gmm test - cnv clones
gmm_final <- readRDS("~/CCLE/cnv clones/gmm_test_filtered.RDS")

### defining clones
clones <- list()
for(i in names(gmm_final)) {
  if(length(gmm_final[[i]]) > 1) {
    a <- colnames(expr[[i]])
    b <- do.call(paste, lapply(gmm_final[[i]], function(x) x[["class"]]))
    b <- data.frame("cell"=a, "class"=b, stringsAsFactors = F)
    c <- sapply(gmm_final[[i]], function(x) x[["uncer"]])
    b <- b[apply(c,1,function(x) length(which(x<=0.1))==length(x)),]
    b <- b[is.element(b$class, names(which(table(b$class)>=5))),]
    d <- data.frame(names(sort(table(b$class), decreasing = T)), 1:length(unique(b$class)), stringsAsFactors = F)
    b$final_class <- d[,2][match(b$class, d[,1])]
    b$final_class <- paste(i, "_", b$final_class, sep="")
    clones[[i]] <- b
  }
  if(length(gmm_final[[i]]) == 1) {
    a <- colnames(expr[[i]])
    b <- data.frame("cell"=a, "class"=gmm_final[[i]][[1]][["class"]], stringsAsFactors = F)
    b <- b[gmm_final[[i]][[1]][["uncer"]] <= 0.1,]
    b$final_class <- paste(i, "_", b$class, sep="")
    clones[[i]] <- b
  }
}


########################## NMF PROGRAMS as continuous
### reading cell scores per cell line - programs nmf
h_coef_6 <- readRDS("~/CCLE/nmf/h_coef_nmf6.rds")
for(i in names(h_coef_6)) {
  colnames(h_coef_6[[i]]) <- paste(i, "_6.", 1:6, sep = "")
}

h_coef_7 <- readRDS("~/CCLE/nmf/h_coef_nmf7.rds")
for(i in names(h_coef_7)) {
  colnames(h_coef_7[[i]]) <- paste(i, "_7.", 1:7, sep = "")
}

h_coef_8 <- readRDS("~/CCLE/nmf/h_coef_nmf8.rds")
for(i in names(h_coef_8)) {
  colnames(h_coef_8[[i]]) <- paste(i, "_8.", 1:8, sep = "")
}

h_coef_9 <- readRDS("~/CCLE/nmf/h_coef_nmf9.rds")
for(i in names(h_coef_9)) {
  colnames(h_coef_9[[i]]) <- paste(i, "_9.", 1:9, sep = "")
}

### combining results
cell_scores <-  unlist(unname(lapply(c(h_coef_6, h_coef_7, h_coef_8, h_coef_9), function(x) apply(x, 2, as.list))), recursive = F)
cell_scores <- lapply(cell_scores, unlist)

### applying redundancy filter
redundancy_filter <- readRDS("~/CCLE/nmf/redundancy_filter_nmf.RDS")
cell_scores <- cell_scores[redundancy_filter]

### removing discrete programs
nmf_discr <- readRDS("~/CCLE/nmf/continous_vs_discr.RDS")
cell_scores <- cell_scores[!is.element(names(cell_scores), names(nmf_discr))]

### association between continuous programs and cnv clones
conti_gmm <- list()

for(i in names(cell_scores)[is.element(gsub(".{4}$", "", names(cell_scores)), names(clones))]) {
  a <- cell_scores[[i]] # gets nmf scores for program i  
  b <- clones[[gsub(".{4}$", "", i)]] # gets clones
  a <- a[b$cell]  
  if(length(unique(b$final_class))==2) {
   conti_gmm[[i]] <- t.test(a~b[,"final_class"])$p.value  
  }
  if(length(unique(b$final_class))>2) {
   conti_gmm[[i]] <- summary(aov(a~b[,"final_class"]))[[1]][1,5]     
  }
}

conti_gmm_final <- conti_gmm[sapply(conti_gmm, function(x) x<0.001)]

100*length(conti_gmm_final)/length(cell_scores)

########################## DISCRETE PROGRAMS eps 1.8
### reading final discrete programs and their signature genes
discr_programs <- readRDS("~/CCLE/discrete programs/discr_clusters_minpt5_eps1.8.RDS") # cluster cell composition
discr_programs <- lapply(discr_programs, function(x) x[["clusters_cells"]])
discr_programs_unlist <- unlist(discr_programs, recursive = F)

### association between discrete programs and cnv clones                         
discr_gmm <- list()

for(i in names(discr_programs_unlist)[is.element(gsub("\\..*", "", names(discr_programs_unlist)), names(clones))]) {
  a <- colnames(expr[[sub("\\..*", "", i)]]) # gets all cells from the selected cell line
  b <- as.character(is.element(a,discr_programs_unlist[[i]])) # checks which cells are part of the discrete cluster
  c <- clones[[sub("\\..*", "", i)]] # gets gmm test for the selected cell line for each chromosome arm
  e <- data.frame("clone" = c[,"final_class"], "discr" = b[match(c[,"cell"], a)]) # combines gmm test results and discrete program classification
  discr_gmm[[i]] <- chisq.test(table(e))$p.value # fisher test - classification vs discrete program classification
}
                         
discr_gmm_final <- discr_gmm[sapply(discr_gmm, function(x) x<0.001)]

100*length(discr_gmm_final)/length(discr_programs_unlist)                         
                         
                                              
### plotting 
plot_cnv <- data.frame("Discrete\nprograms" = c(length(discr_programs_unlist)-length(discr_gmm_final), length(discr_gmm_final)), "Continuous\nprograms" = c(length(cell_scores)-length(conti_gmm_final), length(conti_gmm_final)), check.names = F)
rownames(plot_cnv) <- c("Not linked to CNA sub-clone", "Linked to CNA sub-clone")
plot_cnv <- round(t(100*t(plot_cnv)/colSums(plot_cnv)),0)
plot_cnv <- melt(as.matrix(plot_cnv))

ggplot(plot_cnv, aes(x="", y=value, fill=Var1)) +
  geom_bar(stat="identity", alpha=0.8) +
  coord_polar("y") +
  scale_fill_manual(values=c("gray88", "maroon"), name="") +
  labs(x="", y="") +
  facet_wrap(~ Var2)  +
  theme(panel.background=element_blank() ,axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), strip.background = element_blank(), legend.position = "bottom", legend.justification = "center", strip.text.x=element_text(margin=margin(b=5)))

pdf("percentage_programs_linked_cnvII.pdf", width = 6, height = 5)

