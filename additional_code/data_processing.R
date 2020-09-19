# ----------------------------------------------------------------------------------------------------
# Processing and quality control of UMI counts 
# ----------------------------------------------------------------------------------------------------

# read data UMI counts
expr <- read.table("UMIcount_data.txt", sep="\t", header=T, row.names=1, check.names=F)

# extract annotation
annotation <- expr[1:2,]
expr <- expr[-(1:2),]

# filter cells based on the number of detected genes (i.e. complexity)
complexity <- c()   

for(i in 1:ncol(umi))  {
     complexity <- c(complexity, length(which(expr[,i]!=0)))
}                   

complexity_filter <- complexity > 2000 & complexity < 9000
expr <- expr[,complexity_filter]
annotation <- annotation[,complexity_filter]                    

# filter cell lines based on the number of cells profiled
cells_per_celline <- table(t(annotation[1,]))
cell_line_filter <- is.element(t(annotation[1,]), names(cells_per_celline[cells_per_celline>50]))

expr <- expr[,cell_line_filter]
annotation <- annotation[,cell_line_filter]

# transform UMI counta into CPM 
expr <- t(t(expr)/colSums(expr))*1000000
