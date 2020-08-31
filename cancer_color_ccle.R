# set a color for each cancer type
library(RColorBrewer)
cancer_color <- data.frame("type"=sort(unique(readRDS("CCLE_heterogeneity_Rfiles/CCLE_metadata.RDS") $cancer_type_trunc)), "color"= c(brewer.pal(12, "Set3")[c(1:6,8,7,10:11)], "maroon","gray93", "yellow2", "goldenrod1", "slateblue2"),stringsAsFactors = F)