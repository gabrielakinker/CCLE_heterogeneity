# -------------------------------------------------------------------------------------------
# Function for getting heterogeneity programs using nonnegative matrix factorization (NMF)
# ------------------------------------------------------------------------------------------- 

# - cpm = CPM expression matrix (rows = genes, columns = cells)
# - is.log = indicates if the data is log transformed
# - rank = NMF factorization rank 
# - method = NMF algorithm

# Returns a list with NMF program scores for genes and cells

library(NMF)

nmf_programs <- function(cpm, is.log=F, rank, method="snmf/r", seed=1) {
    
  if(is.log==F) CP100K_log <- log2((cpm/10) + 1) else CP100K_log <- cpm
  CP100K_log <- CP100K_log[apply(CP100K_log, 1, function(x) length(which(x > 3.5)) > ncol(CP100K_log)*0.02),]
  CP100K_log <- CP100K_log - rowMeans(CP100K_log)
  CP100K_log[CP100K_log < 0] <- 0
 
  nmf_programs <- nmf(CP100K_log, rank=rank, method=method, seed=seed)
  
  nmf_programs_scores <- list(w_basis=basis(nmf_programs), h_coef=t(coef(nmf_programs)))
                                 
  return(nmf_programs_scores)
}
    
