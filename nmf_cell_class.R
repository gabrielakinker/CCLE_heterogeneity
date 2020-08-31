# -------------------------------------------------------------------------------------------
# Function for getting heterogeneity programs using nonnegative matrix factorization (NMF)
# ------------------------------------------------------------------------------------------- 

# - nmf_programs_cells = output of nmf_programs (h_coef)
# - rank = NMF factorization ranks used in  nmf_programs 

# Returns a list with cells assigned to each NMF program 

nmf_cell_class <- function(nmf_programs_cells, ranks=6:9) {
    a <- list()
    for(i in paste0("_", ranks)) {
      b <- nmf_programs_cells[,grep(i, colnames(nmf_programs_cells))]
      c <- split(rownames(b), apply(b, 1, function(x) colnames(b)[which.max(x)]))
      a <- c(a, c)                              
  }                                   
  return(a)
}
    
