# ---------------------------------------------------------------------------------------------
# Function for getting a control geneset
# --------------------------------------------------------------------------------------------- 

# - ave_tpm = named vector with the average tpm/cpm expression of genes in the expression matrix
# - program = a character vector of genes in the program 
# - bins = number of bins used to divide ranked genes
# - size = control gene set size scaling factor (e.g. 100x the size of the original program)

# Returns a characther vector with the control gene set 

control_geneset <- function(ave_tpm, program, bins=50, size=100, seed=1) {
  
  if(length(which(!is.element(program, names(ave_tpm)))) != 0) stop("Check the genes in the program. Not all of them are found in the expression matrix.")
  
  # sorting
  ave_tpm_sorted <- sort(ave_tpm) 
  
  # define 50 bins of gene expression
  expr_bins <- list()
  d <- seq(from = 0, to = length(ave_tpm_sorted), length.out = bins+1)

  for(i in 1:bins) {
    expr_bins[[i]] <- names(ave_tpm_sorted)[(1+d[i]):d[i+1]] 
  }
  
  # for each program gene, selects 100 control genes from the correspondent bin 
  control_geneset <- NULL
  
  for(i in 1:length(program)) {
    a <- program[i] # select a program gene
    b <- which(unlist(lapply(expr_bins, function(x) is.element(a, x)))) # indentifies the bin of the selected gene
    if(length(b) != 1) stop("Please, select a different number of bins.")
    c <- expr_bins[[b]] # select the desired bin
    c <- c[c!=a] # prevent the program gene from making into the control geneset
    set.seed(seed)
    d <- sample(c, size, replace = T) # sample 100x from the selected bin
   
    control_geneset <- c(control_geneset,d)
  }
  rm(.Random.seed, envir=.GlobalEnv)
  return(control_geneset)  
  
}




