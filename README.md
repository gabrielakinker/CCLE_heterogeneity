Pan-cancer single cell RNA-seq uncovers recurring programs of cellular heterogeneity
---------------------------------------------------------------------------------------

This resource provides the R code to reproduce key results described in Kinker et al. "Pan-cancer single cell RNA-seq uncovers recurring programs of cellular heterogeneity".

The analyses are divided in 4 main modules:  
**1.** Identifying discrete and continuous patterns of expression heterogeneity within cell lines and human tumors.  
**2.** Definying heterogeneity patterns that are shared between multiple cell lines and between multiple human tumors (i.e. recurrent heterogeneous programs, RHPs).  
**3.** Comparying RHPs found in cell lines to heterogeneity programs observed in human tumor samples.   
**4.** Evaluating the association between expression and genetic heterogeneity in cell lines. 

## Getting Started
**1.** Clone Github repository 
```
git clone https://github.com/gabrielakinker/CCLE_heterogeneity.git
```
**2.** Download the data provided in the [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP542/pan-cancer-cell-line-heterogeneity) (CCLE_heterogeneity_Rfiles.zip) into the ``CCLE_heterogeneity`` GitHub directory. 

**3.** Install required R packages
```
Rscript packages.R
```
**4.** Set the working directory to ``CCLE_heterogeneity`` 

**5.** Run one of the 4 main code modules
* module1_expr_heterogeneity.R
* module2_rhp.R
* module3_invivo_comp.R
* module4_gen_heterogeneity.R 

**6.** Output files will be saved in the ``Output`` directory. 

## General notes

* Please, note that due to the stochastic nature of methods such as T-distributed Stochastic Neighbor Embedding (t-SNE) and nonnegative matrix factorization (NMF) outputs from module1_expr_heterogeneity.r might slightly differ dependending on the version of R/R packages used.

* As all outputs are already provided in the ``Expected_results`` directory, it is possible to run each code module independently.

* Running module1_expr_heterogeneity.r takes several hours. 

## **Requirements**

* R (tested in version 3.6.3).


