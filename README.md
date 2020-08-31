Pan-cancer single cell RNA-seq uncovers recurring programs of cellular heterogeneity
---------------------------------------------------------------------------------------

This resource provides the R code to reproduce key results described in Kinker et al. "Pan-cancer single cell RNA-seq uncovers recurring programs of cellular heterogeneity".

The analyses are divided in 6 main modules:  
**1.** Identifying discrete and continuous patterns of expression heterogeneity within cell lines and human tumors.  
**2.** Definying heterogeneity patterns that are shared between multiple cell lines and between multiple human tumors (i.e. recurrent heterogeneous programs, RHPs).  
**3-4.** Comparying RHPs found in cell lines to RHPs found in human tumor samples.   
**5-6.** Evaluating the association between expression and genetic heterogeneity in cell lines. 

### Getting started
**1.** Clone Github repository 
```
git clone https://github.com/gabrielakinker/CCLE_heterogeneity.git
```

**2.** Set the working directory to ``CCLE_heterogeneity`` 

**3.** Download the data provided in the [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP542/pan-cancer-cell-line-heterogeneity) [CCLE_scRNAseq_github.zip](https://singlecell.broadinstitute.org/single_cell/data/public/SCP542/pan-cancer-cell-line-heterogeneity?filename=CCLE_scRNAseq_github.zip) into the ``CCLE_heterogeneity`` GitHub directory. 

**4.** Extract files 
```
tar -xvf CCLE_scRNAseq_github.tar.gz --strip-components 0
```

**5.** Install required R packages
```
Rscript packages.R
```
**5.** Run one of the 6 code modules in R
* module1_expr_heteroge.R         
* module2_rhp.R                   
* module3_vitro_vs_vivoI.R   
* module4_vitro_vs_vivoII.R 
* module5_cna_subclones.R
* module6_subclones_vs_expr_heteroge.R

**6.** Output files will be saved in the ``Output`` directory. 

### General notes

* Please, note that due to the stochastic nature of methods such as T-distributed Stochastic Neighbor Embedding (t-SNE) and nonnegative matrix factorization (NMF) outputs from module1_expr_heterogeneity.r might slightly differ dependending on the version of R/R packages used.

* As all outputs are already provided in the ``Expected_results`` directory, it is possible to run each code module independently.

* Running module1 takes several hours. 

### **Requirements**

* R (tested in version 3.6.3).


