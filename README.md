Pan-cancer single-cell RNA-seq identifies recurring programs of cellular heterogeneity
---------------------------------------------------------------------------------------

This resource provides the R code to reproduce key results described in Kinker et al., 2020, "Pan-cancer single-cell RNA-seq identifies recurring programs of cellular heterogeneity".

The analyses are divided into 6 main modules:  
**1.** Identifying discrete and continuous patterns of expression heterogeneity within cell lines and human tumors.  
**2.** Definying heterogeneity patterns that are shared between multiple cell lines and between multiple human tumors (i.e. recurrent heterogeneous programs, RHPs).  
**3-4.** Comparying RHPs found in cell lines to RHPs found in human tumor samples.   
**5-6.** Evaluating the association between expression and genetic heterogeneity in cell lines. 

For questions, please contact gabriela.kinker@gmail.com

### Getting started
**1.** Clone Github repository. 
```
git clone https://github.com/gabrielakinker/CCLE_heterogeneity.git
```

**2.** Set the working directory to ``CCLE_heterogeneity``. 

**3.** Download the data provided ([CCLE_scRNAseq_github.zip](https://singlecell.broadinstitute.org/single_cell/study/SCP542/pan-cancer-cell-line-heterogeneity)), which includes processed, post-QC scRNA-seq data, into the ``CCLE_heterogeneity`` GitHub directory. 

**4.** Unzip files. 
```
unzip CCLE_scRNAseq_github.zip && mv CCLE_scRNAseq_github/* . && rm -d CCLE_scRNAseq_github
```

**5.** Install required R packages.
```
Rscript packages.R
```
**5.** Run one of the 6 code modules in R.
* module1_expr_heteroge.R         
* module2_rhp.R                   
* module3_vitro_vs_vivoI.R   
* module4_vitro_vs_vivoII.R 
* module5_cna_subclones.R
* module6_subclones_vs_expr_heteroge.R

**6.** Output files will be saved in the ``Output`` directory. 

### General notes

* Please, note that due to the stochastic nature of methods such as T-distributed Stochastic Neighbor Embedding (t-SNE) and nonnegative matrix factorization (NMF) outputs from module1_expr_heteroge.r might slightly differ dependending on the version of R/R packages used.

* As all outputs are already provided in the ``Expected_results`` directory ([CCLE_scRNAseq_github.zip](https://singlecell.broadinstitute.org/single_cell/study/SCP542/pan-cancer-cell-line-heterogeneity)), it is possible to run each code module independently.

* Running module1 takes several hours. 

* Code used for quality control and processing of UMI counts ([UMIcount_data.txt](https://singlecell.broadinstitute.org/single_cell/study/SCP542/pan-cancer-cell-line-heterogeneity)) can be found in the ``additional_code`` directory. 

### **Requirements**

* R (tested in version 3.6.3).


