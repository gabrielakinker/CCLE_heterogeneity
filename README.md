Pan-cancer single cell RNA-seq uncovers recurring programs of cellular heterogeneity
---------------------------------------------------------------------------------------

This resource provides the R code to reproduce key results described in Kinker et al. "Pan-cancer single cell RNA-seq uncovers recurring programs of cellular heterogeneity".

The analyses are divided in 4 main modules:  
**1.** Identifying discrete and continuous patterns of expression heterogeneity within cell lines.  
**2.** Definying heterogeneity patterns that are shared between multiple cell lines (i.e. recurrent heterogeneous programs, RHPs).  
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
Rscript packages.r
```
**4.** Set the working directory to ``CCLE_heterogeneity`` 

**5.** Run one of the 4 main code modules (i.e. module1_expr_heterogeneity.r, module2_rhp.r, module3_invivo_comp.r, module4_gen_heterogeneity.r). Output files will be saved in the ``Output`` directory. 

## General notes

## **Requirements**

* R (tested in version 3.6.3).


