# HaploVar

## What is HaploVar

`HaploVar` defines local haplotypes and is designed to be used as part of a GWAS or genomic selection pipeline. `HaploVar` takes a VCF file and a LD matrix, calculates local haplotypes, identifies haplotype variants and formats the output to be compatible with a wide range of GWAS and genomic selection tools. Halotypes improve the ability of GWAS to detect QTLs and increase the trait prediction accuracy of genomic selection studies. 

## Installation

`HaploVar` can be installed using the following code:

``` r
#install.packages("devtools")
devtools::install_github("TessaMacNish/HaploVar")
library(HaploVar)
```
