# HaploVar

## What is HaploVar

`HaploVar` defines local haplotypes and is designed to be used as part of a GWAS or genomic selection pipeline. `HaploVar` takes a VCF file and a LD matrix, calculates local haplotypes, identifies haplotype variants and formats the output to be compatible with a wide range of GWAS and genomic selection tools. Halotypes improve the ability of GWAS to detect QTLs and increase the trait prediction accuracy of genomic selection studies. 

## Installation

`HaploVar` can be installed using the following code:

``` r
#Install the development version from github
install.packages("devtools")
devtools::install_github("TessaMacNish/HaploVar")

#or install directly from CRAN
install.packages("HaploVar")

library(HaploVar)

```
## Documentaion

For a basic tutorial on how to use `HaploVar` please refer to [introduction](https://htmlpreview.github.io/?https://github.com/TessaMacNish/HaploVar/blob/main/vignettes/introduction.html)
