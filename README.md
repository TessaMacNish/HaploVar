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

For a tutorial on how to use `HaploVar` please refer to [introduction](https://htmlpreview.github.io/?https://github.com/TessaMacNish/HaploVar/blob/main/vignettes/introduction.html) This tutorial describes how and when to use each function. Each parameter is explained, and example code is provided. 

## Pipeline Overview 

An example of a `HaploVar` workflow for preparing for genomic prediction is below:

1) Collect your genotype data. This should be in VCF format and include SNP genotypes.
2) Calculate linkage disequilibrium for your genotype data using a tool such as PLINK (Purcell et al. 2007, https://zzz.bwh.harvard.edu/plink/ld.shtml) or the ld function of snpStats (Clayton, 2025). For example:

``` bash
plink --allow-extra-chr --r2 square --vcf vcf --out output.ld
```

3) Read your VCF and linkage disequilibrium matrix into R. `crosshap` has an inbuilt function for this purpose (Marsh et al. 2023).

``` r
install.packages("crosshap")
library(crosshap)
VCF <- read_vcf("vcf")
LD <- read_LD("output.ld",vcf = VCF) 
```

4) Calculate haplotypes. The following code outputs haplotype tables displaying which SNPs group together to form haplotypes.

``` r
haplotype_list <- define_haplotypes(VCF, LD, epsilon = 0.8) 
```

5) Identify haplotype variants and format the output for genomic selection pipelines.

``` r
format3 <- haplotype_variants(VCF, LD, epsilon = 0.8, format = 3)
```
6) Save the output

``` r
write.csv(format3, "output.csv", row.names=TRUE, quote = FALSE )
```

7) Run your preferred genomic selection pipeline


NOTE: If you want to prepare haplotype variants for GWAS instead of genomic selection change format to equal 6 and save the data with the following code:

``` r
write.table(format6, "output.vcf", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
```

### References 
Clayton D. snpStats: SnpMatrix and XSnpMatrix classes and methods. R pack-age version 1.59.2. 2025. https://bioconductor.org/packages/snpStats

Marsh JI, Petereit J, Johnston BA et al. crosshap: R package for local haplotype visualization for trait association analysis. Bioinformatics 2023;39. https://doi.org/10.1093/bioinformatics/btad518 

Purcell S, Neale B, Todd-Brown K et al. PLINK: A tool set for whole-genome association and population-based linkage analyses. Am J Hum Genet 2007;81:559–575. https://doi.org/10.1086/519795
