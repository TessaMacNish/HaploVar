## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----installation, eval = FALSE-----------------------------------------------
# install.packages("devtools")
# devtools::install_github("TessaMacNish/HaploVar")

## ----message = FALSE, warning = FALSE, eval = TRUE----------------------------
library(HaploVar)

## ----load-data, eval = TRUE, include = FALSE----------------------------------
# Load data
data("vcf")
data("LD")

## ----load-vcf, eval = TRUE, include = TRUE------------------------------------
head(vcf, c(5,10))

## ----load-LD, eval = TRUE, include = TRUE-------------------------------------
head(LD, c(5,5))

## ----define_haplotypes, message = F, results = 'hide', cache = T--------------
##Run define_haplotypes
haplotype_list <- define_haplotypes(vcf, LD, epsilon = 0.8) #this produces a list of haplotype tables
haplotype_table <- haplotype_list[[1]] #this is the first haplotype tables

## ----head-haplotype-table, eval = TRUE, include = TRUE------------------------
head(haplotype_table, c(5,5))

## ----collate_haplotype_list, message = F, results = 'hide', cache = T---------
##Prepare the data
haplotype_list2 <- haplotype_list #We are copying the haplotype_list so that we have two lists for the demonstration 
list_outputs <- list(haplotype_list, haplotype_list2) #The input must be in list format
##Run collate_haplotype_list
collate_haplotype_list <- collate_define_haplotypes(list_outputs)

## ----define_haplotypes_globally, message = F, results = 'hide', cache = T-----
##Prepare the data
vcf2 <- vcf 
vcf_list <- list(vcf, vcf2) #The vcf files must be in list format
LD2 <- LD
LD_list <- list(LD, LD2) #The LD matrices must be in list format

##Prepare parameters
epsilon_list <- c(0.8, 0.8) #The length of this list must be the same as the number of vcf files.

## ----define_haplotypes_globally_2, message = F, results = 'hide', cache = T----
##Run define_haplotypes_globally
haplotype_list_global <- define_haplotypes_globally(vcf_list, LD_list, epsilon = epsilon_list)

## ----haplotype_variants, message = F, results = 'hide', cache = T-------------
##Run haplotype_variants
format1 <- haplotype_variants(vcf, LD, epsilon = 0.8, format = 1)
format2 <- haplotype_variants(vcf, LD, epsilon = 0.8, format = 2)
format3 <- haplotype_variants(vcf, LD, epsilon = 0.8, format = 3)
format4 <- haplotype_variants(vcf, LD, epsilon = 0.8, format = 4)
format5 <- haplotype_variants(vcf, LD, epsilon = 0.8, format = 5)
format6 <- haplotype_variants(vcf, LD, epsilon = 0.8, format = 6)

## ----format1, eval = TRUE, include = TRUE-------------------------------------
head(format1, c(5,6))

## ----format2, eval = TRUE, include = TRUE-------------------------------------
head(format2, c(5,6))

## ----format3, eval = TRUE, include = TRUE-------------------------------------
head(format3, c(5,3))

## ----format4, eval = TRUE, include = TRUE-------------------------------------
head(format4, c(5,6))

## ----format5, eval = TRUE, include = TRUE-------------------------------------
head(format5, c(5,3))

## ----format6, eval = TRUE, include = TRUE-------------------------------------
head(format6, c(5,10))

## ----collate_haplotype_variants, message = F, results = 'hide', cache = T-----
##format1
format1B <- format1
format1_list <- list(format1, format1B) #The input for collate_haplotype_variants must be in list format.
format1_collate <- collate_haplotype_variants(format1_list, format = 1)
##Format2
format2B <- format2
format2_list <- list(format2, format2B)
format2_collate <- collate_haplotype_variants(format2_list, format = 2)
##format3
format3B <- format3
format3_list <- list(format3, format3B)
format3_collate <- collate_haplotype_variants(format3_list, format = 3)
##Format4
format4B <- format4
format4_list <- list(format4, format4B)
format4_collate <- collate_haplotype_variants(format4_list, format = 4)
##format5
format5B <- format5
format5_list <- list(format5, format5B)
format5_collate <- collate_haplotype_variants(format5_list, format = 5)
##Format6
format6B <- format6
format6_list <- list(format6, format6B)
format6_collate <- collate_haplotype_variants(format6_list, format = 6)

## ----haplotype_variants_global, message = F, results = 'hide', cache = T------
format1_global <- haplotype_variants_global(vcf_list, LD_list, epsilon = epsilon_list, format = 1)
format2_global <- haplotype_variants_global(vcf_list, LD_list, epsilon = epsilon_list, format = 2)
format3_global <- haplotype_variants_global(vcf_list, LD_list, epsilon = epsilon_list, format = 3)
format4_global <- haplotype_variants_global(vcf_list, LD_list, epsilon = epsilon_list, format = 4)
format5_global <- haplotype_variants_global(vcf_list, LD_list, epsilon = epsilon_list, format = 5)
format6_global <- haplotype_variants_global(vcf_list, LD_list, epsilon = epsilon_list, format = 6)

