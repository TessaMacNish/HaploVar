#' Define Haplotypes Globally
#'
#'
#' This function requires a list of VCF files and an LD matrices. The list of VCF files and LD matrices must be the same length.
#' It will then define local haplotypes for each pair of files (VCF and LD matrix) and return a list of tables.
#' Each table within the list represents one haplotype.
#' These haplotype tables display the SNP genotypes within the haplotype.
#' @param vcf_list A list of VCF files.
#' @param LD_list A list LD matrix files.
#' @param epsilon A list of epsilon values the same length as the list of VCF files. The epsilon affects haplotype size. It is a parameter of the DBSCAN clustering tool. The default is 0.6.
#' @param MGmin The minimum number of SNPs within a cluster for it to be defined as a haplotype.The default is 30.
#' @param hetmiss_as Affects how missing data is handled for all instances where one allele in a genotype is missing.If hetmiss_as = "allele" the genotype is assumed to be heterozygous. If hetmiss_as = "miss" the genotype is treated as NA.
#' @param keep_outliers If FALSE, removes SNPs that are determined to be outliers.
#' @return A collated list of haplotype tables for all VCF files provided.
#' @export
#' @importFrom magrittr %>%
define_haplotypes_globally <- function(vcf_list, LD_list, epsilon = NULL, MGmin = 30, hetmiss_as = "allele", keep_outliers = FALSE) {
  # If epsilon is NULL, set it to a vector of 0.6 with the same length as vcf_list
  if (base::is.null(epsilon)) {
    epsilon <- base::rep(0.6, base::length(vcf_list))
  }

  # Ensure the lists of VCF, LD files, and epsilon values are of the same length
  if (base::length(vcf_list) != base::length(LD_list) || base::length(vcf_list) != base::length(epsilon)) {
    base::stop("The number of VCF files, LD files, and epsilon values must be the same.")
  }

  # Initialize an empty list to store all haplotypes
  all_haplotypes <- base::list()

  # Use sapply to apply the define_haplotypes function to each VCF and LD pair with respective epsilon
  base::sapply(base::seq_along(vcf_list), function(i) {
    vcf <- vcf_list[[i]]
    LD <- LD_list[[i]]
    eps <- epsilon[i]  # Use the corresponding epsilon value

    # Get haplotypes for the current chromosome with the specific epsilon
    haplotypes <- define_haplotypes(vcf, LD, epsilon = eps, MGmin = MGmin, hetmiss_as = hetmiss_as, keep_outliers = keep_outliers)

    if (base::length(haplotypes) > 0) {
      # Add these haplotypes to the overall list
      all_haplotypes <<- c(all_haplotypes, haplotypes)
    }

  })

  # Return the combined list of haplotype tables
  return(all_haplotypes)
}
