#' Identify Haplotype Variants Globally
#'
#' This function requires a list of VCF files and an LD matrices.
#' It will then define local haplotypes and identify the variants for each haplotype.
#' The output can be formatted in six ways, to be compatible with a wide range of GWAS and genomic selection tools.
#' @param vcf_list A list of VCF files.
#' @param LD_list A list of LD matrix files.
#' @param epsilon A list of epsilon values the same length as the list of VCF files. The epsilon affects haplotype size. It is a parameter of the DBSCAN clustering tool. The default is 0.6.
#' @param MGmin The minimum number of SNPs within a cluster for it to be defined as a haplotype.The default is 30.
#' @param minFreq The minimum number of individuals a haplotype variant must be present in to be considered a valid haplotype variant. The default is 2.
#' @param hetmiss_as Affects how missing data is handled for all instances where one allele in a genotype is missing.If hetmiss_as = "allele" the genotype is assumed to be heterozygous. If hetmiss_as = "miss" the genotype is treated as NA.
#' @param keep_outliers If FALSE removes SNPs, that are determined to be outliers.
#' @param format The output format. There are six different output formats (1,2,3,4,5,6).
#' @return A table of haplotype genotypes in your chosen format.
#' @export
#' @importFrom magrittr %>%
haplotype_variants_global <- function(vcf_list, LD_list, epsilon = NULL, MGmin = 30, minFreq = 2, hetmiss_as = "allele", keep_outliers = FALSE, format = 1) {
  # If epsilon is NULL, set it to a vector of 0.6 with the same length as vcf_list
  if (base::is.null(epsilon)) {
    epsilon <- base::rep(0.6, base::length(vcf_list))
  }

  # Ensure the lists of VCF, LD files, and epsilon values are of the same length
  if (base::length(vcf_list) != base::length(LD_list) || base::length(vcf_list) != base::length(epsilon)) {
    base::stop("The number of VCF files, LD files, and epsilon values must be the same.")
  }

  # Initialize an empty data frame to store all variants
  all_variants <- base::data.frame()

  # Apply the haplotype_variants function to each VCF and LD pair with the respective epsilon
  base::sapply(seq_along(vcf_list), function(i) {
    vcf <- vcf_list[[i]]
    LD <- LD_list[[i]]
    eps <- epsilon[i]  # Use the corresponding epsilon value

    # Get variants for the current chromosome with the specific epsilon
    variants <- haplotype_variants(vcf, LD, epsilon = eps, MGmin = MGmin, minFreq = minFreq, hetmiss_as = hetmiss_as, keep_outliers = keep_outliers, format = format)

    if (base::length(variants) > 0) {
      if (format == 3 | format == 5) {
        tvariants <- base::t(variants)
        all_variants <<- base::rbind.data.frame(all_variants, tvariants)
      } else {
        # Add these variants to the overall list
        all_variants <<- base::rbind.data.frame(all_variants, variants)
      }
    }
    return(all_variants)
  })

  # Format output based on the format parameter
  if (format == 3 | format == 5) {
    output <- base::t(all_variants)
  } else {
    output <- all_variants
  }

  # Return the combined list of haplotype tables
  return(output)
}
