#' Collate haplotype_variants Tables
#'
#' This function collates a list of output files from haplotype_variants.
#' @param haplotype_variants_list A list of the tables created by the define_haplotypes function.
#' @param format The format you want the output table to be in. This should be the same number you used when running define_haplotypes.
#' @return A collated table of haplotype variants.
#' @examples
#' collate_haplotype_variants(list, format = 1)
#' @export
#' @importFrom magrittr %>%
collate_haplotype_variants <- function(haplotype_variants_list, format = 1) {
  # Initialize an empty data frame
  all_variants <- data.frame()

  # Loop through each table in the list
  for (variants in haplotype_variants_list) {
    if (length(variants) > 0) {
      if (format == 3 | format == 5) {
        tvariants <- t(variants)
        all_variants <- rbind.data.frame(all_variants, tvariants)
      } else {
        # Add these variants to the overall table
        all_variants <- rbind.data.frame(all_variants, variants)
      }
    }
  }
  #Format output
  if (format == 3 | format == 5) {
    output <- t(all_variants)
  } else {
    output <- all_variants
  }

  # Return the combined list of haplotype tables
  return(output)
}
