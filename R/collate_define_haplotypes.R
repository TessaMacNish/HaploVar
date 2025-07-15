#' Collate define_haplotypes Lists
#'
#' This function collates a list of output files from define_haplotypes.
#' @param haplotype_list A list of the lists created by the define_haplotypes function.
#' @return A collated list of all haplotype tables.
#' @export
#' @importFrom magrittr %>%
collate_define_haplotypes <- function(haplotype_list) {
  # Initialize an empty list to store all haplotype tables
  all_haplotypes <- base::list()

  # Loop through each list of haplotypes lists
  for (lists in haplotype_list) {
    if (base::length(lists) > 0) {
      # Add these haplotypes to the overall list
      all_haplotypes <- c(all_haplotypes, lists)
    }
  }

  # Return the combined list of haplotype tables
  return(all_haplotypes)
}
