utils::globalVariables(c("cluster", "rowname", "hap", "row_index", "pair"))

#' Define Haplotypes
#'
#' This function requires a VCF and an LD matrix.
#' It will then define local haplotypes and return a list of tables.
#' Each table within the list represents one haplotype.
#' These haplotype tables display the SNP genotypes within the haplotype.
#' @param vcf A VCF file.
#' @param LD A LD matrix file.
#' @param epsilon Affects haplotype size. It is a parameter of the DBSCAN clustering tool. The default is 0.6.
#' @param MGmin The minimum number of SNPs within a cluster for it to be defined as a haplotype.The default is 30.
#' @param hetmiss_as Affects how missing data is handled for all instances where one allele in a genotype is missing.If hetmiss_as = "allele" the genotype is assumed to be heterozygous. If hetmiss_as = "miss" the genotype is treated as NA.
#' @param keep_outliers If FALSE, removes SNPs that are determined to be outliers.
#' @return A list of haplotype tables.
#' @export
#' @importFrom magrittr %>%
#' @importFrom stats median sd
define_haplotypes <- function(vcf, LD, epsilon = 0.6, MGmin = 30,  hetmiss_as = "allele", keep_outliers = FALSE) {
  ##format the LD file and vcf into usable formats
  bin_vcf <- dplyr::select(vcf, -c(1, 2, 4:9)) %>%
    tibble::column_to_rownames("ID") %>%
    dplyr::mutate_all(function(x) {
      base::switch(hetmiss_as,
                   allele = base::ifelse(x == "1|." , "1|0",
                                         base::ifelse(x == "1/.", "1/0",
                                                      base::ifelse(x == ".|1" , "0|1",
                                                                   base::ifelse(x == "./1", "0/1",
                                                                                base::ifelse(x == "0|." | x == ".|0" , "0|0",
                                                                                             base::ifelse( x == "0/." | x == "./0", "0/0", x)))))),
                   miss = base::ifelse(x == "1|." | x == ".|1" | x == "1/." | x == "./1" |
                                         x == "0|." | x == ".|0" | x == "0/." | x == "./0", NA, x)
      )
    })
  ##Define haplotypes using DBSCAN
  db_out <- dbscan::dbscan(LD, eps = epsilon, minPts = MGmin)

  # Create a data frame with haplotype information and filter out the noise haplotype (0)
  preMGfile <- tibble::tibble(CHROM = vcf$"#CHROM", POS = vcf$POS, ID = rownames(LD), cluster = db_out$cluster) %>%
    dplyr::filter(cluster != 0)  # Remove rows where the haplotype is 0 (noise)

  # Check if there are multiple clusters
  if (base::length(base::unique(preMGfile$cluster)) > 1) {
    # Initialize a list to store the haplotype-specific tables
    cluster_tables <- base::list()

    # Loop through each haplotype and create the SNP table for each haplotype
    for (cluster_id in base::unique(preMGfile$cluster)) {
      # Get the SNPs that belong to the current cluster
      dbscan_cvel <- preMGfile %>% dplyr::filter(cluster == cluster_id) %>% tibble::as_tibble()

      # Filter the VCF data to include only SNPs for this haplotype
      cvel_vcf <- bin_vcf %>%
        tibble::rownames_to_column() %>%
        dplyr::filter(rowname %in% dbscan_cvel$ID) %>%
        tibble::column_to_rownames("rowname")

      # Calculate the start and end positions
      start_pos <- base::min(dbscan_cvel$POS)
      end_pos <- base::max(dbscan_cvel$POS)

      # Create the haplotype name
      hap_name <- base::paste0("hap_", start_pos, "_", end_pos)

      # Store the SNP table for this haplotype in the list
      cluster_tables[[hap_name]] <- cvel_vcf
    }
  } else {
    base::message(base::paste0("No haplotypes identified for epsilon = ", epsilon, ". Please try another value of epsilon or MGmin"))
  }
  ##Smooth the snps
  smooth_cluster_tables <- function(cluster_tables, LD, threshold = 2) {
    smoothed_tables <- base::list()  # Initialize a list to store the smoothed tables

    # Loop through each haplotype table
    for (cluster_id in base::names(cluster_tables)) {
      cluster_table <- cluster_tables[[cluster_id]]

      # Extract LD values for SNPs in this haplotype
      ld_values <- LD[base::rownames(LD) %in% base::rownames(cluster_table), base::colnames(LD) %in% base::rownames(cluster_table)]

      # Calculate the mean r² for each SNP in the haplotype
      mean_r2 <- base::apply(ld_values, 1, mean, na.rm = TRUE)

      # Compute the median and standard deviation of r² values
      median_r2 <- median(mean_r2)
      sd_r2 <- sd(mean_r2)

      # Identify SNPs that are outliers (more than 2 SDs away from the median)
      outliers <- base::which(base::abs(mean_r2 - median_r2) > threshold * sd_r2)

      # Create a smoothed table excluding outliers
      smoothed_table <- cluster_table[-outliers, , drop = FALSE]

      # Store the smoothed SNP table for this haplotype in the list
      smoothed_tables[[cluster_id]] <- smoothed_table
    }
    return(smoothed_tables)  # Return the list of smoothed tables
  }

  # Call the smoothing function
  smoothed_cluster_tables <- smooth_cluster_tables(cluster_tables, LD)

  #return either haplotypes with outliers removed or without, depending on the value of keep_outliers
  if (keep_outliers == FALSE) {
    smoothed_cluster_tables <- smoothed_cluster_tables[base::sapply(smoothed_cluster_tables, nrow) > 0] #remove haplotypes with no variants
    output <- smoothed_cluster_tables
    return(output)
  } else {
    cluster_tables <- cluster_tables[base::sapply(cluster_tables, nrow) > 0]
    output <- cluster_tables
    return(output)
  }
  return(output)
}
