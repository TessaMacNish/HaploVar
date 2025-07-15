utils::globalVariables(c("cluster", "rowname", "hap", "row_index", "pair"))

#' Identify Haplotype Variants
#'
#' This function requires a VCF and an LD matrix.
#' It will then define local haplotypes and identify the variants for each haplotype.
#' The output can be formatted in six ways, to be compatible with a wide range of GWAS and genomic selection tools.
#' @param vcf A VCF file.
#' @param LD A LD matrix file.
#' @param epsilon Affects haplotype size. It is a parameter of the DBSCAN clustering tool. The default is 0.6.
#' @param MGmin The minimum number of SNPs within a cluster for it to be defined as a haplotype.The default is 30.
#' @param minFreq The minimum number of individuals a haplotype variant must be present in to be considered a valid haplotype variant. The default is 2.
#' @param hetmiss_as Affects how missing data is handled for all instances where one allele in a genotype is missing.If hetmiss_as = "allele" the genotype is assumed to be heterozygous. If hetmiss_as = "miss" the genotype is treated as NA.
#' @param keep_outliers If FALSE removes SNPs, that are determined to be outliers.
#' @param format The output format. There are six different output formats (1,2,3,4,5,6).
#' @return A table of haplotype genotypes in your chosen format.
#' @export
#' @importFrom magrittr %>%
#' @importFrom stats median sd setNames
#' @importFrom dplyr across everything
haplotype_variants <- function(vcf, LD, epsilon = 0.6, MGmin = 30, minFreq = 2, hetmiss_as = "allele", keep_outliers = FALSE, format = 1) {
  ##format the LD file and vcf into usable formats
  bin_vcf <- dplyr::select(vcf, -c(1, 2, 4:9)) %>%
    tibble::column_to_rownames("ID") %>%
    dplyr::mutate_all(function(x) {
      switch(hetmiss_as,
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
  preMGfile <- tibble::tibble(CHROM = vcf$"#CHROM",POS = vcf$POS, ID = rownames(LD), cluster = db_out$cluster)
  # Check if there are multiple haplotypes
  if (base::length(base::unique(preMGfile$cluster)) > 1) {
    # Initialize a list to store the haplotype-specific tables
    cluster_tables <- base::list()

    # Loop through each cluster and create the SNP table for each haplotype
    for (cluster_id in base::unique(preMGfile$cluster)) {
      # Get the SNPs that belong to the current cluster
      dbscan_cvel <- preMGfile %>% dplyr::filter(cluster == cluster_id) %>% tibble::as_tibble()

      # Filter the VCF data to include only SNPs for this haplotype
      cvel_vcf <- bin_vcf %>%
        tibble::rownames_to_column() %>%
        dplyr::filter(rowname %in% dbscan_cvel$ID) %>%
        tibble::column_to_rownames("rowname")

      # Store the SNP table for this haplotype in the list
      cluster_tables[[as.character(cluster_id)]] <- cvel_vcf
    }
  } else {
    base::message(paste0("No haplotypes identified for epsilon = ", epsilon, ". Please try another value of epsilon or MGmin"))
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

  ##For each individual, for each haplotype find the alleles (by finding unique rows of SNPs)
  find_alleles <- function(genotype_table) {
    columns_to_split <- base::names(genotype_table)
    sep_char <- ifelse(any(base::grepl("/", base::as.matrix(genotype_table))), "/", "|")

    # Loop through each column and separate
    for (col in columns_to_split) {
      genotype_table <- genotype_table %>%
        tidyr::separate(col,
                        into = c(paste(col, "_strand1", sep = ""), paste(col, "_strand2", sep = "")),
                        sep = "[|/]",  # Use "|" or "/" as the separator
                        remove = TRUE)
    }
    # Transpose the table to make individuals the rows
    transposed_table <- base::t(genotype_table)
    transposed_df <- base::as.data.frame(transposed_table)

    # Generate haplotype variant labels
    labels <- base::LETTERS
    extra_labels <- base::outer(base::LETTERS, base::LETTERS, FUN = "paste0")
    hap_labels <- c(labels, extra_labels)
    hap_labels <- hap_labels[hap_labels != "NA"] #To make sure that no variants are named NA
    # Use dplyr to group by unique genotype patterns
    labeled_haplotypes <- transposed_df %>%
      dplyr::group_by_all() %>%
      dplyr::mutate(hap = hap_labels[dplyr::cur_group_id()]) %>%
      dplyr::ungroup() %>%
      dplyr::select(hap)  # Keep only the haplotype labels

    # Calculate haplotype frequencies
    hap_freq <- base::table(labeled_haplotypes$hap)

    # Identify valid haplotypes with frequency >= minFreq
    valid_haps <- base::names(hap_freq[hap_freq >= minFreq])

    # Create a named vector for valid haplotypes to assign consecutive labels
    hap_name_mapping <- setNames(hap_labels[1:base::length(valid_haps)], valid_haps)

    # Rename haplotypes so the variants are consecutive
    labeled_haplotypes$hap <- ifelse(labeled_haplotypes$hap %in% valid_haps,
                                     hap_name_mapping[labeled_haplotypes$hap],
                                     0)

    # Create a row index for grouping every two rows
    labeled_haplotypes <- labeled_haplotypes %>%
      dplyr::mutate(row_index = dplyr::row_number())

    # Create pairs of haplotypes
    paired_haplotypes <- labeled_haplotypes %>%
      dplyr::group_by(pair = (row_index - 1) %/% 2) %>%  # Group every two rows
      dplyr::summarize(hap = paste(hap, collapse = sep_char), .groups = 'drop') %>% # Combine with /
      dplyr::select(-pair)  # Drop the grouping column if not needed
    return(paired_haplotypes)
  }
  ## Create final haplotype table
  create_haplotype_table <- function(cluster_tables, preMGfile) {
    haplotype_table <- base::list()

    for (cluster_id in names(cluster_tables)) {
      cluster_vcf <- cluster_tables[[cluster_id]]
      cluster_snp_info <- preMGfile %>% dplyr::filter(cluster == as.numeric(cluster_id))

      if (nrow(cluster_snp_info) > 1) {
        hap_name <- paste0("hap_", base::min(cluster_snp_info$POS), "_", base::max(cluster_snp_info$POS))
        chromosome <- base::unique(cluster_snp_info$CHROM)[1]
        start_pos <- base::min(cluster_snp_info$POS)

        haplotype_alleles <- find_alleles(cluster_vcf)

        individual_alleles <- tibble::tibble(IDIV = base::colnames(cluster_vcf), hap_alleles = haplotype_alleles)
        hap_row <- tibble::tibble(MARKER = hap_name, CHROM = chromosome, POS = start_pos)
        hap_row <- base::cbind(hap_row, base::as.data.frame(base::t(individual_alleles$hap_alleles), stringsAsFactors = FALSE))
        base::colnames(hap_row)[4:base::ncol(hap_row)] <- individual_alleles$IDIV

        haplotype_table[[cluster_id]] <- hap_row
      }
    }


    return(haplotype_table)
  }

  ##Apply it to the cluster table or smoothed cluster table depending on the value of keep_outliers
  if (keep_outliers == FALSE) {
    smoothed_cluster_tables <- smoothed_cluster_tables[base::names(smoothed_cluster_tables) != "0"] #remove noise
    smoothed_cluster_tables <- smoothed_cluster_tables[base::sapply(smoothed_cluster_tables, nrow) > 0] #remove haplotypes with no variants
    haplotype_alleles <- base::lapply(smoothed_cluster_tables, find_alleles)
    haplotype_table_list <- create_haplotype_table(smoothed_cluster_tables, preMGfile)
  } else {
    cluster_tables <- cluster_tables[base::names(cluster_tables) != "0"]
    cluster_tables <- cluster_tables[base::sapply(cluster_tables, nrow) > 0]
    haplotype_alleles <- base::lapply(cluster_tables, find_alleles)
    haplotype_table_list <- create_haplotype_table(cluster_tables, preMGfile)
  }
  ##make a table with the individual, and the haplotype and variants
  final_haplotype_table <- base::data.frame()
  for (haplotype_table_id in base::names(haplotype_table_list)) {
    haplotype_table_entry <- haplotype_table_list[[haplotype_table_id]]

    # Append the new entry to the final table
    final_haplotype_table <- base::rbind(final_haplotype_table, haplotype_table_entry)
  }
  # Remove rows where the IDIV column starts with "IDIV"
  final_haplotype_table <- final_haplotype_table[!base::grepl("^IDIV", base::rownames(final_haplotype_table)), ]
  base::rownames(final_haplotype_table) <- NULL
  final_haplotype_table[base::is.na(final_haplotype_table)] <- 0

  ##Define function to find the dosage of each haplotype variant
  # Function to process the table
  dosage_haplotype_variants <- function(haplotype_table) {
    # Initialize an empty data frame for the result
    dosage <- base::data.frame(MARKER = character(),
                               CHROM = integer(),
                               POS = integer(),
                               stringsAsFactors = FALSE)

    # Extract unique individual IDs
    individuals <- base::colnames(haplotype_table)[-(1:3)]

    # Gather all unique variants across the entire table, excluding "0"
    all_variants <- base::setdiff(base::unique(base::unlist(base::strsplit(base::unlist(haplotype_table[individuals]), "[|/]"))), "0")

    # Loop over each row
    for (i in 1:base::nrow(haplotype_table)) {
      # Get marker information
      marker <- haplotype_table$MARKER[i]
      chrom <- haplotype_table$CHROM[i]
      pos <- haplotype_table$POS[i]

      # Create an empty list to store variant counts per individual
      variant_counts <- base::list()

      # Loop over each individual
      for (indiv in individuals) {
        # Split the variant by | or /
        variants <- base::unlist(base::strsplit(haplotype_table[i, indiv], "[|/]"))
        variants <- variants[variants != "0"]

        # Count occurrences of each unique variant found in the table
        count_table <- base::table(base::factor(variants, levels = all_variants))

        # Store counts with individual ID suffix
        for (variant in base::names(count_table)) {
          variant_counts[[paste0(indiv, "_", variant)]] <- base::as.integer(count_table[variant])
        }
      }

      # Create rows for each detected variant
      for (variant in all_variants) {
        # Only include variants that are present in the row
        if (base::sum(base::unlist(variant_counts[base::grep(paste0("_", variant), base::names(variant_counts))])) > 0) {
          # Create a new row
          new_row <- base::data.frame(
            MARKER = paste(marker, variant, sep = "_"),
            CHROM = chrom,
            POS = pos,
            stringsAsFactors = FALSE
          )

          # Add counts for each individual for this variant
          for (indiv in individuals) {
            new_row[[indiv]] <- variant_counts[[paste0(indiv, "_", variant)]]
          }

          # Append the row to the result
          dosage <- base::rbind(dosage, new_row)
        }
      }
    }

    return(dosage)
  }

  # Function to process the table in strand-specific pseudo VCF format
  vcf_format <- function(final_haplotype_table, vcf) {
    # Initialize an empty data frame for the result
    pseudo_vcf <- base::data.frame(
      CHROM = character(),
      POS = integer(),
      ID = character(),
      REF = character(),
      ALT = character(),
      QUAL = character(),
      FILTER = character(),
      INFO = character(),
      FORMAT = character(),
      stringsAsFactors = FALSE
    )

    # Extract unique individual IDs
    individuals <- base::colnames(final_haplotype_table)[-(1:3)]

    # Gather all unique variants across the entire table
    all_variants <- base::setdiff(base::unique(base::unlist(base::strsplit(base::unlist(final_haplotype_table[individuals]), "[|/]"))), "0")

    # Loop over each row
    for (i in 1:base::nrow(final_haplotype_table)) {
      # Get marker information
      marker <- final_haplotype_table$MARKER[i]
      chrom <- final_haplotype_table$CHROM[i]
      pos <- final_haplotype_table$POS[i]

      # Get QUAL, FILTER, INFO, and FORMAT from the VCF using POS in the final_haplotype_table
      vcf_row <- vcf[vcf$POS == pos, ]
      if (base::nrow(vcf_row) > 0) {
        qual <- vcf_row$QUAL
        filter <- vcf_row$FILTER
        info <- vcf_row$INFO
        format <- vcf_row$FORMAT
      } else {
        # If no matching POS, assign NA or default values
        qual <- "."
        filter <- "."
        info <- "."
        format <- "."
      }

      # Loop over each variant in the row
      for (variant in all_variants) {
        # Create an empty list
        presence_absence <- base::list()
        variant_present <- FALSE  # Track if variant is present in any individual

        # Loop over each individual
        for (indiv in individuals) {
          # Get the genotype format (| or /) from the original data
          genotype <- final_haplotype_table[i, indiv]
          delimiter <- ifelse(base::grepl("\\|", genotype), "|", "/")

          # Split the genotype into its haplotypes based on | or /
          haplotypes <- base::strsplit(genotype, "[|/]")[[1]]
          strand_encoding <- c(0, 0)

          # Check each strand for the specific variant and set to 1 if the variant is present
          for (strand in seq_along(haplotypes)) {
            if (haplotypes[strand] == variant) {
              strand_encoding[strand] <- 1
              variant_present <- TRUE
            }
          }

          # Combine strand encodings into a single value separated by the original delimiter
          presence_absence[[indiv]] <- paste(strand_encoding, collapse = delimiter)
        }

        # Only add the row if the variant is present in at least one individual
        if (variant_present) {
          new_row <- base::data.frame(
            CHROM = chrom,
            POS = pos,
            ID = paste(marker, variant, sep = "_"),
            REF = ".",
            ALT = ".",
            QUAL = qual,
            FILTER = filter,
            INFO = info,
            FORMAT = format,
            stringsAsFactors = FALSE
          )

          # Add individual presence-absence data
          for (indiv in individuals) {
            new_row[[indiv]] <- presence_absence[[indiv]]
          }

          # Append the row to the result
          pseudo_vcf <- base::rbind(pseudo_vcf, new_row)
        }
      }
    }
    base::colnames(pseudo_vcf)[base::colnames(pseudo_vcf) == "CHROM"] <- "#CHROM"
    return(pseudo_vcf)
  }


  ##Format the output
  if (format == 2){

    # Apply the function
    formatted_haplotype_table <- dosage_haplotype_variants(final_haplotype_table)

  } else if (format == 3){
    dosage_haplotype_table <- dosage_haplotype_variants(final_haplotype_table)
    dosage_haplotype_table2 <- dplyr::select(dosage_haplotype_table, -c(1:3))
    dosage_haplotype_table3 <- base::t(dosage_haplotype_table2)
    formatted_haplotype_table <- as.data.frame(dosage_haplotype_table3)
    base::colnames(formatted_haplotype_table) <- dosage_haplotype_table$MARKER
  } else if (format == 4){
    # Apply the function
    formatted_haplotype_table <- dosage_haplotype_variants(final_haplotype_table)
    # Map values: 0 -> -1, 1 -> 0, 2 -> 1
    formatted_haplotype_table <- formatted_haplotype_table %>%
      dplyr::mutate(across(-c(1:3), ~ dplyr::case_when(
        . == 0 ~ -1,
        . == 1 ~ 0,
        . == 2 ~ 1,
        TRUE ~ .
      )))
  } else if (format == 5){
    dosage_haplotype_table <- dosage_haplotype_variants(final_haplotype_table)
    dosage_haplotype_table2 <- dplyr::select(dosage_haplotype_table, -c(1:3))
    dosage_haplotype_table3 <- base::t(dosage_haplotype_table2)
    formatted_haplotype_table <- as.data.frame(dosage_haplotype_table3)
    colnames(formatted_haplotype_table) <- dosage_haplotype_table$MARKER
    # Map values: 0 -> -1, 1 -> 0, 2 -> 1
    formatted_haplotype_table <- formatted_haplotype_table %>%
      dplyr::mutate(across( everything(), ~ dplyr::case_when(
        . == 0 ~ -1,
        . == 1 ~ 0,
        . == 2 ~ 1,
        TRUE ~ .
      )))
  } else if (format == 6){
    formatted_haplotype_table <- vcf_format(final_haplotype_table, vcf)
  }
  else {
    formatted_haplotype_table <- final_haplotype_table
  }

  ##print results
  return(formatted_haplotype_table)

}
