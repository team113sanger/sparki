
get_name_at_rank <- function(report, species, rank) {

    # Get the species position (row number) in the report (first appearance only). 
    species_position <- which(report[, COLNAME_STD_TAXON] == species)[1]

    # Find the corresponding rank position (row number) in the report (first appearance only).
    i <- 1
    while (i > 0) {

        # If the desired rank has been found, break the loop.
        if (report[, COLNAME_STD_RANK, drop = TRUE][(species_position - i)] == rank) {
            rank_position <- (species_position - i)
            break
        } else {
            i <- i + 1
        }
    }

    subset_report <- report |> dplyr::slice(rank_position)
    
    return(subset_report[[COLNAME_STD_TAXON]])
}

subsetReports <- function(report, species, verbose) {

    if (!(species %in% report[[COLNAME_STD_TAXON]])) {
        stop(paste0(
            "The organism ", species, " is not present ",
            "in the dataframe. Please review your input!"
        ))
    }

    # Determine the genus and family of an user-defined species.
    genus <- get_name_at_rank(report, species, "G")
    family <- get_name_at_rank(report, species, "F")

    # Create vector with species, genus and family to be filtered out.
    taxa_to_remove <- c(species, genus, family)

    if (verbose == TRUE) {
        warning(paste0(
            "Filtering out ", species, ", ", genus, ", and ", family, "..."
        ))
    }

    report <- report |> 
        # Select rows corresponding to species, genus and family.
        dplyr::filter(!!as.name(COLNAME_STD_RANK) %in% c("F", "G", "S")) |>
        # Filter out the user-defined organism.
        dplyr::filter(!(!!as.name(COLNAME_STD_TAXON) %in% taxa_to_remove))
    
    return(report)
}

#' MERGE STANDARD AND MPA-STYLE REPORTS
#' 
#' This function takes a standard and an MPA-style dataframe, loaded with load_STDreports() and
#' load_MPAreports() respectively, and merges the latter into the former.
#' 
#' @param std_reports A dataframe (tibble) with standard reports, loaded with load_STDreports().
#' @param mpa_reports A dataframe (tibble) with MPA-style reports, loaded with load_MPAreports().
#' @return An updated version of the dataframe (tibble) with standard reports, now containing the
#' information from the dataframe (tibble) with MPA-style reports.
#' @export
#'
mergeReports <- function(std_reports, mpa_reports, samples_to_remove) {

    mpa_reports <- mpa_reports |> 
        # Rename column for consistency with the standard report.
        dplyr::rename(!!COLNAME_STD_TAXON := !!COLNAME_MPA_TAXON_LEAF) |>
        # Drop columns that are already present in the standard report.
        dplyr::select(!c(!!COLNAME_MPA_N_FRAG_CLADE, !!COLNAME_MPA_RANK)) 
        
    # The primary keys that will be used for joining are "taxon leaf"/"name" and "sample_id"
    merged_reports <- dplyr::left_join( 
        std_reports, mpa_reports,
        by = dplyr::join_by(!!COLNAME_STD_TAXON, !!COLNAME_STD_SAMPLE)
    )

    if (!missing(samples_to_remove)) {
        merged_reports <- merged_reports |>
            dplyr::filter(!(!!as.name(COLNAME_STD_SAMPLE) %in% samples_to_remove))
    } else {
        warning("No samples were filtered out from the merged reports table.")
    }

    return(merged_reports)
}

assessMinimiserRatio <- function(report) {

    if (is_mpa(report)) {
        stop(paste0("This function does not support MPA-style reports. Please provide a standard report."))
    } 
    
    report[, COLNAME_STD_RATIO_TAXON] <- (report[, COLNAME_STD_UNIQ_MINIMISERS] / report[, COLNAME_STD_DB_MINIMISERS_TAXON])
    report[, COLNAME_STD_RATIO_CLADE] <- (report[, COLNAME_STD_UNIQ_MINIMISERS] / report[, COLNAME_STD_DB_MINIMISERS_CLADE])

    return(report)
}


###############################################
## HELPER FUNCTIONS FOR STATISTICAL ANALYSES ##
#######################################################################################################

calculate_p_value <- function(sample_n_uniq_minimisers_taxon, db_n_minimisers_taxon, sample_size, ref_db) {

    # Get proportion of clade-level minimisers of a given taxon in the reference database (DB).
    # This is the same as the probability of getting this taxon from the database.
    p_clade_in_db <- (db_n_minimisers_taxon / sum(ref_db[, COLNAME_REF_DB_MINIMISERS_CLADE]))
    p_clade_in_db <- as.numeric(p_clade_in_db)

    # Mean is the total number of reads analysed from the sample (sample size) times the probability of success which 
    # is equal to the proportion of clade-level minimisers of the taxon out of the total available in the reference DB.
    mean <- (sample_size * p_clade_in_db)
    mean <- as.numeric(mean)

    # Standard deviation = sqrt(n*P*(1-p))
    sdev <- sqrt(sample_size * p_clade_in_db * (1 - p_clade_in_db)) 
    sdev <- as.numeric(sdev)

    # Calculate p-values.
    pval <- pnorm(
        q = as.numeric(sample_n_uniq_minimisers_taxon), 
        mean = mean, 
        sd = sdev, 
        lower.tail = FALSE
    )

    pval <- as.numeric(pval)
            
    return(pval)
}

#' ASSESS THE STATISTICAL SIGNIFICANCE OF RESULTS
#' 
#' @param report
#' @param ref_db 
#' @return An updated version of the input dataframe, with new columns containing statistical significance results.
assessStatistics <- function(report, ref_db, verbose = TRUE) {

    if (is_mpa(report)) stop(paste0("This function does not support MPA-style reports."))

    if (verbose) cat("Calculating p-values...\n")

    # Calculate p-values.
    report[, COLNAME_STD_PVALUE] <- calculate_p_value(
        report[[COLNAME_STD_UNIQ_MINIMISERS]], 
        report[[COLNAME_STD_DB_MINIMISERS_CLADE]], 
        report[[COLNAME_STD_SAMPLE_SIZE]],
        ref_db
    )

    # Add number of taxa identified in a given rank for a given sample.
    report <- add_nTaxaInRank(report)

    if (verbose) cat("Adjusting p-values...\n")

    # Perform p-value correction.
    report[, COLNAME_STD_PADJ] <- sapply(seq_len(nrow(report)), function(x) {

        p.adjust(
            report[[COLNAME_STD_PVALUE]][x], 
            method = "BH",
            n = report[["n_taxa_in_rank"]][x]
        )
    })

    # Add column stating whether results are significant or not based on the adjusted p-value.
    report[, COLNAME_STD_SIGNIF] <- ifelse(
        report[, COLNAME_STD_PADJ] <= 0.05,
        "Significant",
        "Non-significant"
    )

    if (verbose) cat("Successfully completed.\n")

    return(report)
}



#####################################################
## UTILITY FUNCTIONS FOR ADDING COLUMNS TO REPORTS ##
#######################################################################################################

addRank <- function(taxon) {

    # Define rank prefixes and corresponding letters.
    rank_prefixes <- c("d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__")
    rank_letters <- c("D", "K", "P", "C", "O", "F", "G", "S")

    # Create a single regex pattern.
    pattern <- paste0("(", paste(rank_prefixes, collapse = "|"), ")")

    # Extract all rank prefixes from each taxon string.
    all_matches <- stringr::str_extract_all(taxon, pattern)

    # Find the last (most specific) rank for each taxon.
    last_ranks <- sapply(all_matches, function(x) tail(x, 1))

    # Map rank prefixes to letters.
    rank <- rank_letters[match(last_ranks, rank_prefixes)]

    return(rank)
}

#' CALCULATE SAMPLE SIZES
#' 
#' @param report A report.
#' @return An updated version of the input dataframe, with a new column containing sample sizes.
#' @export
#' 
addSampleSize <- function(report) {

    # If report is in MPA format...
    if (is_mpa(report)) {
        stop(paste0(
            "This function is not applicable to MPA-style reports. The purpose of this function ",
            "is to assess the total number of reads (classified + unclassified) assessed in each ",
            "sample (i.e., sample size), but the MPA-style reports do not contain information on ",
            "unclassified reads. If you would like to have the total number of reads in your MPA-style ",
            "report, please do:\n\nstd_reports <- add_nReads(std_reports)\n", "mpa_reports <- ",
            "transfer_nReads(mpa_reports, std_reports)\n\n"
        ))
    }

    subset <- report |> 
        # Select columns of interest.
        dplyr::select(!!COLNAME_STD_SAMPLE, !!COLNAME_STD_TAXON, !!COLNAME_STD_N_FRAG_CLADE) |>
        # Select rows of interest.
        dplyr::filter(!!as.name(COLNAME_STD_TAXON) %in% c("unclassified", "root")) |>
        # Reformat dataframe into a more convenient shape.
        tidyr::pivot_wider(
            names_from = !!COLNAME_STD_TAXON, 
            values_from = !!COLNAME_STD_N_FRAG_CLADE
        ) |>
        # Rename columns.
        dplyr::rename(
            Unclassified = unclassified,
            Classified = root
        ) |>
        # Create column with sample sizes (i.e. sum of unclasified and classified reads).
        dplyr::mutate(!!COLNAME_STD_SAMPLE_SIZE := rowSums(dplyr::across((where(is.numeric))))) |>
        # Select final columns.
        dplyr::select(!!COLNAME_STD_SAMPLE, !!COLNAME_STD_SAMPLE_SIZE)

    # Add sample sizes to report.
    report <- dplyr::left_join(report, subset, by = dplyr::join_by(!!COLNAME_STD_SAMPLE)) |>
        # Reorder columns.
        dplyr::relocate(!!COLNAME_STD_SAMPLE_SIZE, .after = !!COLNAME_STD_SAMPLE)

    return(report)
}

addMinimiserData <- function(report, ref_db) {

    if (is_mpa(report)) stop(paste0("This function does not support MPA-style reports."))

    # Select relevant columns.
    ref_db <- ref_db |> 
        dplyr::select(
            !!COLNAME_REF_DB_NCBI_ID,
            !!COLNAME_REF_DB_MINIMISERS_TAXON,
            !!COLNAME_REF_DB_MINIMISERS_CLADE
        )

    report <- report |>
        # Add DB info to the report.
        dplyr::left_join(ref_db, by = dplyr::join_by(!!COLNAME_STD_NCBI_ID)) |>
        # Rename columns.
        dplyr::rename(
            !!COLNAME_STD_DB_MINIMISERS_CLADE := !!COLNAME_REF_DB_MINIMISERS_CLADE,
            !!COLNAME_STD_DB_MINIMISERS_TAXON := !!COLNAME_REF_DB_MINIMISERS_TAXON
        )

    return(report)
}

add_nTaxaInRank <- function(report) {

    n_taxa_in_rank <- report |>
        # Create column with sample+rank combinations.
        dplyr::mutate(
            temp = paste(!!as.name(COLNAME_STD_SAMPLE), !!as.name(COLNAME_STD_RANK), sep = "_")
        ) |>
        # Count sample+rank occurrences (this will tell how many taxa are present in each rank
        # for each sample).
        dplyr::count(temp)

    report <- report |> 
        # Create temporary column with sample+rank combinations.
        dplyr::mutate(
            temp = paste(!!as.name(COLNAME_STD_SAMPLE), !!as.name(COLNAME_STD_RANK), sep = "_")
        ) |>
        # Merge dataframes.
        dplyr::left_join(n_taxa_in_rank, by = dplyr::join_by(temp)) |>
        # Remove the temporary column.
        dplyr::select(!temp) |>
        # Rename column.
        dplyr::rename(!!COLNAME_STD_N_TAXA_RANK := n)

    report[[COLNAME_STD_N_TAXA_RANK]] <- as.numeric(report[[COLNAME_STD_N_TAXA_RANK]])

    return(report)
}
