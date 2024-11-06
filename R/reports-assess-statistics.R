
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
