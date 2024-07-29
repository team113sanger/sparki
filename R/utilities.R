merge_sample <- function(std_subset, mpa_subset) {

    std_subset[, COLNAME_STD_HIERARCHY] <- sapply(seq_len(nrow(std_subset)), function(x) {

        # Get rank and taxon at a given line.
        rank_std <- std_subset[, COLNAME_STD_RANK][x]
        taxon_std <- std_subset[, COLNAME_STD_TAXON][x]

        # For subranks (i.e. ranks followed by a number), return NA,
        # as subranks are not present in MPA-style reports.
        if (!(rank_std %in% mpa_subset[, COLNAME_MPA_RANK])) return(NA)

        taxcol_mpa <- names(get_association(rank_std))

        pos <- which(
            (mpa_subset[, COLNAME_MPA_RANK] == rank_std) & 
            (mpa_subset[, taxcol_mpa] == taxon_std)
        )

        return(mpa_subset[, COLNAME_MPA_TAXON][pos])
    })

    return(std_subset)
}

# Look at last rank and taxon name to compare
mergeReports <- function(std_report, mpa_report, verbose) {

    if (verbose) {
        pb <- txtProgressBar(min = 0, max = length(unique(std_report[, COLNAME_STD_SAMPLE])), style = 3)
        i <- 1
    }

    updated_std_report <- data.frame(matrix(nrow = 0, ncol = (ncol(std_report) + 1)))

    for (sample in unique(std_report[, COLNAME_STD_SAMPLE])) {

        if (verbose) {
            setTxtProgressBar(pb, i)
            i <- i + 1
        }

        std_subset <- std_report[std_report[, COLNAME_STD_SAMPLE] == sample, ]
        mpa_subset <- mpa_report[mpa_report[, COLNAME_MPA_SAMPLE] == sample, ]

        updated_std_subset <- merge_sample(std_subset, mpa_subset)

        updated_std_report <- rbind(updated_std_report, updated_std_subset)
    }
    if (verbose) cat("\n")

    return(updated_std_report)
}

#########################################################
## UTILITY FUNCTIONS FOR CREATING AUXILIARY DATAFRAMES ##
#######################################################################################################

get_nDomainReads <- function(report) {

    # Identify if report follows a standard or MPA format; there are some slight 
    # differences in how domains are written depending on the format!
    if (is_mpa(report)) {

        # In an MPA-style report, the domain will have a prefix ("d__"). Additionally,
        # since the MPA format shows the hierarchy of each taxon, we also need to select
        # specifically the lines that have only the domain-level information.
        all_domains <- "d__Eukaryota$|d__Viruses$|d__Archaea$|d__Bacteria$"
        subset_domains <- "d__Viruses$|d__Archaea$|d__Bacteria$"
        colname_sample <- COLNAME_MPA_SAMPLE
        
    } else {

        # In a standard report, each taxon is shown independently from their hierarchy,
        # therefore it is easier to select the lines that contain domain-level information.
        all_domains <- "Eukaryota|Viruses|Archaea|Bacteria"
        subset_domains <- "Viruses|Archaea|Bacteria"
        colname_sample <- COLNAME_STD_SAMPLE
    }

    n_domainReads <- data.frame(matrix(nrow = 0, ncol = 4))

    # Iterate over samples in report...
    for (sample in unique(report[, colname_sample])) {

        # Subset report to keep only a given sample.
        report_sample <- report[report[, colname_sample] == sample,]

        # Here we create a named vector to establish an association between the domain vectors
        # we created above and what is their scope (all domains or only non-eukaryote domains).
        scopes <- c(all_domains, subset_domains)
        names(scopes) <- c("all", "non-eukaryote")

        # Sum all sample-level reads classified under all domains ("all") 
        # or under non-eukaryote domains ("non-eukaryote").
        for (scope in scopes) {

            n_domainReads <- rbind(
                n_domainReads, 
                c(sample, names(scopes)[scopes == scope], sum_domainReads(report = report_sample, domains = scope))
            )
        }
    }

    colnames(n_domainReads) <- c(
        COLNAME_DOMAIN_READS_SAMPLE, 
        COLNAME_DOMAIN_READS_SCOPE, 
        COLNAME_DOMAIN_READS_N_READS
    )
    
    n_domainReads[, COLNAME_DOMAIN_READS_N_READS] <- as.numeric(
        n_domainReads[, COLNAME_DOMAIN_READS_N_READS]
    )

    return(n_domainReads)
}

get_nTaxa <- function(report, ranks) {

    if (!("domain" %in% colnames(report))) {
        stop(paste0(
            "The column 'domain' was not found in the report provided. If your report has MPA format, ",
            "please make sure you run addConciseTaxon(report, ranks = c('D')) before running get_nTaxa(). ",
            "Alternatively, if your report has standard format, please make sure you run addDomain() before ",
            "running get_nTaxa()."
        ))
    }

    if (is_mpa(report)) {
        colname_sample <- COLNAME_MPA_SAMPLE
        colname_taxon <- COLNAME_MPA_TAXON
        colname_rank <- COLNAME_MPA_RANK
        colname_domain <- "domain"

    } else {
        colname_sample <- COLNAME_STD_SAMPLE
        colname_taxon <- COLNAME_STD_TAXON
        colname_rank <- COLNAME_STD_RANK
        colname_domain <- COLNAME_STD_DOMAIN
    }

    n_taxa_in_samples <- data.frame(matrix(nrow = 0, ncol = 4))

    # Iterate over samples...
    for (sample in unique(report[, colname_sample])) {

        # Subset report to keep only a given sample.
        report_sample <- report[report[, colname_sample] == sample,]

        # Iterate over domains...
        for (domain in unique(report[, colname_domain])) {

            # Iterate over ranks...
            for (rank in ranks) {

                n_taxa <- length(unique(
                    report_sample[, colname_taxon][report[, colname_rank] == rank & report[, colname_domain] == domain]
                ))
                n_taxa_in_samples <- rbind(n_taxa_in_samples, c(sample, domain, rank, n_taxa))

            }
        }
    }

    colnames(n_taxa_in_samples) <- c("sample", "domain", "rank", "value")

    return(n_taxa_in_samples)
}


getClassificationSummary <- function(report) {

    if (is_mpa(report)) {
        stop(paste0(
            "This function is not applicable to MPA-style reports. The purpose of this function ",
            "is to assess the number of classified and unclassified reads for each sample, but the ",
            "MPA-style reports do not contain information on unclassified reads."
        ))
    } else {
        colname_n_frag_clade <- COLNAME_STD_N_FRAG_CLADE
        colname_taxon <- COLNAME_STD_TAXON
        colname_sample <- COLNAME_STD_SAMPLE
    }

    class_unclass_df <- data.frame(matrix(nrow = 0, ncol = 3))

    classifications <- list(
        "Unclassified" = "unclassified",
        "Classified" = "root"
    )

    for (sample in unique(report[, colname_sample])) {
        for (classification in names(classifications)) {

            value <- as.numeric(
                report[, colname_n_frag_clade][report[, colname_taxon] == classifications[[classification]] & report[, colname_sample] == sample]
            )

            class_unclass_df <- rbind(
                class_unclass_df, 
                c(sample, classification, value)
            )
        }
    }

    colnames(class_unclass_df) <- c(
        COLNAME_CLASSIF_SUMMARY_SAMPLE, 
        COLNAME_CLASSIF_SUMMARY_READ_TYPE, 
        COLNAME_CLASSIF_SUMMARY_N_READS
    )

    class_unclass_df[, COLNAME_CLASSIF_SUMMARY_N_READS] <- as.numeric(class_unclass_df[, COLNAME_CLASSIF_SUMMARY_N_READS])

    return(class_unclass_df)
}


getProportion <- function(report, taxon) {

    if (is_mpa(report)) {
        stop(paste0(
            "This function is not applicable to MPA-style reports. The purpose of this function ",
            "is to assess the proportion of classified reads relative to the total number of reads ",
            "assessed, but the MPA-style reports do not contain information on unclassified reads."
        ))
    } else {
        colname_n_frag_clade <- COLNAME_STD_N_FRAG_CLADE
        colname_taxon <- COLNAME_STD_TAXON
        colname_sample <- COLNAME_STD_SAMPLE
    }

    class_unclass_df <- getClassificationSummary(report)

    for (sample in unique(class_unclass_df[, COLNAME_CLASSIF_SUMMARY_SAMPLE])) {
        subset <- class_unclass_df[, class_unclass_df[, COLNAME_CLASSIF_SUMMARY_SAMPLE] == sample, ]
        
        classified_reads <- subset[, COLNAME_CLASSIF_SUMMARY_N_READS][subset[, COLNAME_CLASSIF_SUMMARY_READ_TYPE] == "Classified"]
        unclassified_reads <- subset[, COLNAME_CLASSIF_SUMMARY_N_READS][subset[, COLNAME_CLASSIF_SUMMARY_READ_TYPE] == "Unclassified"]

        total_reads <- classified_reads + unclassified_reads

    }



    return(class_unclass_df)
}


#############################################################################
## UTILITY FUNCTIONS FOR ASSESSING THE STATISTICAL SIGNIFICANCE OF RESULTS ##
#######################################################################################################

assess_ratioMinimisers <- function(report) {

    if (is_mpa(report)) {
        stop(paste0("This function does not support MPA-style reports. Please provide a standard report."))
    } 
    
    report[, COLNAME_STD_RATIO_TAXON] <- (report[, COLNAME_STD_UNIQ_MINIMISERS] / report[, COLNAME_STD_DB_MINIMISERS_TAXON])
    report[, COLNAME_STD_RATIO_CLADE] <- (report[, COLNAME_STD_UNIQ_MINIMISERS] / report[, COLNAME_STD_DB_MINIMISERS_CLADE])

    return(report)
}


#' ASSESS THE STATISTICAL SIGNIFICANCE OF RESULTS
#' 
#' @param report
#' @param ref_db 
#' @return An updated version of the input dataframe, with new columns containing statistical significance results.
assess_statSig <- function(report, ref_db, verbose = TRUE) {

    if (is_mpa(report)) stop(paste0("This function does not support MPA-style reports."))

    # Calculate p-values.
    report[, COLNAME_STD_PVALUE] <- calculate_p_value(
        report[[COLNAME_STD_UNIQ_MINIMISERS]], 
        report[[COLNAME_STD_DB_MINIMISERS_CLADE]], 
        report[[COLNAME_STD_SAMPLE_SIZE]]
    )

    # Add number of taxa identified in a given rank for a given sample.
    report <- add_nTaxaInRank(report)

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

    return(report)
}
