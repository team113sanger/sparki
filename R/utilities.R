merge_sample <- function(std_subset, mpa_subset) {

    std_subset[, COLNAME_STD_HIERARCHY] <- sapply(seq_len(nrow(std_subset)), function(x) {

        rank_std <- std_subset[, COLNAME_STD_RANK][x]

        if (!(rank_std %in% mpa_subset[, COLNAME_MPA_RANK])) return(NA)

        taxon_std <- std_subset[, COLNAME_STD_TAXON][x]
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
mergeReports <- function(std_report, mpa_report) {

    pb <- txtProgressBar(min = 0, max = length(unique(std_report[, COLNAME_STD_SAMPLE])), style = 3)

    updated_std_report <- data.frame(matrix(nrow = 0, ncol = (ncol(std_report) + 1)))

    i <- 1
    for (sample in unique(std_report[, COLNAME_STD_SAMPLE])) {

        setTxtProgressBar(pb, i)

        std_subset <- std_report[std_report[, COLNAME_STD_SAMPLE] == sample, ]
        mpa_subset <- mpa_report[mpa_report[, COLNAME_MPA_SAMPLE] == sample, ]

        updated_std_subset <- merge_sample(std_subset, mpa_subset)

        updated_std_report <- rbind(updated_std_report, updated_std_subset)

        i <- i + 1
    }
    cat("\n")

    return(updated_std_report)
}

#####################################################
## UTILITY FUNCTIONS FOR ADDING COLUMNS TO REPORTS ##
#######################################################################################################

addRank <- function(report, verbose = TRUE) {

    if (!(is_mpa(report))) {
        stop(paste0(
            "This function is not applicable to a standard report (i.e. not MPA-style). ",
            "The standard report format already has the column ", COLNAME_STD_RANK, ", which ", 
            "indicates the taxonomic rank of each line."
        ))   
    }

    if (verbose) cat("Adding", COLNAME_MPA_RANK, "column to the MPA-style report...\n")

    # Define rank prefixes and corresponding letters.
    rank_prefixes <- c("d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__")
    rank_letters <- c("D", "K", "P", "C", "O", "F", "G", "S")

    # Create a single regex pattern.
    pattern <- paste0("(", paste(rank_prefixes, collapse = "|"), ")")

    # Extract all rank prefixes from each taxon string.
    all_matches <- stringr::str_extract_all(report[, COLNAME_MPA_TAXON], pattern)

    # Find the last (most specific) rank for each taxon.
    last_ranks <- sapply(all_matches, function(x) tail(x, 1))

    # Map rank prefixes to letters.
    report[, COLNAME_MPA_RANK] <- rank_letters[match(last_ranks, rank_prefixes)]

    if (verbose) cat("\tRank column added successfully.\n")

    return(report)
}

addConciseTaxon <- function(report, verbose = TRUE) {

    if(!(is_mpa(report))) {
        stop(paste0(
            "This function is not applicable to a standard report (i.e. not MPA-style). The standard report ",
            "format already has the column ", COLNAME_STD_TAXON, ", which has concise taxon names instead of ",
            "taxon hierarchies."
        )) 
    }

    if (verbose) cat("Adding columns with taxon names for all ranks...\n")

    ranks <- c("D", "K", "P", "C", "O", "F", "G", "S")
    ranks <- get_association(ranks)
    
    # Iterate over ranks... 
    for (rank in ranks) {

        # Get name for column that will be added to the MPA-style report.
        colname <- names(ranks)[ranks == rank]

        if (verbose) cat("\tProcessing rank:", colname, "...\n")

        # Create column with concise taxon name at a given rank.
        # Examples:
        #
        # "d__Eukaryota" -> "Eukaryota" (under "domain")
        #
        # "(...)f__Papillomaviridae|g__Alphapapillomavirus" -> "Papillomaviridae" (under "family") and 
        # "Alphapapillomavirus" (under "genus")
        #
        # "(...)g__Homo|s__Homo sapiens" -> "Homo" (under "genus") and "Homo sapiens" (under "species")
        report[, colname] <- sapply(seq_len(nrow(report)), function(x) {

            ifelse(
                report[, COLNAME_MPA_RANK][x] == rank,
                return(extract_taxon(report[, COLNAME_MPA_TAXON][x], rank = rank, last_in_hierarchy = TRUE)),
                return(extract_taxon(report[, COLNAME_MPA_TAXON][x], rank = rank, last_in_hierarchy = FALSE))
            )
        })
    }

    if (verbose) cat("\tAll columns added successfully.\n")

    return(report)
}

#' ASSESS THE TOTAL NUMBER OF READS ANALYSED FOR EACH SAMPLE
#' 
#' @param merged_reports
#' @return An updated version of the input dataframe, with a new column containing sample sizes.
add_nReads <- function(report) {

    # If report is in MPA format...
    if (is_mpa(report)) {
        stop(paste0(
            "This function is not applicable to MPA-style reports. The purpose of this function ",
            "is to assess the total number of reads (classified + unclassified) assessed in each ",
            "sample, but the MPA-style reports do not contain information on unclassified reads. ",
            "If you would like to have the total number of reads in your MPA-style report, please do:\n\n",
            "std_reports <- add_nReads(std_reports)\n", "mpa_reports <- transfer_nReads(mpa_reports, ",
            "std_reports)\n\n"
        ))
    }

    report[, COLNAME_STD_TOTAL_READS] <- NA

    # Get numbers of classified/unclassified reads for each sample.
    subset <- report[report[, COLNAME_STD_TAXON] %in% c("unclassified", "root"), ]

    for (sample in unique(subset[, COLNAME_STD_SAMPLE])) {

        subset_sample <- subset[subset[, COLNAME_STD_SAMPLE] == sample, ]

        n_unclassified_reads <- as.numeric(
            subset_sample[, COLNAME_STD_N_FRAG_CLADE][subset_sample[, COLNAME_STD_TAXON] == "unclassified"]
        )
        n_classified_reads <- as.numeric(
            subset_sample[, COLNAME_STD_N_FRAG_CLADE][subset_sample[, COLNAME_STD_TAXON] == "root"]
        )

        report[, COLNAME_STD_TOTAL_READS][report[, COLNAME_STD_SAMPLE] == sample] <- (n_unclassified_reads + n_classified_reads)
    }

    return(report)
}

add_DBinfo <- function(report, ref_db) {

    if (is_mpa(report)) {

        colname_db_minimisers_taxon <- COLNAME_MPA_DB_MINIMISERS_TAXON
        colname_db_minimisers_clade <- COLNAME_MPA_DB_MINIMISERS_CLADE
        colname_ncbi_id <- COLNAME_MPA_NCBI_ID

    } else {

        colname_db_minimisers_taxon <- COLNAME_STD_DB_MINIMISERS_TAXON
        colname_db_minimisers_clade <- COLNAME_STD_DB_MINIMISERS_CLADE
        colname_ncbi_id <- COLNAME_STD_NCBI_ID

    }

    report[, colname_db_minimisers_taxon] <- ref_db[, COLNAME_REF_DB_MINIMISERS_TAXON][match(
        report[, colname_ncbi_id], ref_db[, COLNAME_REF_DB_NCBI_ID]
    )]
    report[, colname_db_minimisers_clade] <- ref_db[, COLNAME_REF_DB_MINIMISERS_CLADE][match(
        report[, colname_ncbi_id], ref_db[, COLNAME_REF_DB_NCBI_ID]
    )]

    return(report)
}


#################################################################################
## UTILITY FUNCTIONS FOR TRANSFERRING COLUMNS BETWEEN DIFFERENT REPORT FORMATS ##
#######################################################################################################

transfer_nReads <- function(report_mpa, report_std) {

    if (is_mpa(report_std)) {
        stop(paste0(
            "The report you provided as 'report_std' is not actually in standard format. ",
            "Please review your input!"
        ))
    } else if (!(is_mpa(report_mpa))) {
        stop(paste0(
            "The report you provided as 'report_mpa' is not actually in MPA format. ",
            "Please review your input!"
        ))
    }

    report_mpa[, COLNAME_MPA_TOTAL_READS] <- report_std[, COLNAME_STD_TOTAL_READS][match(
        report_mpa[, COLNAME_MPA_SAMPLE], report_std[, COLNAME_STD_SAMPLE]
    )]
    report_mpa[, COLNAME_MPA_TOTAL_READS] <- as.numeric(report_mpa[, COLNAME_MPA_TOTAL_READS])

    return(report_mpa)
}


transfer_ncbiID <- function(report_mpa, report_std) {

    if (is_mpa(report_std)) {
        stop(paste0(
            "The report you provided as 'report_std' is not actually in standard format. ",
            "Please review your input!"
        ))
    } else if (!(is_mpa(report_mpa))) {
        stop(paste0(
            "The report you provided as 'report_mpa' is not actually in MPA format. ",
            "Please review your input!"
        ))
    }

    ranks <- unique(report_std[, COLNAME_STD_RANK])
    ranks <- ranks[ranks %in% report_mpa[, COLNAME_MPA_RANK]]
    ranks <- get_association(ranks)

    report_mpa[["temp_column"]] <- seq_len(nrow(report_mpa))

    final_report <- data.frame(matrix(nrow = 0, ncol = (ncol(report_mpa) + 1)))

    for (rank in ranks) {

        subset <- report_mpa[report_mpa[, COLNAME_MPA_RANK] == rank, ]

        subset[, COLNAME_MPA_NCBI_ID] <- report_std[, COLNAME_STD_NCBI_ID][match(
            subset[, names(ranks)[ranks == rank]], report_std[, COLNAME_STD_TAXON]
        )]

        final_report <- rbind(final_report, subset)

    }

    final_report <- final_report[order(final_report[["temp_column"]], decreasing = FALSE), ]
    final_report[["temp_column"]] <- NULL

    return(final_report)
}

transferDomains <- function(report_std, report_mpa, verbose = TRUE, inference = TRUE) {

    if (is_mpa(report_std)) {
        stop(paste0(
            "The report you provided as 'report_std' is not actually in standard format. ",
            "Please review your input!"
        ))
    } else if (!(is_mpa(report_mpa))) {
        stop(paste0(
            "The report you provided as 'report_mpa' is not actually in MPA format. ",
            "Please review your input!"
        ))
    }

    if (verbose) {
        cat("Adding domain info from the MPA-style report to the standard report\n")
    }

    report_std <- retrieve_rankDomains(report_std, report_mpa, verbose)
    #report_std <- retrieve_subrankDomains(report_std, report_mpa, verbose, inference)
    
    return(report_std)
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
        
        stop(paste0(
            "This function does not support MPA-style reports. Please provide a standard report."
        ))

    } else {

        colname_uniq_minimisers <- COLNAME_STD_UNIQ_MINIMISERS
        colname_db_minimisers_clade <- COLNAME_STD_DB_MINIMISERS_CLADE
        colname_db_minimisers_taxon <- COLNAME_STD_DB_MINIMISERS_TAXON
        colname_ratio_clade <- COLNAME_STD_RATIO_CLADE
        colname_ratio_taxon <- COLNAME_STD_RATIO_TAXON
    
        report[, colname_ratio_taxon] <- (report[, colname_uniq_minimisers] / report[, colname_db_minimisers_taxon])
        report[, colname_ratio_clade] <- (report[, colname_uniq_minimisers] / report[, colname_db_minimisers_clade])
    }

    return(report)
}


#' ASSESS THE STATISTICAL SIGNIFICANCE OF RESULTS
#' 
#' @param report
#' @param ref_db 
#' @return An updated version of the input dataframe, with new columns containing statistical significance results.
assess_statSig <- function(report, ref_db, verbose = TRUE) {

    if (is_mpa(report)) {
        stop(paste0("This function does not support MPA-style reports. Please provide a standard report."))
    } 

    p_clade_in_db <- numeric(nrow(report))
    pval <- numeric(nrow(report))

    if (verbose) {
        cat("Assessing the statistical significance of results at the family, genus and species levels\n")
        pb <- txtProgressBar(min = 0, max = nrow(report), style = 3)
    }

    for (i in seq_len(nrow(report))) {

        if (verbose) setTxtProgressBar(pb, i)

        # Get proportion of clade-level minimisers of a given taxon in the reference database (DB).
        p_clade_in_db_ <- (report[, COLNAME_STD_DB_MINIMISERS_CLADE][i] / sum(ref_db[, COLNAME_REF_DB_MINIMISERS_CLADE]))

        # Mean is the total number of reads analysed from the sample (sample size) times the probability of success which 
        # is equal to the proportion of clade-level minimisers of the taxon out of the total available in the reference DB.
        mean_ <- (report[, COLNAME_STD_TOTAL_READS][i] * p_clade_in_db_)

        # Standard deviation = sqrt(n*P*(1-p))
        sdev <- sqrt(report[, COLNAME_STD_TOTAL_READS][i] * p_clade_in_db_ * (1 - p_clade_in_db_)) 

        # Get p-values.
        pval_ <- pnorm(
            q = report[, COLNAME_STD_UNIQ_MINIMISERS][i], 
            mean = mean_, 
            sd = sdev, 
            lower.tail = FALSE
        )
            
        # Saving probability of success.
        p_clade_in_db[i] <- p_clade_in_db_

        # Saving p-value.
        pval[i] <- pval_
    }

    report[, COLNAME_STD_P_CLADE_IN_DB] <- p_clade_in_db
    report[, COLNAME_STD_PVALUE] <- pval

    # Create vector to store adjusted p-values.
    padj <- rep(0, times = nrow(report))

    if (verbose) {
        cat("\nCorrecting p-values\n")
        pb <- txtProgressBar(min = 0, max = nrow(report), style = 3)
    }

    # Iterate over dataframe rows...
    for (i in seq_len(nrow(report))){

        if (verbose) setTxtProgressBar(pb, i)
        
        # Identify how many families, genuses or species were identified in a given sample.
        n_taxa_rank <- NA

        if (report[, COLNAME_STD_RANK][i] %in% c("F", "G", "S")) {
            n_taxa_rank <- get_nTaxaInRank(
            report = report, 
                rank = report[, COLNAME_STD_RANK][i], 
                sample = report[, COLNAME_STD_SAMPLE][i]
            )
        } else {
            stop(paste0("Invalid rank value in row: ", i))
        }

        # Perform p-value correction.
        padj[i] <- p.adjust(report[, COLNAME_STD_PVALUE][i], method = "BH", n = n_taxa_rank)
    }

    # Add adjusted p-values to dataframe.
    report[, COLNAME_STD_PADJ] <- padj

    # Add column stating whether results are significant or not based on the adjusted p-value.
    report[, COLNAME_STD_SIGNIF] <- ifelse(
        report[, COLNAME_STD_PADJ] <= 0.05,
        "Significant",
        "Non-significant"
    )

    cat("\n")

    return(report)
}
