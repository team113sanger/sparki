addRank <- function(report, verbose = TRUE) {

    if (is_mpa(report)) {

        if (verbose) {
            cat("Adding", COLNAME_MPA_RANK, "column to the MPA-style report\n")
            pb <- txtProgressBar(min = 0, max = nrow(report), style = 3)
        }

        report[, COLNAME_MPA_RANK] <- sapply(seq_len(nrow(report)), function(x) {

            if (verbose) setTxtProgressBar(pb, x)

            # Get scientific name (in this case it will have a prefix indicating the taxonomic rank).
            rank <- tail(unlist(stringr::str_split(report[, COLNAME_MPA_TAXON][x], "\\|")), n = 1)

            # Identify which taxonomic rank the scientific name belongs to.
            if(grepl("s__", rank)) return("S")
            if(grepl("g__", rank)) return("G")
            if(grepl("f__", rank)) return("F")
            if(grepl("o__", rank)) return("O")
            if(grepl("c__", rank)) return("C")
            if(grepl("p__", rank)) return("P")
            if(grepl("k__", rank)) return("K")
            if(grepl("d__", rank)) return("D")
        })
        cat("\n")

    } else {
        stop(paste0(
            "This function is not applicable to a standard report (i.e. not MPA-style). ",
            "The standard report format already has the column ", COLNAME_STD_RANK, ", which indicates the taxonomic rank of each line."
        ))   
    }

    return(report)
}

addConciseTaxon <- function(report, verbose = TRUE) {

    # If report is in MPA format...
    if (is_mpa(report)) {

        ranks <- c("D", "K", "P", "C", "O", "F", "G", "S")
        ranks <- get_association(ranks)

        # Iterate over ranks... 
        for (rank in ranks) {

            if (verbose) {
                cat("Adding column with taxon names for rank:", names(ranks)[ranks == rank], "\n")
                pb <- txtProgressBar(min = 0, max = nrow(report), style = 3)
            }

            # Get name for column that will be added to the MPA-style report.
            colname <- names(ranks)[ranks == rank]

            # Create column with concise taxon name at a given rank.
            # Examples:
            # "d__Eukaryota" -> "Eukaryota" (under "domain")
            # "(...)f__Papillomaviridae|g__Alphapapillomavirus" -> "Papillomaviridae" (under "family") and "Alphapapillomavirus" (under "genus")
            # "(...)g__Homo|s__Homo sapiens" -> "Homo" (under "genus") and "Homo sapiens" (under "species")
            report[, colname] <- sapply(seq_len(nrow(report)), function(x) {

                if (verbose) setTxtProgressBar(pb, x)

                if (report[, COLNAME_MPA_RANK][x] == rank) {

                    return(extract_taxon(report[, COLNAME_MPA_TAXON][x], rank = rank, last_in_hierarchy = TRUE))

                } else {

                    return(extract_taxon(report[, COLNAME_MPA_TAXON][x], rank = rank, last_in_hierarchy = FALSE))

                }
            })

            cat("\n")

        }

    } else {
        stop(paste0(
            "This function is not applicable to a standard report (i.e. not MPA-style). The standard report ",
            "format already has the column ", COLNAME_STD_TAXON, ", which has concise taxon names instead of ",
            "taxon hierarchies."
        ))   
    }

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
    } else {

        report[, COLNAME_STD_TOTAL_READS] <- NA

        # Get numbers of classified/unclassified reads for each sample.
        subset <- report[report[, COLNAME_STD_TAXON] %in% c("unclassified", "root"), ]

        for (sample in unique(subset[, COLNAME_STD_SAMPLE])) {

            subset_sample <- subset[subset[, COLNAME_STD_SAMPLE] == sample, ]
            n_unclassified_reads <- as.numeric(subset_sample[, COLNAME_STD_N_FRAG_CLADE][subset_sample[, COLNAME_STD_TAXON] == "unclassified"])
            n_classified_reads <- as.numeric(subset_sample[, COLNAME_STD_N_FRAG_CLADE][subset_sample[, COLNAME_STD_TAXON] == "root"])

            report[, COLNAME_STD_TOTAL_READS][report[, COLNAME_STD_SAMPLE] == sample] <- (n_unclassified_reads + n_classified_reads)
        }
    }

    return(report)
}

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

    report_mpa[, COLNAME_MPA_TOTAL_READS] <- report_std[, COLNAME_STD_TOTAL_READS][match(report_mpa[, COLNAME_MPA_SAMPLE], report_std[, COLNAME_STD_SAMPLE])]
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

    report_mpa$temp_column <- seq_len(nrow(report_mpa))

    final_report <- data.frame(matrix(nrow = 0, ncol = (ncol(report_mpa) + 1)))

    for (rank in ranks) {

        subset <- report_mpa[report_mpa[, COLNAME_MPA_RANK] == rank, ]

        subset[, COLNAME_MPA_NCBI_ID] <- report_std[, COLNAME_STD_NCBI_ID][match(
            subset[, names(ranks)[ranks == rank]], report_std[, COLNAME_STD_TAXON]
        )]

        final_report <- rbind(final_report, subset)

    }

    final_report <- final_report[order(final_report$temp_column, decreasing = FALSE), ]
    final_report$temp_column <- NULL

    return(final_report)
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

    report[, colname_db_minimisers_taxon] <- ref_db[, COLNAME_REF_DB_MINIMISERS_TAXON][match(report[, colname_ncbi_id], ref_db[, COLNAME_REF_DB_NCBI_ID])]
    report[, colname_db_minimisers_clade] <- ref_db[, COLNAME_REF_DB_MINIMISERS_CLADE][match(report[, colname_ncbi_id], ref_db[, COLNAME_REF_DB_NCBI_ID])]

    return(report)
}

transferDomains <- function(report_std, report_mpa, verbose = TRUE) {

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

    if (verbose) cat("Adding domain info from the MPA-style report to the standard report.\n")

    ranks <- unique(report_std[, COLNAME_STD_RANK])
    ranks <- get_association(ranks)

    pb <- txtProgressBar(min = 0, max = nrow(report_std), style = 3)

    report_std[, COLNAME_STD_DOMAIN] <- sapply(seq_len(nrow(report_std)), function(x) {

        setTxtProgressBar(pb, x)

        rank <- names(ranks)[ranks == report_std[, COLNAME_STD_RANK][x]]

        if (rank %in% colnames(report_mpa)) {

            domain <- unique(report_mpa[, NAME_RANK_DOMAIN][which(
                report_mpa[, rank] == report_std[, COLNAME_STD_TAXON][x]
            )])
            
            return(domain)
        
        } else {
            return(NA)
        }

    })

    return(report_std)
}


get_nDomainReads <- function(report) {

    n_domainReads <- data.frame(matrix(nrow = 0, ncol = 4))

    # Iterate over samples in report...
    for (sample in unique(report$sample)) {

        # Subset report to keep only a given sample.
        report_sample <- report[report$sample == sample,]

        # Identify if report follows a standard or MPA format; there are some slight 
        # differences in how domains are written depending on the format!
        if (is_mpa(report)) {

            # In an MPA-style report, the domain will have a prefix ("d__"). Additionally,
            # since the MPA format shows the hierarchy of each taxon, we also need to select
            # specifically the lines that have only the domain-level information.
            all_domains <- "d__Eukaryota$|d__Viruses$|d__Archaea$|d__Bacteria$"
            subset_domains <- "d__Viruses$|d__Archaea$|d__Bacteria$"
        
        } else {

            # In a standard report, each taxon is shown independently from their hierarchy,
            # therefore it is easier to select the lines that contain domain-level information.
            all_domains <- "Eukaryota|Viruses|Archaea|Bacteria"
            subset_domains <- "Viruses|Archaea|Bacteria"

        }

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
    
    n_domainReads[, COLNAME_DOMAIN_READS_N_READS] <- as.numeric(n_domainReads[, COLNAME_DOMAIN_READS_N_READS])

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
    } else {
        colname_sample <- COLNAME_STD_SAMPLE
        colname_taxon <- COLNAME_STD_TAXON
        colname_rank <- COLNAME_STD_RANK
    }

    n_taxa_in_samples <- data.frame(matrix(nrow = 0, ncol = 4))

    # Iterate over samples...
    for (sample in unique(report[, colname_sample])) {

        # Subset report to keep only a given sample.
        report_sample <- report[report[, colname_sample] == sample,]

        # Iterate over domains...
        for (domain in unique(report$domain)) {

            # Iterate over ranks...
            for (rank in ranks) {

                n_taxa <- length(unique(report_sample[, colname_taxon][report[, colname_rank] == rank & report$domain == domain]))
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

    classifications <- c("unclassified", "root")
    names(classifications) <- c("Unclassified", "Classified")

    for (sample in unique(report[, colname_sample])) {
        for (classification in classifications) {

            value <- as.numeric(report[, colname_n_frag_clade][report[, colname_taxon] == classification & report[, colname_sample] == sample])
            class_unclass_df <- rbind(class_unclass_df, c(sample, names(classifications)[classifications == classification], value))
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


report_assess_n_taxa_per_sample <- function(merged_reports) {

    n_taxa_per_sample <- table(merged_reports$sample)
    merged_reports$n_taxa_per_sample <- n_taxa_per_sample[merged_reports$sample]

    return(merged_reports)
}

#' This function takes a report and assesses the number of taxa available at a specific
#' rank (e.g. family, genus, species etc) for each sample. 
#' 
#' @param merged_reports A dataframe containing Kraken2 standard reports for one or more samples.
#' @param rank A vector containing one or more ranks that should be looked into.
#' @return An updated version of the input merged_reports dataframe.
#' @export
report_assess_n_taxa_specific_rank_per_sample <- function(merged_reports, rank) {


    for (r in rank) {

        

        column_name <- paste0("n", r, "_sample")

        merged_reports[, column_name] <- sapply(seq_len(nrow(merged_reports)), function(x) {
            n_in_sample <- length(
                unique(merged_reports$scientific_name[merged_reports$rank == r & merged_reports$sample == merged_reports$sample[x]])
            )
            return(n_in_sample)
        })
    }

    return(merged_reports)
}




#' ASSESS THE STATISTICAL SIGNIFICANCE OF RESULTS
#' 
#' @param merged_reports
#' @param reference_db 
#' @return An updated version of the input dataframe, with new columns containing statistical significance results.
report_assess_statistical_significance <- function(merged_reports, reference_db) {

    p_clade_in_db <- numeric(nrow(merged_reports))
    pval <- numeric(nrow(merged_reports))

    for (i in 1:nrow(merged_reports)){

        # Get proportion of clade-level minimisers of a given taxon in the reference database.
        p_clade_in_db_ <- (merged_reports$db_number_minimisers_clade[i] / sum(reference_db$number_minimisers_clade))

        # Mean is the total number of reads analysed from the sample times the probability of success which is equal 
        # to the proportion of clade-level minimisers of the taxon out of the total available in the reference database.
        mean_ <- (merged_reports$sample_size[i] * p_clade_in_db_)

        # Standard deviation = sqrt(n*P*(1-p))
        sdev <- sqrt(merged_reports$sample_size[i] * p_clade_in_db_ * (1 - p_clade_in_db_)) 

        # Get p-values.
        pval_ <- pnorm(
            q = merged_reports$number_distinct_minimisers[i], 
            mean = mean_, 
            sd = sdev, 
            lower.tail = FALSE
        )

        if (is.na(pval_)) {
            cat(merged_reports$scientific_name[i], "\n")
            cat(mean_, "\n")
            cat(sdev, "\n")
            cat(merged_reports$number_distinct_minimisers[i], "\n")
            cat("\n")

        } else {
            print("No NA")
        }
        
        # Saving probability of success.
        p_clade_in_db[i] <- p_clade_in_db_
        pval[i] <- pval_
    }

    merged_reports$p_clade_in_db <- p_clade_in_db
    merged_reports$pval <- pval

    # Create vector to store adjusted p-values.
    padj <- rep(0, times = nrow(merged_reports))

    # Iterate over dataframe rows...
    for (i in 1:nrow(merged_reports)){
        
        # Identify how many families, genuses or species were identified in a given sample.
        n_taxons_rank <- NULL
        if (merged_reports$rank[i] == "F") {
            n_taxons_rank <- merged_reports$nF_sample[i]  
        } else if (merged_reports$rank[i] == "G") {
            n_taxons_rank <- merged_reports$nG_sample[i] 
        } else if (merged_reports$rank[i] == "S") {
            n_taxons_rank <- merged_reports$nS_sample[i] 
        } else {
            stop(paste0("Invalid rank value in row: ", i))
        }

        # Perform p-value correction.
        padj[i] <- p.adjust(merged_reports$pval[i], method = "BH", n = n_taxons_rank)
    }

    # Add adjusted p-values to dataframe.
    merged_reports$padj <- padj

    return(merged_reports)
}
