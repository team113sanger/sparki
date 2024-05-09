addRank <- function(report) {

    if (is_mpa(report)) {

        report[, COLNAME_RANK_MPA] <- sapply(seq_len(nrow(report)), function(x) {

            # Get scientific name (in this case it will have a prefix indicating the taxonomic rank).
            rank <- tail(unlist(stringr::str_split(report[, COLNAME_TAXON_MPA][x], "\\|")), n = 1)

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

    } else {
        stop(paste0(
            "This function is not applicable to a standard report (i.e. not MPA-style). ",
            "The standard report format already has the column ", COLNAME_RANK_STD, ", which indicates the taxonomic rank of each line."
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

            if (verbose) cat("Adding column with taxon names for rank:", names(ranks)[ranks == rank], "\n")

            # Get name for column that will be added to the MPA-style report.
            colname <- names(ranks)[ranks == rank]

            # Create column with concise taxon name at a given rank.
            # Examples:
            # "d__Eukaryota" -> "Eukaryota" (under "domain")
            # "(...)f__Papillomaviridae|g__Alphapapillomavirus" -> "Papillomaviridae" (under "family") and "Alphapapillomavirus" (under "genus")
            # "(...)g__Homo|s__Homo sapiens" -> "Homo" (under "genus") and "Homo sapiens" (under "species")
            report[, colname] <- sapply(seq_len(nrow(report)), function(x) {

                if (report[, COLNAME_RANK_MPA][x] == rank) {

                    return(extract_taxon(report[, COLNAME_TAXON_MPA][x], rank = rank, last_in_hierarchy = TRUE))

                } else {

                    return(extract_taxon(report[, COLNAME_TAXON_MPA][x], rank = rank, last_in_hierarchy = FALSE))

                }
            })
        }
    } else {
        stop(paste0(
            "This function is not applicable to a standard report (i.e. not MPA-style). The standard report ",
            "format already has the column ", COLNAME_TAXON_STD, ", which has concise taxon names instead of ",
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
            "is to assess the number of classified and unclassified reads for each sample, but the ",
            "MPA-style reports do not contain information on unclassified reads."
        ))
    } else {

        report$n_total_reads_in_sample <- NA

        # Get numbers of classified/unclassified reads for each sample.
        subset <- report[report$scientific_name %in% c("unclassified", "root"), ]

        for (sample in unique(subset$sample)) {

            subset_sample <- subset[subset$sample == sample, ]
            n_unclassified_reads <- as.numeric(subset_sample$number_fragments_clade[subset_sample$scientific_name == "unclassified"])
            n_classified_reads <- as.numeric(subset_sample$number_fragments_clade[subset_sample$scientific_name == "root"])

            report$n_total_reads_in_sample[report$sample == sample] <- (n_unclassified_reads + n_classified_reads)
        }
    }

    return(report)
}

add_ncbiID <- function(report_std, report_mpa) {

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

    ranks <- unique(report_std[, COLNAME_RANK_STD])
    ranks <- ranks[ranks %in% report_mpa[, COLNAME_RANK_MPA]]
    ranks <- get_association(ranks)

    for (rank in ranks) {

        report_mpa[, COLNAME_NCBI_ID_MPA] <- report_std[, COLNAME_NCBI_ID_STD][match(
            report_mpa[, names(ranks)[ranks == rank]], report_std[, COLNAME_TAXON_STD]
        )]

    }

    return(report_mpa)
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
        COLNAME_SAMPLE_DOMAIN_READS, 
        COLNAME_SCOPE_DOMAIN_READS, 
        COLNAME_N_READS_DOMAIN_READS
    )
    
    n_domainReads[, COLNAME_N_READS_DOMAIN_READS] <- as.numeric(n_domainReads[, COLNAME_N_READS_DOMAIN_READS])

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
        colname_sample <- COLNAME_SAMPLE_MPA
        colname_taxon <- COLNAME_TAXON_MPA
        colname_rank <- COLNAME_RANK_MPA
    } else {
        colname_sample <- COLNAME_SAMPLE_STD
        colname_taxon <- COLNAME_TAXON_STD
        colname_rank <- COLNAME_RANK_STD
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



add_DBinfo <- function(report, ref_db) {

    if (is_mpa(report)) {

        colname_db_minimisers_taxon <- COLNAME_DB_MINIMISERS_TAXON_MPA
        colname_db_minimisers_clade <- COLNAME_DB_MINIMISERS_CLADE_MPA
        colname_ncbi_id <- COLNAME_NCBI_ID_MPA

    } else {

        colname_db_minimisers_taxon <- COLNAME_DB_MINIMISERS_TAXON_STD
        colname_db_minimisers_clade <- COLNAME_DB_MINIMISERS_CLADE_STD
        colname_ncbi_id <- COLNAME_NCBI_ID_STD

    }

    report[, colname_db_minimisers_taxon] <- ref_db[, COLNAME_MINIMISERS_TAXON_REF_DB][match(report[, colname_ncbi_id], ref_db[, COLNAME_NCBI_ID_REF_DB])]
    report[, colname_db_minimisers_clade] <- ref_db[, COLNAME_MINIMISERS_CLADE_REF_DB][match(report[, colname_ncbi_id], ref_db[, COLNAME_NCBI_ID_REF_DB])]

    return(report)
}

getClassificationSummary <- function(report) {

    if (is_mpa(report)) {
        stop(paste0(
            "This function is not applicable to MPA-style reports. The purpose of this function ",
            "is to assess the number of classified and unclassified reads for each sample, but the ",
            "MPA-style reports do not contain information on unclassified reads."
        ))
    } else {
        colname_n_frag_clade <- COLNAME_N_FRAG_CLADE_STD
        colname_taxon <- COLNAME_TAXON_STD
        colname_sample <- COLNAME_SAMPLE_STD
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
        COLNAME_SAMPLE_CLASSIF_SUMMARY, 
        COLNAME_READ_TYPE_CLASSIF_SUMMARY, 
        COLNAME_N_READS_CLASSIF_SUMMARY
    )

    class_unclass_df[, COLNAME_N_READS_CLASSIF_SUMMARY] <- as.numeric(class_unclass_df[, COLNAME_N_READS_CLASSIF_SUMMARY])

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
        colname_n_frag_clade <- COLNAME_N_FRAG_CLADE_STD
        colname_taxon <- COLNAME_TAXON_STD
        colname_sample <- COLNAME_SAMPLE_STD
    }

    class_unclass_df <- getClassificationSummary(report)

    for (sample in unique(class_unclass_df[, COLNAME_SAMPLE_CLASSIF_SUMMARY])) {
        subset <- class_unclass_df[, class_unclass_df[, COLNAME_SAMPLE_CLASSIF_SUMMARY] == sample, ]
        
        classified_reads <- subset[, COLNAME_N_READS_CLASSIF_SUMMARY][subset[, COLNAME_READ_TYPE_CLASSIF_SUMMARY] == "Classified"]
        unclassified_reads <- subset[, COLNAME_N_READS_CLASSIF_SUMMARY][subset[, COLNAME_READ_TYPE_CLASSIF_SUMMARY] == "Unclassified"]

        total_reads <- classified_reads + unclassified_reads

    }



    return(class_unclass_df)
}

addDomain <- function(std_report, mpa_report) {

    ranks <- c("D", "K", "P", "C", "O", "F", "G", "S")

    report$domain <- sapply(seq_len(nrow(report)), function(x) {





        for (rank in ranks) {
            
            if (report[, COLNAME_RANK_STD][x] == rank) {
                return()
            }
        }

    })


    merged_reports$domain <- sapply(seq_len(nrow(merged_reports)), function(x) {

        domain <- NA

        if (merged_reports$rank[x] == "F") {
            domain <- unique(merged_mpa$domain[which(merged_mpa$family == merged_reports$scientific_name[x])])
        } else if (merged_reports$rank[x] == "G") {
            domain <- unique(merged_mpa$domain[which(merged_mpa$genus == merged_reports$scientific_name[x])])
        } else if (merged_reports$rank[x] == "S") {
            domain <- unique(merged_mpa$domain[which(merged_mpa$species == merged_reports$scientific_name[x])])
        }

        return(domain)
    })

    return(merged_reports)
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
