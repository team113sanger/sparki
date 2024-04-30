assess_n_taxa_per_sample <- function(merged_reports) {

    n_taxa_per_sample <- table(merged_reports$sample)
    merged_reports$n_taxa_per_sample <- n_taxa_per_sample[merged_reports$sample]

    return(merged_reports)
}

#' ASSESS THE TOTAL NUMBER OF READS ANALYSED FOR EACH SAMPLE
#' 
#' @param merged_reports
#' @return An updated version of the input dataframe, with a new column containing sample sizes.
assess_n_reads_per_sample <- function(merged_reports) {

    # Get numbers of classified/unclassified reads for each sample.
    merged_reports$sample_size <- sapply(1:nrow(merged_reports), function(x) {

        sample <- merged_reports$sample[x]
        n_unclassified_reads <- as.numeric(merged_reports$number_fragments_clade[(merged_reports$scientific_name == "unclassified") & (merged_reports$sample == sample)])
        n_classified_reads <- as.numeric(merged_reports$number_fragments_clade[(merged_reports$scientific_name == "root") & (merged_reports$sample == sample)])

        sample_size <- (n_unclassified_reads + n_classified_reads)

        return(sample_size)
    })

    return(merged_reports)
}

#' ASSESS THE STATISTICAL SIGNIFICANCE OF RESULTS
#' 
#' @param merged_reports
#' @param reference_db 
#' @return An updated version of the input dataframe, with new columns containing statistical significance results.
assess_statistical_significance <- function(merged_reports, reference_db) {

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
            n_taxons_rank <- merged_reports$nFamilies_sample[i]  
        } else if (merged_reports$rank[i] == "G") {
            n_taxons_rank <- merged_reports$nGenuses_sample[i] 
        } else if (merged_reports$rank[i] == "S") {
            n_taxons_rank <- merged_reports$nSpecies_sample[i] 
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
