#######################################################
## UTILITY FUNCTION FOR LOADING A REFERENCE DATABASE ##
#######################################################################################################

#' LOAD INFORMATION FROM KRAKEN2'S REFERENCE DATABASE
#' 
#' This function takes a path to an inspect.txt file inside a Kraken2 reference database
#' and reads the table. The first 7 lines are skipped, as they correspond to header lines.
#' Column names are added before the reference database dataframe is returned.
#' 
#' @param reference_path Path to an inspect.txt file inside a Kraken2 reference database.
#' @return 
#' 
loadReference <- function(reference_path) {

    # Read Kraken2 reference file (inspect.txt).
    ref <- read.table(
        reference_path,
        row.names = NULL, 
        header = FALSE,
        sep = "\t",
        skip = 7, # Skip header lines.
        quote = "",
        stringsAsFactors = FALSE
    )
    # Add column names.
    colnames(ref) <- c(
        COLNAME_REF_DB_PCT_FRAG_CLADE, 
        COLNAME_REF_DB_MINIMISERS_CLADE,
        COLNAME_REF_DB_MINIMISERS_TAXON,
        COLNAME_REF_DB_RANK, 
        COLNAME_REF_DB_NCBI_ID, 
        COLNAME_REF_DB_TAXON
    )

    ref[, COLNAME_REF_DB_TAXON] <- gsub("^[[:space:]]\\s*(.*?)", "", ref[, COLNAME_REF_DB_TAXON], perl = TRUE)

    return(ref)
}

#############################################
## UTILITY FUNCTIONS FOR HANDLING METADATA ##
#######################################################################################################

#' LOAD METADATA TABLE
#' 
#' This function takes a path to a metadata file and reads the table.
#' and reads the table. The first 7 lines are skipped, as they correspond to header lines.
#' Column names are added before the reference database dataframe is returned.
#' 
#' @param mdata_path Path to a metadata table.
#' @return 
#' 
loadMetadata <- function(mdata_path) {

    # Read metadata file.
    mdata <- readr::read_delim(mdata_path)

    return(mdata)
}

#' ADD COLUMNS FROM ONE DATAFRAME TO ANOTHER
#' 
#' This function takes as input a dataframe (df1) and adds metadata to it from another dataframe (df2). Both 
#' dataframes must have the same sample IDs present, as these IDs will be used as an "anchor" between the tables. 
#' Sample IDs can be specified as either row names or a standard column in the dataframes. The columns from df2
#' that should be added to df1 are specified by the user.
#' 
#' @param df1 Dataframe 1.
#' @param df1_sample_col Name of column in df1 that contains sample IDs; alternatively, this can be "rownames".
#' @param df2 Dataframe 2.
#' @param df2_sample_col Name of column in df2 that contains sample IDs; alternatively, this can be "rownames".
#' @param categories Categories of interest from df2 that should be added to df1.
#' @return An updated version of df1.
#' @export
addMetadata <- function(report, metadata, metadata_sample_col, metadata_columns) {

    # Check report format.
    report_colname_sample <- ifelse(is_mpa(report), COLNAME_MPA_SAMPLE, COLNAME_STD_SAMPLE)

    # Process metadata dataframe (tibble).
    metadata <- metadata |> 
        # Select metadata columns that will be added to the report.
        dplyr::select(metadata_sample_col, metadata_columns) |>
        # Rename column with sample IDs to keep its name consistent
        # with the report.
        dplyr::rename(!!report_colname_sample := metadata_sample_col)
    
    # Replace spaces (if any) with underscores.
    colnames(metadata) <- gsub(" ", "_", colnames(metadata))
    metadata_columns <- gsub(" ", "_", metadata_columns)

    # Add metadata columns to the report.
    report <- dplyr::left_join(
        report, metadata,
        by = dplyr::join_by(!!report_colname_sample)
    ) |>
        # Get names of columns that contain results (and not sample names / metadata) and
        # reorder the columns so that the metadata stays between the sample IDs and the
        # Kraken2 results.
        dplyr::relocate(metadata_columns, .after = !!report_colname_sample)

    return(report)
}

##################################################################
## HELPER FUNCTIONS FOR READING, MERGING AND SUBSETTING REPORTS ##
#######################################################################################################

#' LOAD MPA-STYLE REPORTS 
#' 
#' This function takes a path to a directory containing MPA-style reports, reading and processing
#' all reports into a single dataframe.
#' 
#' @param mpa_reports_dir Path to a directory containing MPA-style reports.
#' @param verbose Whether to output log messages.
#' @return A dataframe (tibble) with the content of all MPA-style reports from the specified directory.
#' @export
#'
load_MPAreports <- function(mpa_reports_dir, verbose = TRUE) {

    # Get paths to MPA-style reports in a specified directory.
    mpa_files <- fs::dir_ls(mpa_reports_dir, glob = "*.mpa$")

    # Check if directory really has any MPA-style reports...
    if (length(mpa_files) == 0) {
        stop(paste0("No MPA-style reports were found at ", mpa_reports_dir, ". Please review your input."))
    } 

    # Create vector with taxonomic rank names.
    taxonomy <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")

    # Create a dataframe (tibble).
    mpa_reports <- readr::read_tsv(
        mpa_files, 
        col_names = c(COLNAME_MPA_TAXON_HIERARCHY, COLNAME_MPA_N_FRAG_CLADE), 
        id = get("COLNAME_MPA_SAMPLE")
    ) 

    # Add rank column.
    mpa_reports[, COLNAME_MPA_RANK] <- lapply(mpa_reports[, COLNAME_MPA_TAXON_HIERARCHY], addRank)

    # Further process the dataframe.
    mpa_reports <- mpa_reports |>
        # Split rows by "|" into primary taxa.
        tidyr::separate(col = COLNAME_MPA_TAXON_HIERARCHY, into = taxonomy, sep = "\\|") |>
        # Cleanup the names in the taxonomy columns.
        dplyr::mutate(dplyr::across(taxonomy, ~ stringr::str_remove(.x, pattern = "[a-z]__"))) |>
        # Simplify sample IDs.
        dplyr::mutate(!!COLNAME_MPA_SAMPLE := stringr::str_remove(basename(sample), ".kraken.mpa")) |>
        # Collect the rightmost non-NA item in each row.
        dplyr::mutate(
            !!COLNAME_MPA_TAXON_LEAF := dplyr::coalesce(species, genus, family, order, class, phylum, kingdom, domain),
            .before = "domain"
        ) |>
        # Replace underscores with spaces in taxon names.
        dplyr::mutate(
            !!COLNAME_MPA_TAXON_LEAF := stringr::str_replace_all(taxon_leaf, pattern = "_", replacement = " ")
        ) |>
        # Put rank column right after the taxon leaf column.
        dplyr::relocate(!!COLNAME_MPA_RANK, !!COLNAME_MPA_N_FRAG_CLADE, .after = !!COLNAME_MPA_TAXON_LEAF)

    if (verbose) cat("MPA-style reports loaded successfully.\n")
    
    return(mpa_reports) 
}

#' LOAD STANDARD REPORTS 
#' 
#' This function takes a path to a directory containing standard reports, reading and processing
#' all reports into a single dataframe.
#' 
#' @param mpa_reports_dir Path to a directory containing standard reports.
#' @param verbose Whether to output log messages.
#' @return A dataframe (tibble) with the content of all standard reports from the specified directory.
#' @export
#'
load_STDreports <- function(std_reports_dir, verbose = TRUE) {

    # Get paths to standard reports in a specified directory.
    std_files <- fs::dir_ls(std_reports_dir, glob = "*.kraken$")

    # Check if directory really has any standard reports...
    if (length(std_files) == 0) {
        stop(paste0("No standard reports were found at ", std_reports_dir, ". Please review your input."))
    }
    # Create a dataframe (tibble) and process.
    std_reports <- readr::read_tsv(
        std_files,
        col_names = c(COLNAME_STD_PCT_FRAG_CLADE, COLNAME_STD_N_FRAG_CLADE, COLNAME_STD_N_FRAG_TAXON, 
            COLNAME_STD_MINIMISERS, COLNAME_STD_UNIQ_MINIMISERS, COLNAME_STD_RANK, COLNAME_STD_NCBI_ID, 
            COLNAME_STD_TAXON),
        id = get("COLNAME_STD_SAMPLE")
    ) |>
        # Remove the subranks.
        dplyr::filter(!grepl(!!as.name(COLNAME_STD_RANK), pattern = "[0-9]")) |>
        # Simplify sample IDs.
        dplyr::mutate(!!COLNAME_STD_SAMPLE := stringr::str_remove(basename(sample), ".kraken"))

    if (verbose) cat("Standard reports loaded successfully.\n")

    return(std_reports) 
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
mergeReports <- function(std_reports, mpa_reports) {

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

    return(merged_reports)
}

subsetReports <- function(report, include_human = FALSE) {

    # Select relevant rows.
    report <- report |> dplyr::filter(!!as.name(COLNAME_STD_RANK) %in% c("F", "G", "S"))
    
    # Filter out human results.    
    if (include_human == FALSE) {
        report <- report |> 
            dplyr::filter(!(!!as.name(COLNAME_STD_TAXON) %in% c("Hominidae", "Homo", "Homo sapiens")))
    }

    return(report)
}

#########################################################
## UTILITY FUNCTIONS FOR CREATING AUXILIARY DATAFRAMES ##
#######################################################################################################

# Number of reads under each domain for each sample.
get_nDomainReads <- function(report, include_eukaryotes) {

    colname_sample <- ifelse(is_mpa(report), COLNAME_MPA_SAMPLE, COLNAME_STD_SAMPLE)
    colname_taxon <- ifelse(is_mpa(report), COLNAME_MPA_TAXON_LEAF, COLNAME_STD_TAXON)
    colname_n_frag_clade <- ifelse(is_mpa(report), COLNAME_MPA_N_FRAG_CLADE, COLNAME_STD_N_FRAG_CLADE)

    domains <- "Viruses|Bacteria|Archaea"
    if (include_eukaryotes == TRUE) domains <- paste0("Eukaryota|", domains)

    domain_reads <- report |>
        # Select relevant columns.
        dplyr::select(!!colname_sample, !!colname_taxon, !!colname_n_frag_clade) |>
        # Select rows that contain the relevant taxa.
        dplyr::filter(grepl(!!as.name(colname_taxon), pattern = domains)) |>
        # Create grouping for sample + taxon.
        dplyr::group_by(!!as.name(colname_sample), !!as.name(colname_taxon)) |>
        # Sum classified reads for each sample + taxon group.
        dplyr::summarise(
            !!COLNAME_DOMAIN_READS_N_FRAG := sum(!!as.name(colname_n_frag_clade)),
            .groups = "keep"
        ) |>
        # Rename columns.
        dplyr::rename(
            !!COLNAME_DOMAIN_READS_TAXON := !!colname_taxon,
            !!COLNAME_DOMAIN_READS_SAMPLE := !!colname_sample
        )

    return(domain_reads)
}

getClassificationSummary <- function(report) {

    if (is_mpa(report)) {
        stop(paste0(
            "This function is not applicable to MPA-style reports. The purpose of this function ",
            "is to assess the number of classified and unclassified reads for each sample, but the ",
            "MPA-style reports do not contain information on unclassified reads."
        ))
    }

    summary <- report |> 
        # Select columns of interest.
        dplyr::select(!!COLNAME_STD_SAMPLE, !!COLNAME_STD_TAXON, !!COLNAME_STD_N_FRAG_CLADE) |>
        # Select rows of interest.
        dplyr::filter(!!as.name(COLNAME_STD_TAXON) %in% c("unclassified", "root")) |>
        # Rename column.
        dplyr::rename(
            !!COLNAME_CLASSIF_SUMMARY_SAMPLE := !!COLNAME_STD_SAMPLE,
            !!COLNAME_CLASSIF_SUMMARY_READ_TYPE := !!COLNAME_STD_TAXON,
            !!COLNAME_CLASSIF_SUMMARY_N_FRAG := !!COLNAME_STD_N_FRAG_CLADE
        ) |>
        # Rename elements in column.
        dplyr::mutate(!!COLNAME_CLASSIF_SUMMARY_READ_TYPE := dplyr::case_when(
            !!as.name(COLNAME_CLASSIF_SUMMARY_READ_TYPE) == "unclassified" ~ "Unclassified",
            !!as.name(COLNAME_CLASSIF_SUMMARY_READ_TYPE) == "root" ~ "Classified"
        )) |>
        # Create column with log-transformed number of clade-level fragments.
        dplyr::mutate(
            !!COLNAME_CLASSIF_SUMMARY_LOG_N_FRAG := log10(!!as.name(COLNAME_CLASSIF_SUMMARY_N_FRAG))
        ) |>
        # Reorder columns.
        dplyr::relocate(
            !!COLNAME_CLASSIF_SUMMARY_LOG_N_FRAG, .after = !!COLNAME_CLASSIF_SUMMARY_N_FRAG
        )
        
    return(summary)
}

getClassificationProportion <- function(report, taxon) {

    if (is_mpa(report)) {
        stop(paste0(
            "This function is not applicable to MPA-style reports. The purpose of this function ",
            "is to assess the proportion of classified reads relative to the total number of reads ",
            "assessed, but the MPA-style reports do not contain information on unclassified reads."
        ))
    }

    proportion <- report |>
        # Select relevant columns.
        dplyr::select(
            !!COLNAME_STD_SAMPLE,
            !!COLNAME_STD_SAMPLE_SIZE,
            !!COLNAME_STD_TAXON,
            !!COLNAME_STD_N_FRAG_CLADE
        ) |>
        # Select relevant rows.
        dplyr::filter(!!as.name(COLNAME_STD_TAXON) == "root") |>
        # Create column with proportion of classified reads.
        dplyr::mutate(
            !!COLNAME_PROP_CLASSIFIED := (!!as.name(COLNAME_STD_N_FRAG_CLADE) / !!as.name(COLNAME_STD_SAMPLE_SIZE))
        ) |>
        # Remove columns.
        dplyr::select(!c(!!COLNAME_STD_SAMPLE_SIZE, !!COLNAME_STD_TAXON, !!COLNAME_STD_N_FRAG_CLADE)) |>
        # Rename columns.
        dplyr::rename(!!COLNAME_PROP_SAMPLE := !!COLNAME_STD_SAMPLE)

    return(proportion)
}

getSignificanceSummary <- function(report) {

    if (is_mpa(report)) stop(paste0("This function does not support MPA-style reports."))

    summary <- report |>
        # Select relevant columns.
        dplyr::select(
            !!COLNAME_STD_SAMPLE,
            !!COLNAME_STD_RANK,
            !!COLNAME_STD_TAXON,
            !!COLNAME_STD_SIGNIF
        ) |>
        # Create grouping based on sample ID, rank and statistical significance.
        dplyr::group_by(
            !!as.name(COLNAME_STD_SAMPLE),
            !!as.name(COLNAME_STD_RANK),
            !!as.name(COLNAME_STD_SIGNIF)
        ) |>
        # Count number of taxa per group.
        dplyr::summarise(
            !!COLNAME_SIGNIF_SUMMARY_N_TAXA := length(!!as.name(COLNAME_STD_TAXON)),
            .groups = "keep"
        ) |>
        # Rename columns.
        dplyr::rename(
            !!COLNAME_SIGNIF_SUMMARY_SAMPLE := !!COLNAME_STD_SAMPLE,
            !!COLNAME_SIGNIF_SUMMARY_RANK := !!COLNAME_STD_RANK,
            !!COLNAME_SIGNIF_SUMMARY_SIGNIF := !!COLNAME_STD_SIGNIF
        ) |>
        # Replace values in column.
        dplyr::mutate(!!COLNAME_SIGNIF_SUMMARY_RANK := dplyr::case_when(
            !!as.name(COLNAME_SIGNIF_SUMMARY_RANK) == "F" ~ "Family",
            !!as.name(COLNAME_SIGNIF_SUMMARY_RANK) == "G" ~ "Genus",
            !!as.name(COLNAME_SIGNIF_SUMMARY_RANK) == "S" ~ "Species"
        ))

    return(summary)

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

#############################################################################
## UTILITY FUNCTIONS FOR ASSESSING THE STATISTICAL SIGNIFICANCE OF RESULTS ##
#######################################################################################################

assessMinimiserRatio <- function(report) {

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
