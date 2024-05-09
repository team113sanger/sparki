########################################
## HELPER FUNCTIONS FOR READING FILES ##
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
load_reference <- function(reference_path) {

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
        COLNAME_PCT_FRAG_CLADE_REF_DB, 
        COLNAME_MINIMISERS_CLADE_REF_DB,
        COLNAME_MINIMISERS_TAXON_REF_DB,
        COLNAME_RANK_REF_DB, 
        COLNAME_NCBI_ID_REF_DB, 
        COLNAME_TAXON_REF_DB
    )

    ref[, COLNAME_TAXON_REF_DB] <- gsub("^[[:space:]]\\s*(.*?)", "", ref[, COLNAME_TAXON_REF_DB], perl = TRUE)

    return(ref)

}

#' LOAD METADATA TABLE
#' 
#' This function takes a path to a metadata file and reads the table.
#' and reads the table. The first 7 lines are skipped, as they correspond to header lines.
#' Column names are added before the reference database dataframe is returned.
#' 
#' @param mdata_path Path to a metadata table.
#' @return 
#' 
load_metadata <- function(mdata_path) {

    # Read metadata file.
    mdata <- read.csv(
        mdata_path,
        row.names = 1,
        header = TRUE
    )

    return(mdata)
}

#' PREPARE DATAFRAME WITH CONCATENATED MPA-STYLE REPORTS 
#' 
#' This function takes a path to a metadata file and reads the table.
#' and reads the table. The first 7 lines are skipped, as they correspond to header lines.
#' Column names are added before the reference database dataframe is returned.
#' 
#' @param mdata_path Path to a metadata table.
#' @return 
#' 
load_mpa_reports <- function(mpa_reports_dir, verbose = TRUE) {

    # Get paths to MPA-style reports in a specified directory.
    mpa_file_names <- list.files(mpa_reports_dir, pattern = "mpa$")

    # Check if directory really has any MPA-style reports...
    if (length(mpa_file_names) == 0) {
        stop(paste0("No MPA-style reports were found at ", mpa_reports_dir, ". Please review your input."))
    
    # If it does, read MPA-style reports.
    } else {

        # Create empty dataframe to store MPA-style reports.
        mpa_reports <- data.frame(matrix(nrow = 0, ncol = 2))

        # Iterate over files in specified MPA directory (ideally each file should be a different sample)...
        for (file_name in mpa_file_names) {

            # Get sample name.
            sample_name <- gsub("\\..*", "", basename(file_name))

            if (verbose == TRUE) message(paste0("Reading MPA-style report for sample ", sample_name))

            # Read MPA-style report.
            mpa_report <- read.csv(file = paste0(mpa_reports_dir, "/", file_name), header = FALSE, sep = "\t")

            mpa_report_colnames <- c(
                COLNAME_TAXON_MPA, 
                COLNAME_N_FRAG_CLADE_MPA
            )

            # Add column names.
            colnames(mpa_report) <- mpa_report_colnames

            # Add sample name to MPA dataframe.
            mpa_report$sample <- rep(sample_name, times = nrow(mpa_report))

            # Reorder columns.
            mpa_report <- mpa_report[, c("sample", mpa_report_colnames)]

            # Add sample results to dataframe where all results will be stored.
            mpa_reports <- rbind(mpa_reports, mpa_report)
        } 

        return(mpa_reports) 
    }
}

load_std_reports <- function(std_reports_dir, verbose = TRUE) {

    # Get paths to reports in a specified directory.
    std_file_names <- list.files(std_reports_dir, pattern = "kraken$")

    # Check if directory really has any reports...
    if (length(std_file_names) == 0) {
        stop(paste0("No reports were found at ", std_reports_dir, ". Please review your input."))
    
    # If it does, read reports.
    } else {

        # Create empty dataframe to store reports.
        std_reports <- data.frame(matrix(nrow = 0, ncol = 2))

        # Iterate over files in specified reports directory (ideally each file should be a different sample)...
        for (file_name in std_file_names) {

            # Get sample name.
            sample_name <- gsub("\\..*", "", basename(file_name))               
            
            if (verbose == TRUE) message(paste0("Reading standard report for sample ", sample_name))

            # Read report.
            std_report <- read.csv(file = paste0(std_reports_dir, "/", file_name), header = FALSE, sep = "\t")
            
            std_report_colnames <- c(
                COLNAME_PCT_FRAG_CLADE_STD, 
                COLNAME_N_FRAG_CLADE_STD, 
                COLNAME_N_FRAG_TAXON_STD, 
                COLNAME_MINIMISERS_STD, 
                COLNAME_UNIQ_MINIMISERS_STD,
                COLNAME_RANK_STD, 
                COLNAME_NCBI_ID_STD, 
                COLNAME_TAXON_STD
            )

            # Add column names.
            colnames(std_report) <- std_report_colnames

            # Add sample name to report dataframe.
            std_report[, COLNAME_SAMPLE_STD] <- rep(sample_name, times = nrow(std_report))

            # Reorder columns.
            std_report <- std_report[, c(COLNAME_SAMPLE_STD, std_report_colnames)]

            # Remove identation.
            std_report[, COLNAME_TAXON_STD] <- gsub("^[[:space:]]\\s*(.*?)", "", std_report[, COLNAME_TAXON_STD], perl = TRUE)

            # Add sample results to dataframe where all results will be stored.
            std_reports <- rbind(std_reports, std_report)
        } 

        return(std_reports) 
    }
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
addMetadata <- function(dataframe, metadata, columns) {

    # Create temporary sample columns to facilitate the matching between dataframes.
    metadata$sample <- rownames(metadata)

    # Iterate over categories...
    for (column in columns) {

        column <- gsub(" ", "_", column)

        # Add metadata from one dataframe to another based on sample IDs.
        dataframe[, column] <- metadata[, column][match(dataframe$sample, metadata$sample)]
    }

    # Remove temporary sample columns.
    metadata$sample <- NULL

    # Get names of columns that contain results (and not sample names / metadata).
    results_cols <- colnames(dataframe)[!(colnames(dataframe) %in% c("sample", columns))]

    dataframe <- dataframe[, c("sample", columns, results_cols)]

    return(dataframe)
}


############################################
## HELPER FUNCTIONS FOR MPA-STYLE REPORTS ##
#######################################################################################################

is_mpa <- function(report) {
    
    if (COLNAME_TAXON_MPA == COLNAME_TAXON_STD ) {
        stop("The column names are the same in both report formats and therefore it will not
             be able to distinguish between them.")
    } else {
        # If the column "taxon" is present in the report, then it is an MPA-style report.
        if (COLNAME_TAXON_MPA %in% colnames(report)) { 
            return(TRUE)

        # If the column "taxon" is not present...
        } else {

            # If the column "scientific_name" is present, then it is a standard report.
            if (COLNAME_TAXON_STD %in% colnames(report)) {
                return(FALSE)
            
            # If none of those columns are present then it is not an MPA-style or standard report.
            } else {
                stop(paste0(
                    "There is no support for the report format that has been provided. ",
                    "Please review your input. If your report is in standard or MPA format, ",
                    "make sure you load it using load_std_reports() or load_mpa_reports()."
                ))
            }
        }
    }
}

extract_taxon <- function(line, rank, last_in_hierarchy) {

    taxon <- NA

    # Domain.
    if (rank == "D") {

        if (last_in_hierarchy) {
            taxon <- stringr::str_match(line, "^d__(.+)")[2]
        }
        else if (!(last_in_hierarchy)) {
            taxon <- stringr::str_match(line, "^d__\\s*(.*?)\\s*[|](.*?)")[2]
        }

    # Kingdom, phylum, class, order, family, or genus.
    } else if (rank %in% c("K", "P", "C", "O", "F", "G")) { 

        if (last_in_hierarchy) {
            taxon <- stringr::str_match(line, paste0(tolower(rank), "__(.+)"))[2]
        }
        else if (!(last_in_hierarchy)) {
            taxon <- stringr::str_match(line, paste0("(.*?)[|]\\s*", tolower(rank), "__(.*?)\\s*[|](.*?)"))[3]
        }

    # Species.
    } else if (rank == "S") { 

        taxon <- stringr::str_match(line, "s__(.+)")[2] # Will always be the last.

    } else {
        stop(paste0(
            "This function does not support the rank ", rank, ". Supported ranks are 'D' (domain), ", 
            "'K' (kingdom), 'P' (phylum), 'C' (class), 'O' (order), 'F' (family), 'G' (genus) and 'S' (species)."
        ))
    }

    return(taxon)
}

get_association <- function(ranks) {
    names(ranks) <- dplyr::case_when(
        ranks == "D" ~ "domain",
        ranks == "K" ~ "kingdom",
        ranks == "P" ~ "phylum",
        ranks == "C" ~ "class",
        ranks == "O" ~ "order",
        ranks == "F" ~ "family",
        ranks == "G" ~ "genus",
        ranks == "S" ~ "species",
    )
    return(ranks)
}

get_association_2 <- function(ranks) {
    names(ranks) <- dplyr::case_when(
        ranks == "domain" ~ "D",
        ranks == "kingdom" ~ "K",
        ranks == "phylum" ~ "P",
        ranks == "class" ~ "C",
        ranks == "order" ~ "O",
        ranks == "family" ~ "F",
        ranks == "genus" ~ "G",
        ranks == "species" ~ "S",
    )
    return(ranks)
}

sum_domainReads <- function(report, domains) {

    if (is_mpa(report)) {
        
        subset <- report[grep(domains, report[, COLNAME_TAXON_MPA]),]
        sum_domains <- sum(subset[, COLNAME_N_FRAG_CLADE_MPA])

    } else {

        subset <- report[grep(domains, report[, COLNAME_TAXON_STD]),]
        sum_domains <- sum(subset[, COLNAME_N_FRAG_CLADE_STD])

    }

    return(as.numeric(sum_domains))
}

prepare_for_plotDomainReads <- function(report, include_eukaryotes) {

    if(include_eukaryotes) domains <- "^Eukaryota$|^Viruses$|^Archaea$|^Bacteria$"
    else domains <- "^Viruses$|^Archaea$|^Bacteria$"

    if (is_mpa(report)) {
           
        colname_taxon <- COLNAME_TAXON_MPA
        colname_n_frag_clade <- COLNAME_N_FRAG_CLADE_MPA

        report[, COLNAME_TAXON_MPA] <- gsub("d__", "", report[, COLNAME_TAXON_MPA])
        
    } else {

        colname_taxon <- COLNAME_TAXON_STD
        colname_n_frag_clade <- COLNAME_N_FRAG_CLADE_STD

    }

    report <- report[grep(domains, report[, colname_taxon]), ]

    report$colname_taxon <- report[, colname_taxon]
    report$colname_n_frag_clade <- report[, colname_n_frag_clade]

    return(report)
}




#subset_mpa <- function(merged_mpa, include_human = TRUE) {

    # Subsetting family-, genus- and species-level results
#    merged_mpa <- merged_mpa[merged_mpa$last_rank_in_line %in% c("family", "genus", "species"),]

#    if (include_human) {
#        return(merged_mpa)
#    } else {
#        merged_mpa <- merged_mpa[!(merged_mpa$family %in% c("Hominidae", "Homo", "Homo sapiens")), ]
#        return(merged_mpa)
#    }
#}

##################################
## HELPER FUNCTIONS FOR REPORTS ##
#######################################################################################################

get_ClassificationSummary <- function(report) {

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

    colnames(class_unclass_df) <- c("sample", "type", "value")
    class_unclass_df$value <- as.numeric(class_unclass_df$value)

    return(class_unclass_df)
}

get_ProportionClassifiedReads <- function(report) {

    if (is_mpa(report)) {
        stop(paste0(
            "This function is not applicable to MPA-style reports. The purpose of this function ",
            "is to assess the proportion of classified reads relative to the total number of reads, ",
            "but the MPA-style reports do not contain information on unclassified reads."
        ))
    } else {
        colname_n_frag_clade <- COLNAME_N_FRAG_CLADE_STD
        colname_taxon <- COLNAME_TAXON_STD
        colname_sample <- COLNAME_SAMPLE_STD
    }

    proportion_df <- data.frame(matrix(nrow = 0, ncol = 3))

    for (sample in unique(report[, colname_sample])) {
        for (classification in classifications) {

            value <- as.numeric(report[, colname_n_frag_clade][report[, colname_taxon] == classification & report[, colname_sample] == sample])
            class_unclass_df <- rbind(class_unclass_df, c(sample, names(classifications)[classifications == classification], value))
        }
    }

    colnames(class_unclass_df) <- c("sample", "type", "value")
    class_unclass_df$value <- as.numeric(class_unclass_df$value)

    return(class_unclass_df)    
}





subset_reports <- function(merged_reports, include_human = TRUE) {

    # Subsetting family-, genus- and species-level results
    merged_reports <- merged_reports[merged_reports$rank %in% c("F", "G", "S"),]

    if (include_human) {
        return(merged_reports)
    } else {
        merged_reports <- merged_reports[!(merged_reports$scientific_name %in% c("Hominidae", "Homo", "Homo sapiens")),]
        return(merged_reports)
    }
}

report_add_domains <- function(merged_reports, merged_mpa) {

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



#' DETERMINE WIDTH FOR PDF FILE
#' 
#' This function takes a number of elements (e.g. number of samples) and determines the width
#' that should be used for a plot to be saved in PDF format. This function is tailored for the
#' plots generated by plot_mahalanobis_with_metadata() and plot_mahalanobis_with_driver_genes(),
#' so should be used with caution for any other ends.
#' 
#' @param n_elements Number of elements on the x-axis.
#' @param factor Value to help adjust the width.
#' @return Plot width for PDF file.
#' @export
determine_pdf_width <- function(n_elements, factor = 1) {

    # Basic plot width when there is only one element on the x-axis.
    pdf_width <- 1

    # Increment width if there are 2 elements or more.
    for (i in 1:n_elements) {
        pdf_width <- pdf_width + (0.25 * factor)
    }

    return(pdf_width)
}

#' DETERMINE HEIGHT FOR PDF FILE
#' 
#' This function takes a number of elements (e.g. number of genes) and determines the height
#' that should be used for a plot to be saved in PDF format. This function is tailored for the
#' plots generated by plot_mahalanobis_with_metadata() and plot_mahalanobis_with_driver_genes(),
#' so should be used with caution for any other ends.
#' 
#' @param n_elements Number of elements on the y-axis.
#' @param factor Value to help adjust the height.
#' @return Plot height for PDF file.
#' @export
determine_pdf_height <- function(n_elements, factor = 1) {

    # Basic plot height when there is only one element on the y-axis.
    pdf_height <- 2

    # Increment height if there are 2 elements or more.
    for (i in 1:n_elements) {
        pdf_height <- pdf_height + (1 * factor)
    }

    return(pdf_height)
}
