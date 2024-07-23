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
load_MPAreports <- function(mpa_reports_dir, verbose = TRUE) {

    # Get paths to MPA-style reports in a specified directory.
    mpa_file_names <- list.files(mpa_reports_dir, pattern = "mpa$")

    # Check if directory really has any MPA-style reports...
    if (length(mpa_file_names) == 0) {
        stop(paste0("No MPA-style reports were found at ", mpa_reports_dir, ". Please review your input."))
    } 

    # Create empty dataframe to store MPA-style reports.
    mpa_reports <- data.frame(matrix(nrow = 0, ncol = 2))

    # Iterate over files in specified MPA directory (ideally each file should be a different sample)...
    for (file_name in mpa_file_names) {

        # Get sample name.
        sample_name <- gsub("\\..*", "", basename(file_name))

        if (verbose == TRUE) message(paste0("Reading MPA-style report for sample ", sample_name))

        # Read MPA-style report.
        mpa_report <- read.csv(file = paste0(mpa_reports_dir, "/", file_name), header = FALSE, sep = "\t")

        mpa_report_colnames <- c(COLNAME_MPA_TAXON, COLNAME_MPA_N_FRAG_CLADE)

        # Add column names.
        colnames(mpa_report) <- mpa_report_colnames

        # Add sample name to MPA dataframe.
        mpa_report[, COLNAME_MPA_SAMPLE] <- rep(sample_name, times = nrow(mpa_report))

        # Reorder columns.
        mpa_report <- mpa_report[, c(COLNAME_MPA_SAMPLE, mpa_report_colnames)]

        # Add sample results to dataframe where all results will be stored.
        mpa_reports <- rbind(mpa_reports, mpa_report) 

    }    
    
    return(mpa_reports) 
}

load_STDreports <- function(std_reports_dir, verbose = TRUE) {

    # Get paths to reports in a specified directory.
    std_file_names <- list.files(std_reports_dir, pattern = "kraken$")

    # Check if directory really has any reports...
    if (length(std_file_names) == 0) {
        stop(paste0("No reports were found at ", std_reports_dir, ". Please review your input."))
    }

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
            COLNAME_STD_PCT_FRAG_CLADE, 
            COLNAME_STD_N_FRAG_CLADE, 
            COLNAME_STD_N_FRAG_TAXON, 
            COLNAME_STD_MINIMISERS, 
            COLNAME_STD_UNIQ_MINIMISERS,
            COLNAME_STD_RANK, 
            COLNAME_STD_NCBI_ID, 
            COLNAME_STD_TAXON
        )

        # Add column names.
        colnames(std_report) <- std_report_colnames

        # Add sample name to report dataframe.
        std_report[, COLNAME_STD_SAMPLE] <- rep(sample_name, times = nrow(std_report))

        # Reorder columns.
        std_report <- std_report[, c(COLNAME_STD_SAMPLE, std_report_colnames)]

        # Remove identation.
        std_report[, COLNAME_STD_TAXON] <- gsub(
            "^[[:space:]]\\s*(.*?)", "", 
            std_report[, COLNAME_STD_TAXON], 
            perl = TRUE
        )

        # Add sample results to dataframe where all results will be stored.
        std_reports <- rbind(std_reports, std_report)
    } 

    return(std_reports) 
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
addMetadata <- function(report, metadata, columns) {

    if (is_mpa(report)) {
        colname_sample <- COLNAME_MPA_SAMPLE
    } else {
        colname_sample <- COLNAME_STD_SAMPLE
    }

    # Create temporary sample columns to facilitate the matching between dataframes.
    metadata[["sample"]] <- rownames(metadata)

    # Iterate over categories...
    for (column in columns) {

        column <- gsub(" ", "_", column)

        # Add metadata from one dataframe to another based on sample IDs.
        report[, column] <- metadata[, column][match(report[, colname_sample], metadata[["sample"]])]
    }

    # Remove temporary sample columns.
    metadata$sample <- NULL

    # Get names of columns that contain results (and not sample names / metadata).
    results_cols <- colnames(report)[!(colnames(report) %in% c(colname_sample, columns))]

    report <- report[, c(colname_sample, columns, results_cols)]

    return(report)
}


############################################
## HELPER FUNCTIONS FOR MPA-STYLE REPORTS ##
#######################################################################################################

is_mpa <- function(report) {
    
    if (COLNAME_MPA_TAXON == COLNAME_STD_TAXON ) {
        stop("The column names are the same in both report formats and therefore it will
             not be possible to distinguish between them.")
    } else if (
        !(COLNAME_MPA_TAXON %in% colnames(report)) &&
        !(COLNAME_STD_TAXON %in% colnames(report))
    ) {
        stop(paste0("There is no support for the report format that has been provided. ",
            "Please review your input. If your report is in standard or MPA format, ",
            "make sure you load it using load_STDreports() or load_MPAreports()."))
    }

    # If the column "taxon" is present in the report, then it is an MPA-style report.
    if (COLNAME_MPA_TAXON %in% colnames(report)) return(TRUE)
    else if (COLNAME_STD_TAXON %in% colnames(report)) return(FALSE)
}

extract_taxon <- function(line, rank, last_in_hierarchy) {

    if (!(rank %in% c("D", "K", "P", "C", "O", "F", "G", "S"))) {
        stop(paste0(
            "This function does not support the rank ", rank, ". Supported ranks are 'D' (domain), ", 
            "'K' (kingdom), 'P' (phylum), 'C' (class), 'O' (order), 'F' (family), 'G' (genus) and 'S' (species)."
        ))
    }

    # Domain.
    if (rank == "D") {

        taxon <- ifelse(
            last_in_hierarchy,
            stringr::str_match(line, "^d__(.+)")[2],
            stringr::str_match(line, "^d__\\s*(.*?)\\s*[|](.*?)")[2]
        )

    # Kingdom, phylum, class, order, family, or genus.
    } else if (rank %in% c("K", "P", "C", "O", "F", "G")) { 

        taxon <- ifelse(
            last_in_hierarchy,
            stringr::str_match(line, paste0(tolower(rank), "__(.+)"))[2],
            stringr::str_match(line, paste0("(.*?)[|]\\s*", tolower(rank), "__(.*?)\\s*[|](.*?)"))[3]
        )

    # Species.
    } else if (rank == "S") taxon <- stringr::str_match(line, "s__(.+)")[2] # Will always be the last.

    return(taxon)
}

get_association <- function(ranks) {
    names(ranks) <- dplyr::case_when(
        ranks == "U" ~ NAME_RANK_UNCLASS,
        ranks == "R" ~ NAME_RANK_ROOT,
        ranks == "R1" ~ NAME_RANK_SUBROOT_1,
        ranks == "R2" ~ NAME_RANK_SUBROOT_2,
        ranks == "R3" ~ NAME_RANK_SUBROOT_3,
        ranks == "D" ~ NAME_RANK_DOMAIN,
        ranks == "D1" ~ NAME_RANK_SUBDOMAIN_1,
        ranks == "D2" ~ NAME_RANK_SUBDOMAIN_2,
        ranks == "D3" ~ NAME_RANK_SUBDOMAIN_3,
        ranks == "D4" ~ NAME_RANK_SUBDOMAIN_4,
        ranks == "D5" ~ NAME_RANK_SUBDOMAIN_5,
        ranks == "K" ~ NAME_RANK_KINGDOM,
        ranks == "K1" ~ NAME_RANK_SUBKINGDOM_1,
        ranks == "K2" ~ NAME_RANK_SUBKINGDOM_2,
        ranks == "K3" ~ NAME_RANK_SUBKINGDOM_3,
        ranks == "P" ~ NAME_RANK_PHYLUM,
        ranks == "P1" ~ NAME_RANK_SUBPHYLUM_1,
        ranks == "P2" ~ NAME_RANK_SUBPHYLUM_2,
        ranks == "P3" ~ NAME_RANK_SUBPHYLUM_3,
        ranks == "P4" ~ NAME_RANK_SUBPHYLUM_4,
        ranks == "P5" ~ NAME_RANK_SUBPHYLUM_5,
        ranks == "P6" ~ NAME_RANK_SUBPHYLUM_6,
        ranks == "P7" ~ NAME_RANK_SUBPHYLUM_7,
        ranks == "P8" ~ NAME_RANK_SUBPHYLUM_8,
        ranks == "P9" ~ NAME_RANK_SUBPHYLUM_9,
        ranks == "C" ~ NAME_RANK_CLASS,
        ranks == "C1" ~ NAME_RANK_SUBCLASS_1,
        ranks == "C2" ~ NAME_RANK_SUBCLASS_2,
        ranks == "C3" ~ NAME_RANK_SUBCLASS_3,
        ranks == "C4" ~ NAME_RANK_SUBCLASS_4,
        ranks == "C5" ~ NAME_RANK_SUBCLASS_5,
        ranks == "C6" ~ NAME_RANK_SUBCLASS_6,
        ranks == "O" ~ NAME_RANK_ORDER,
        ranks == "O1" ~ NAME_RANK_SUBORDER_1,
        ranks == "O2" ~ NAME_RANK_SUBORDER_2,
        ranks == "O3" ~ NAME_RANK_SUBORDER_3,
        ranks == "O4" ~ NAME_RANK_SUBORDER_4,
        ranks == "F" ~ NAME_RANK_FAMILY,
        ranks == "F1" ~ NAME_RANK_SUBFAMILY_1,
        ranks == "F2" ~ NAME_RANK_SUBFAMILY_2,
        ranks == "F3" ~ NAME_RANK_SUBFAMILY_3,
        ranks == "F4" ~ NAME_RANK_SUBFAMILY_4,
        ranks == "F5" ~ NAME_RANK_SUBFAMILY_5,
        ranks == "F6" ~ NAME_RANK_SUBFAMILY_6,
        ranks == "F7" ~ NAME_RANK_SUBFAMILY_7,
        ranks == "G" ~ NAME_RANK_GENUS,
        ranks == "G1" ~ NAME_RANK_SUBGENUS_1,
        ranks == "G2" ~ NAME_RANK_SUBGENUS_2,
        ranks == "G3" ~ NAME_RANK_SUBGENUS_3,
        ranks == "S" ~ NAME_RANK_SPECIES,
        ranks == "S1" ~ NAME_RANK_SUBSPECIES_1,
        ranks == "S2" ~ NAME_RANK_SUBSPECIES_2,
        ranks == "S3" ~ NAME_RANK_SUBSPECIES_3,
        ranks == "S4" ~ NAME_RANK_SUBSPECIES_4
    )
    return(ranks)
}

sum_domainReads <- function(report, domains) {

    if (is_mpa(report)) {
        subset <- report[grep(domains, report[, COLNAME_MPA_TAXON]),]
        sum_domains <- sum(subset[, COLNAME_MPA_N_FRAG_CLADE])
    } else {
        subset <- report[grep(domains, report[, COLNAME_STD_TAXON]),]
        sum_domains <- sum(subset[, COLNAME_STD_N_FRAG_CLADE])
    }

    return(as.numeric(sum_domains))
}

is_subrank <- function(rank) {

    if (nchar(rank) > 1) return(TRUE)
    else return(FALSE)

}

###################################
## HELPER FUNCTIONS FOR PLOTTING ##
#######################################################################################################

# Private function.
exportPlot <- function(plot, filename, outdir, fig_width, fig_height) {

    if (missing(fig_width) && missing(fig_height)) {
        fig_width <- 10
        fig_height <- 10
    }

    ggpubr::ggexport(
        plot, 
        filename = paste0(outdir, filename),
        width = fig_width,
        height = fig_height
    )
}

# Private function.
handlePrefix <- function(filename, prefix) {

    if (prefix != "") prefix <- paste0(prefix, "_")
    filename <- paste0(prefix, filename)

    return(filename)
}

# Private function.
handlePlot <- function(plot, prefix, filename, outdir, fig_width, fig_height, return_plot) {

    # If outdir is not NA, it means it has been provided by the user and therefore
    # the plot should be saved to an output directory.
    if (!(is.na(outdir))) {

        # Since the plot will be saved to an output directory, we need to deal with 
        # the file name and ensure a prefix is added in case that was specified by
        # the user.
        filename <- handlePrefix(filename = filename, prefix = prefix)

        # Now we need to export the plot.
        exportPlot(
            plot = plot, filename = filename, outdir = outdir, 
            fig_width = fig_width, fig_height = fig_height
        )

        # The user may want the plot to be returned by the function even if it
        # has been saved to an output directory.
        if (return_plot == TRUE) return(plot)

    # If outdir is NA, it means the user has not provided an output directory and
    # therefore the plot will not be saved. In this case, the plot will be returned
    # regardless of the value of return_plot.
    } else if (is.na(outdir)) { return(plot) }

}

prepare_for_plotDomainReads <- function(report, include_eukaryotes) {

    if(include_eukaryotes) domains <- "^Eukaryota$|^Viruses$|^Archaea$|^Bacteria$"
    else domains <- "^Viruses$|^Archaea$|^Bacteria$"

    if (is_mpa(report)) {
           
        colname_taxon <- COLNAME_MPA_TAXON
        colname_n_frag_clade <- COLNAME_MPA_N_FRAG_CLADE

        report[, COLNAME_MPA_TAXON] <- gsub("d__", "", report[, COLNAME_MPA_TAXON])
        
    } else {

        colname_taxon <- COLNAME_STD_TAXON
        colname_n_frag_clade <- COLNAME_STD_N_FRAG_CLADE

    }

    report <- report[grep(domains, report[, colname_taxon]), ]

    report$colname_taxon <- report[, colname_taxon]
    report$colname_n_frag_clade <- report[, colname_n_frag_clade]

    return(report)
}

# private function
process_barplot_orientation <- function(plot, n_samples, orientation, include_sample_names, factor) {

    # proportion (x-axis) vs samples (y-axis)
    if (orientation %in% c("vertical", "v")) { 

        base_size <- 3
        if (include_sample_names == TRUE) base_size <- 12
        
        pdf_width <- 5
        pdf_height <- determine_pdf_height(
            n_elements = n_samples, 
            base_size = 3,
            factor = factor
        )

    # samples (x-axis) vs proportion (y-axis)
    } else if (orientation %in% c("horizontal", "h")) {

        base_size <- 4
        if (include_sample_names == TRUE) base_size <- 16

        plot <- plot + ggplot2::coord_flip() 

        pdf_height <- 3.5
        pdf_width <- determine_pdf_width(
            n_elements = n_samples, 
            base_size = 4,
            factor = factor
        )

    } else {
        stop(paste0(
            "The orientation provided (", orientation, ") is not valid. ",
            "Please choose from the following options: 'horizontal' or 'h' ",
            "for horizontal plots; 'vertical' or 'v' for vertical plots." 
        ))
    }

    return(list(plot, pdf_width, pdf_height))
}

process_barplot_ticknames <- function(plot, orientation, include_sample_names) {

    # proportion (x-axis) vs samples (y-axis)
    if (orientation %in% c("vertical", "v")) { 

        if (include_sample_names == TRUE) {

            plot <- plot + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 3))

        } else {

            plot <- plot + ggplot2::theme(
                axis.text.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank()
            )

        }

    } else if (orientation %in% c("horizontal", "h")) {

        if (include_sample_names == TRUE) {

            plot <- plot + ggplot2::theme(
                axis.text.x = ggplot2::element_text(size = 3, angle = 90, vjust = 1, hjust = 1)
            )
            
        } else {

            plot <- plot + ggplot2::theme(
                axis.text.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank()
            )

        }

    } else {

        stop(paste0(
            "The orientation provided (", orientation, ") is not valid. ",
            "Please choose from the following options: 'horizontal' or 'h' ",
            "for horizontal plots; 'vertical' or 'v' for vertical plots." 
        ))

    }

    return(plot)
}

# private function
adjust_barplot <- function(plot, n_samples, include_sample_names, orientation, filename) {

    if (include_sample_names) {

        plot <- plot + ggplot2::geom_bar(position = "fill", stat = "identity", width = 0.9) 
        pdf_factor <- 0.6
        if (!(is.na(filename))) filename <- paste0(filename, orientation, "_withSampleNames.pdf")

    } else {

        plot <- plot + ggplot2::geom_bar(position = "fill", stat = "identity", width = 1)
        pdf_factor <- 0.3
        if (!(is.na(filename))) filename <- paste0(filename, orientation, "_withoutSampleNames.pdf")
            
    }

    plot <- process_barplot_ticknames(
        plot = plot, 
        orientation = orientation, 
        include_sample_names = include_sample_names
    )

    processed_plot <- process_barplot_orientation(
        plot = plot,
        n_samples = n_samples,
        orientation = orientation, 
        include_sample_names = include_sample_names,
        factor = pdf_factor
    )

    plot <- processed_plot[[1]]
    pdf_width <- processed_plot[[2]]
    pdf_height <- processed_plot[[3]]

    return(list(plot, filename, pdf_width, pdf_height))
}

prepare_for_plotMinimisers <- function(std_reports, domain) {

    # Subset report to keep only a given domain.
    subset <- std_reports[std_reports[, COLNAME_STD_DOMAIN] == domain, ]

    df <- data.frame(matrix(nrow = 0, ncol = 8))

    # Iterate over all samples....
    for (sample in unique(subset[, COLNAME_STD_SAMPLE])) {

        # Iterate over all taxa in a given domain...
        for (taxon in unique(subset[, COLNAME_STD_TAXON])) {

            # Identify whether a given taxon was identified in a given sample.
            n_res <- length(subset[, COLNAME_STD_RATIO_CLADE][((subset[, COLNAME_STD_TAXON] == taxon) & 
                (subset[, COLNAME_STD_SAMPLE] == sample))])

            # If a given taxon WAS NOT identified in a given sample...
            if(n_res == 0) {
                
                df <- rbind(
                    df, 
                    c(sample, domain, taxon, unique(subset[, COLNAME_STD_RANK][subset[, COLNAME_STD_TAXON] == taxon]),
                      NA, NA, "Non-significant", NA)
                )

            # If a given taxon WAS identified in a given sample...
            } else {

                sample_taxon <- subset[(subset[, COLNAME_STD_TAXON] == taxon) & (subset[, COLNAME_STD_SAMPLE] == sample), ]

                df <- rbind(
                    df, 
                    c(sample, domain, taxon, sample_taxon[, COLNAME_STD_RANK], sample_taxon[, COLNAME_STD_RATIO_CLADE],
                      sample_taxon[, COLNAME_STD_PADJ], sample_taxon[, COLNAME_STD_SIGNIF], sample_taxon[, COLNAME_STD_N_FRAG_CLADE])
                )
            }
        }
    }

    colnames(df) <- c(
        COLNAME_STD_SAMPLE,
        COLNAME_STD_DOMAIN, 
        COLNAME_STD_TAXON, 
        COLNAME_STD_RANK,
        COLNAME_STD_RATIO_CLADE,
        COLNAME_STD_PADJ,
        COLNAME_STD_SIGNIF,
        COLNAME_STD_N_FRAG_CLADE
    )

    df[, COLNAME_STD_RATIO_CLADE] <- as.numeric(df[, COLNAME_STD_RATIO_CLADE])
    df[, COLNAME_STD_PADJ] <- as.numeric(df[, COLNAME_STD_PADJ])
    df[, COLNAME_STD_N_FRAG_CLADE] <- as.numeric(df[, COLNAME_STD_N_FRAG_CLADE])

    return(df)
}

############################################
## HELPER FUNCTIONS FOR transferDomains() ##
#######################################################################################################

# Will be a private function.
retrieve_rankDomains <- function(report_std, report_mpa, verbose = TRUE) {

    # Create vector with ranks present in the standard report provided.
    ranks <- unique(report_std[, COLNAME_STD_RANK])
    ranks <- get_association(ranks)

    # Initialise progress bar.
    if (verbose == TRUE) {
        pb <- txtProgressBar(min = 0, max = nrow(report_std), style = 3)
    }

    # Iterate over each line of the standard report... 
    report_std[, COLNAME_STD_DOMAIN] <- sapply(seq_len(nrow(report_std)), function(x) {

        # Update progress bar.
        if (verbose) setTxtProgressBar(pb, x)

        # Get rank in line.
        rank <- names(ranks)[ranks == report_std[, COLNAME_STD_RANK][x]]

        # If rank is also a column in the MPA-style report...
        # (Note that sub-ranks [e.g. C2, F5 etc] are not present in the MPA format!)
        if (rank %in% colnames(report_mpa)) {

            # Obtain domain.
            domain <- unique(report_mpa[, NAME_RANK_DOMAIN][which(
                report_mpa[, rank] == report_std[, COLNAME_STD_TAXON][x]
            )])
            return(domain)
        
        # As stated above, sub-ranks are not present in MPA-style reports so for
        # these the value returned will be NA.
        } else return(NA)
    })

    if (verbose) cat("\n") # Print line after progress bar is finished.

    return(report_std) # Return the updated version of the standard report!
}

infer_fromAdjacencies <- function(report_std, line) {

    domain_upstream <- report_std[, COLNAME_STD_DOMAIN][line - 1]
    domain_downstream <- report_std[, COLNAME_STD_DOMAIN][line + 1]
                
    if (domain_upstream == domain_downstream) domain <- domain_upstream

    return(domain)
}

determine_subrankDomain <- function(report_std, line, inference) {

    # Initialise domain. 
    domain <- NA

    # Identify sub-rank.
    sub_rank <- report_std[, COLNAME_STD_RANK][line]

    # Create a counter to go up one line at a time to find the nearest rank.
    go_up_one_line <- 1

    # Keep looping until a domain is found...
    while (is.na(domain)) {
        
        # Break loop in case we have reached the beginning of the dataframe
        # and there are no more lines to look at.
        if ((line - go_up_one_line) == 0) {
                   
            warning(paste0(
                "The loop has reached the beginning of the dataframe while ",
                "trying to assign a domain to line ", line, ". Exiting loop..."
            ))
            break
        } 

        # Get rank/sub-rank value above a given sub-rank.
        attempt <- report_std[, COLNAME_STD_RANK][line - go_up_one_line]

        # If the value of "attempt" is a sub-rank, it means we have not reached
        # a rank yet.
        if (is_subrank(attempt)) { 

            # Update counter. 
            go_up_one_line <- go_up_one_line + 1 

        # Alternatively, if the value of "attempt" is not a sub-rank, then we
        # have gotten to a rank!
        } else if (!(is_subrank(attempt))) { 
                    
            # We also need to double check we are looking at the right rank for
            # the sub-rank (e.g. if sub-rank is C2, then the rank should be C,
            # not O, F or anything else).
            if (grepl(attempt, sub_rank)) {

                # Assign new value to domain (i.e. it will no longer be NA).
                domain <- report_std[, COLNAME_STD_DOMAIN][line - go_up_one_line] 

            } else {

                expected_rank <- substring(sub_rank, 1, 1)
                actual_rank <- report_std[, COLNAME_STD_RANK][line - go_up_one_line]

                warning_msg <- paste0(
                    "The nearest rank found for ", sub_rank, " at line ", line, " was ", actual_rank,
                    " while it should have been ", expected_rank, ". "
                )

                if (inference == FALSE) {

                    warning(paste0(warning_msg, "Exiting loop..."))
                    break

                } else {

                    warning(paste0(
                        warning_msg, "Trying to infer the domain based on the sub-rank's adjacencies..."
                    ))

                    #domain <- infer_fromAdjacencies()
                }
            }
        } 
    }

    return(domain)
}

# Will be a private function.
retrieve_subrankDomains <- function(report_std, report_mpa, inference = TRUE, verbose = TRUE) {

    # Initialise progress bar.
    if (verbose) {
        cat("\nChecking sub-ranks\n")
        pb <- txtProgressBar(min = 0, max = nrow(report_std), style = 3)
    }

    # Iterate over lines in standard report...
    for (i in seq_len(nrow(report_std))) {

        # Update progress bar.
        if (verbose) setTxtProgressBar(pb, i) 

        if (!is.na(report_std[, COLNAME_STD_DOMAIN][i]) || report_std[, COLNAME_STD_RANK][i] %in% c("U", "R", "R1", "R2", "R3")) {
            next
        }
        # Every time domain = NA (except when rank is related to root or unclassified), this means that
        # the line corresponds to a sub-rank. The chunck of code below will look for the nearest rank upstream 
        # (i.e. if sub-rank is F2, the nearest rank upstream will be F) and assign the nearest rank's domain to
        # the sub-rank.

        # Determine sub-rank's domain.
        domain <- determine_subrankDomain(report_std, i, inference)

        # Assign domain value to sub-rank.
        report_std[, COLNAME_STD_DOMAIN][i] <- domain

    }
    
    cat("\n")

    return(report_std)

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

get_ProportionClassifiedReads <- function(report) {

    if (is_mpa(report)) {
        stop(paste0(
            "This function is not applicable to MPA-style reports. The purpose of this function ",
            "is to assess the proportion of classified reads relative to the total number of reads, ",
            "but the MPA-style reports do not contain information on unclassified reads."
        ))
    } else {
        colname_n_frag_clade <- COLNAME_STD_N_FRAG_CLADE
        colname_taxon <- COLNAME_STD_TAXON
        colname_sample <- COLNAME_STD_SAMPLE
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


subset_STDreport <- function(report, include_human = TRUE) {

    # Subsetting family-, genus- and species-level results.
    report <- report[report[, COLNAME_STD_RANK] %in% c("F", "G", "S"), ]

    # Filter out human results.    
    if (!(include_human)) {
        report <- report[!(report[, COLNAME_STD_TAXON] %in% c("Hominidae", "Homo", "Homo sapiens")), ]
    }

    return(report)
}

get_nTaxaInRank <- function(report, rank, sample) {

    n_taxa_in_rank <- length(unique(
        report[, COLNAME_STD_TAXON][report[, COLNAME_STD_RANK] == rank & report[, COLNAME_STD_SAMPLE] == sample]
    ))

    return(n_taxa_in_rank)
}

#' DETERMINE WIDTH FOR PDF FILE
#' 
#' This function takes a number of elements (e.g. number of samples) and determines the width
#' that should be used for a plot to be saved in PDF format. This function is tailored for the
#' plots generated by plot_mahalanobis_with_metadata() and plot_mahalanobis_with_driver_genes(),
#' so should be used with caution for any other ends.
#' 
#' @param n_elements Number of elements on the x-axis.
#' @param base_size Base plot width when there is only one element on the x-axis.
#' @param factor Value to help adjust the width.
#' @return Plot width for PDF file.
#' @export
determine_pdf_width <- function(n_elements, base_size = 1, factor = 1) {

    # Increment width if there are 2 elements or more.
    for (i in seq_len(n_elements)) {
        pdf_width <- base_size + (0.25 * factor)
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
#' @param base_size Base plot height when there is only one element on the y-axis.
#' @param factor Value to help adjust the height.
#' @return Plot height for PDF file.
#' @export
determine_pdf_height <- function(n_elements, base_size = 2, factor = 1) {

    # Increment height if there are 2 elements or more.
    for (i in seq_len(n_elements)) {
        pdf_height <- base_size + (1 * factor)
    }

    return(pdf_height)
}

#' PARSE SYMBOL-DELIMITED LIST TO GET INDIVIDUAL ELEMENTS
#' 
#' This function takes a symbol-delimited list of elements and splits it up into individual elements.
#' Symbols can be anything, e.g. ",", ";", "/", "//", "@", "." etc.
#' 
#' @param del_list A list with symbol-delimited values; the symbol can be a comma, for example.
#' @return A vector with individual elements.
#' @export
parse_delimited_list <- function(del_list, delimiter) {

    # Split comma-separated list.
    elements <- unlist(stringr::str_split(del_list, delimiter))

    return(elements)
}
