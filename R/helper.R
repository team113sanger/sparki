###############################################
## HELPER FUNCTIONS FOR STATISTICAL ANALYSES ##
#######################################################################################################

calculate_p_value <- function(sample_n_uniq_minimisers_taxon, db_n_minimisers_taxon, sample_size) {

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

################################################################################
## HELPER FUNCTIONS FOR DISTINGUISHING BETWEEN MPA-STYLE AND STANDARD REPORTS ##
#######################################################################################################

is_mpa <- function(report) {

    ifelse(
        COLNAME_MPA_TAXON_LEAF %in% colnames(report),
        return(TRUE),
        return(FALSE)
    )
}

######################################################
## HELPER FUNCTIONS FOR PREPARING DATA FOR PLOTTING ##
#######################################################################################################

prepare_for_plotDistribution <- function(report) {

    if (is_mpa(report)) stop(paste0("This function does not support MPA-style reports."))

    report <- report |>
        # Filter out rows.
        dplyr::filter(!(!!as.name(COLNAME_STD_TAXON) %in% c("root", "unclassified"))) |>
        # Rename elements in column.
        dplyr::mutate(!!COLNAME_STD_RANK := dplyr::case_when(
            !!as.name(COLNAME_STD_RANK) == "D" ~ "Domain",
            !!as.name(COLNAME_STD_RANK) == "K" ~ "Kingdom",
            !!as.name(COLNAME_STD_RANK) == "P" ~ "Phylum",
            !!as.name(COLNAME_STD_RANK) == "C" ~ "Class",
            !!as.name(COLNAME_STD_RANK) == "O" ~ "Order",
            !!as.name(COLNAME_STD_RANK) == "F" ~ "Family",
            !!as.name(COLNAME_STD_RANK) == "G" ~ "Genus",
            !!as.name(COLNAME_STD_RANK) == "S" ~ "Species"
        ))

    return(report)
}

prepare_for_plotDomainReads <- function(report, include_eukaryotes) {

    if (is_mpa(report)) stop(paste0("This function does not support MPA-style reports."))

    domains <- "Viruses|Archaea|Bacteria"
    if (include_eukaryotes == TRUE) domains <- paste0("Eukaryota|", domains)

    report <- report |> 
        # Select relevant columns.
        dplyr::select(
            !!COLNAME_STD_SAMPLE,
            !!COLNAME_STD_TAXON,
            !!COLNAME_STD_N_FRAG_CLADE) |>
        # Select relevant rows.
        dplyr::filter(grepl(!!as.name(COLNAME_STD_TAXON), pattern = domains)) |>
        # Rename columns.
        dplyr::rename(
            !!COLNAME_DOMAIN_READS_SAMPLE := !!COLNAME_STD_SAMPLE,
            !!COLNAME_DOMAIN_READS_TAXON := !!COLNAME_STD_TAXON,
            !!COLNAME_DOMAIN_READS_N_FRAG := !!COLNAME_STD_N_FRAG_CLADE
        ) |>
        # Create column with log-transformed number of clade-level fragments.
        dplyr::mutate(
            !!COLNAME_DOMAIN_READS_LOG_N_FRAG := log10(!!as.name(COLNAME_DOMAIN_READS_N_FRAG))
        )

    return(report)
}

prepare_for_plotMinimisers <- function(report, domain_of_interest) {

    if (is_mpa(report)) stop(paste0("This function does not support MPA-style reports."))

    subset <- report |>
        # Select relevant rows.
        dplyr::filter(domain == domain_of_interest) |>
        # Select relevant columns.
        dplyr::select(
            !!COLNAME_STD_SAMPLE,
            !!COLNAME_STD_TAXON,
            !!COLNAME_STD_RANK,
            domain,
            !!COLNAME_STD_N_FRAG_CLADE,
            !!COLNAME_STD_RATIO_CLADE,
            !!COLNAME_STD_PADJ,
            !!COLNAME_STD_SIGNIF
        ) |>
        # Create column with log-transformed number of clade-level fragments.
        dplyr::mutate(
            !!COLNAME_STD_LOG_N_FRAG_CLADE := log10(!!as.name(COLNAME_STD_N_FRAG_CLADE))
        ) |>
        # Reorder columns.
        dplyr::relocate(
            !!COLNAME_STD_LOG_N_FRAG_CLADE, .after = !!COLNAME_STD_N_FRAG_CLADE
        ) |>
        # Replace values in column.
        dplyr::mutate(!!COLNAME_STD_RANK := dplyr::case_when(
            !!as.name(COLNAME_STD_RANK) == "F" ~ "Family",
            !!as.name(COLNAME_STD_RANK)  == "G" ~ "Genus",
            !!as.name(COLNAME_STD_RANK) == "S" ~ "Species"
        ))

    return(subset)
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
