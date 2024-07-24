#' PLOT NUMBERS OF CLASSIFIED/UNCLASSIFIED READS PER SAMPLE (VIOLIN PLOT)
#' 
#' This function takes a standard Kraken2 report and creates a violin plot showing
#' the number of classified and unclassified reads per sample.
#' 
#' @param report Standard Kraken2 report.
#' @param return_plot Whether plot should be returned.
#' @param outdir Output directory where the plot should be saved.
#' @return 
#' 
plotClassificationSummary_violin <- function(report, return_plot, outdir, prefix = "") {

    # Stop if report provided is in MPA-style format.
    if (is_mpa(report)) stop(paste0("This function is not applicable to MPA-style reports."))

    # Assign NA to outdir in case it has not been provided by the user.
    if (missing(outdir)) outdir <- NA

    # Prepare data for plotting.
    class_unclass_df <- getClassificationSummary(report)

    # Create violin plot.
    plot <- ggplot2::ggplot(
            class_unclass_df, 
            ggplot2::aes(x = type, y = log10(n_reads))
        ) +
        ggplot2::geom_violin(scale = "width", fill = "white", color = "black") +
        ggplot2::geom_line(ggplot2::aes(group = sample), alpha = 0.25) +
        ggplot2::geom_point(ggplot2::aes(color = type), alpha = 0.5) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            # x-axis
            axis.text.x = ggplot2::element_text(size = 12, vjust = 0.5),
            axis.title.x = ggplot2::element_text(size = 14),
            # y-axis
            axis.text.y = ggplot2::element_text(size = 12),
            axis.title.y = ggplot2::element_text(size = 14, angle = 90),
            # legend
            legend.position = "none"
        ) +
        ggplot2::xlab("\nRead classification") +
        ggplot2::ylab(expression("log"[10]~"(# reads)")) +
        ggplot2::scale_color_manual(
            values = c("indianred2", "royalblue")
        )

    # Decide what to do with plot based on user-defined options.
    handlePlot(
        plot = plot, prefix = prefix, return_plot = return_plot, 
        filename = "nReads_classified_vs_unclassified_absNumbers_violinPlot.pdf", 
        outdir = outdir, fig_width = 3, fig_height = 4
    )
}

#' PLOT NUMBERS OF CLASSIFIED/UNCLASSIFIED READS PER SAMPLE (BAR PLOT)
#' 
#' This function takes a standard Kraken2 report and creates a bar plot showing
#' the number of classified and unclassified reads per sample.
#' 
#' @param report Standard Kraken2 report.
#' @param include_sample_names Whether to include sample names in the plot.
#' @param orientation Whether plot should be horizontally or vertically oriented.
#' @param return_plot Whether plot should be returned.
#' @param outdir Output directory where the plot should be saved.
#' @param prefix Prefix to be added to output plot name.
#' @return 
#' 
plotClassificationSummary_barplot <- function(
    report, include_sample_names = FALSE, orientation = "vertical", 
    return_plot, outdir, prefix = ""
) {

    # Stop if report provided is in MPA-style format.
    if (is_mpa(report)) stop(paste0("This function is not applicable to MPA-style reports."))

    # Assign NA to outdir in case it has not been provided by the user.
    if (missing(outdir)) outdir <- NA

    # Prepare data for plotting.
    class_unclass_df <- getClassificationSummary(report)

    # Create bar plot.
    plot <- ggplot2::ggplot(
            class_unclass_df, 
            ggplot2::aes(x = n_reads, y = sample, fill = type)
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
            # x-axis
            axis.text.x = ggplot2::element_text(size = 12, vjust = 0.5),
            axis.title.x = ggplot2::element_text(size = 14),
            # y-axis
            axis.text.y = ggplot2::element_text(size = 12),
            axis.title.y = ggplot2::element_text(size = 14, angle = 90),
            # legend
            legend.text = ggplot2::element_text(size = 12),
            legend.title = ggplot2::element_text(size = 14),
            legend.position = "right",
            # title
            plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = "bold"),
            # margin
            plot.margin = ggplot2::unit(c(0, 0, 0, 0), "pt")
        ) +
        ggplot2::scale_fill_brewer(
            palette = "Paired",
            name = "Read\nclassification",
            labels = c("Classified", "Unclassified")
        ) +
        ggplot2::xlab("Proportion of reads\n") +
        ggplot2::ylab("Sample\n")

    adjusted_plot <- adjust_barplot(
        plot = plot, n_samples = length(unique(report[, COLNAME_STD_SAMPLE])),
        include_sample_names = include_sample_names, orientation = orientation, 
        filename = "nReads_classified_vs_unclassified_proportion_perSample_barPlot"
    ) 

    # Decide what to do with plot based on user-defined options.
    handlePlot(
        plot = adjusted_plot[[1]], prefix = prefix, return_plot = return_plot, 
        filename = adjusted_plot[[2]], outdir = outdir, fig_width = adjusted_plot[[3]], 
        fig_height = adjusted_plot[[4]]
    )
}


#' PLOT NUMBERS OF CLASSIFIED READS PER DOMAIN (VIOLIN PLOT)
#' 
#' This function takes a Kraken2 report, either in standard or MPA-style format, and creates
#' a violin plot showing the number of classified reads per domain.
#' 
#' @param report Kraken2 report, either in standard or MPA-style format.
#' @param include_eukaryotes Whether to include eukaryotes in the plot.
#' @param return_plot Whether plot should be returned.
#' @param outdir Output directory where the plot should be saved.
#' @param prefix Prefix to be added to output plot name.
#' @return 
#' 
plotDomainReads_violin <- function(
    report, include_eukaryotes = FALSE, return_plot, outdir, 
    prefix = "", fig_width, fig_height
) {

    if (include_eukaryotes == TRUE) {
        x_lab <- "\nProportion of classified reads\n(all domains)"
        filename <- "nReadsDomains_violin_with_eukaryotes.pdf"
        plot_width <- 5
        plot_height <- 5
    } else if (include_eukaryotes == FALSE) {
        x_lab <- "\nProportion of classified reads\n(non-eukaryotes only)"
        filename <- "nReadsDomains_violin_without_eukaryotes.pdf"
        plot_width <- 5
        plot_height <- 3.5
    } else {
        stop(paste0(
            "The value of include_eukaryotes should be either TRUE or FALSE, but ",
            include_eukaryotes, " was provided. Please review your input."
        ))
    }

    report <- prepare_for_plotDomainReads(report, include_eukaryotes)

    # Create violin plot.
    plot <- ggplot2::ggplot(
            report, 
            ggplot2::aes(x = colname_taxon, y = log10(colname_n_frag_clade), fill = colname_taxon)
        ) +
        ggplot2::geom_violin(scale = "width") +
        ggplot2::geom_jitter(size = 0.5) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            # x-axis
            axis.text.x = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            strip.text.x = ggplot2::element_text(size = 15),
            # y-axis
            axis.text.y = ggplot2::element_text(size = 15),
            axis.title.y = ggplot2::element_text(size = 16),
            # legend
            legend.position = "none"
        ) +
        ggplot2::ylab(expression("log"[10]~"(# classified reads)")) +
        ggplot2::scale_fill_manual(values = c("gold1", "royalblue", "snow2", "indianred2")) +
        ggplot2::facet_wrap(~colname_taxon, scales = "free")

    # Decide what to do with plot based on user-defined options.
    handlePlot(
        plot = plot, prefix = prefix, return_plot = return_plot, filename = filename,
        outdir = outdir, fig_width = plot_width, fig_height = plot_height
    )
}

#' PLOT NUMBERS OF CLASSIFIED READS PER DOMAIN (BAR PLOT)
#' 
#' This function takes a Kraken2 report, either in standard or MPA-style format, and creates
#' a bar plot showing the number of classified reads per domain.
#' 
#' @param report Kraken2 report, either in standard or MPA-style format.
#' @param include_sample_names Whether sample names should be displayed.
#' @param orientation Whether plot should be horizontally or vertically oriented.
#' @param include_eukaryotes Whether to include eukaryotes in the plot.
#' @param return_plot Whether plot should be returned.
#' @param outdir Output directory where the plot should be saved.
#' @param prefix Prefix to be added to output plot name.
#' @return 
#' 
plotDomainReads_barplot <- function(
    report, include_sample_names = FALSE, orientation = "vertical", 
    include_eukaryotes = FALSE, return_plot, outdir, prefix = ""
) {

    # Assign NA to outdir in case it has not been provided by the user.
    if (missing(outdir)) outdir <- NA

    report <- prepare_for_plotDomainReads(report, include_eukaryotes)

    if (include_eukaryotes) {
        x_lab <- "\nProportion of classified reads\n(all domains)"
        filename <- "nReadsDomains_barplot_with_eukaryotes"
        colours <- c("gold1", "royalblue", "snow4", "indianred2")
    } else {
        x_lab <- "\nProportion of classified reads\n(non-eukaryotes only)"
        filename <- "nReadsDomains_barplot_without_eukaryotes"
        colours <- c("gold1", "royalblue", "indianred2")
    }

    plot <- ggplot2::ggplot(
            report, ggplot2::aes(x = colname_n_frag_clade, y = sample, fill = colname_taxon)
        ) +
        ggplot2::geom_bar(position = "fill", stat = "identity") +
        ggplot2::theme_void() +
        ggplot2::theme(
            # x-axis
            axis.text.x = ggplot2::element_text(size = 15),
            axis.title.x = ggplot2::element_text(size = 16),
            # y-axis
            axis.text.y = ggplot2::element_text(size = 15),
            axis.title.y = ggplot2::element_text(size = 16, angle = 90),
            # legend
            legend.text = ggplot2::element_text(size = 15),
            legend.title = ggplot2::element_text(size = 16),
            legend.justification = "top"
        ) +
        ggplot2::xlab(x_lab) +
        ggplot2::ylab("Sample\n") +
        ggplot2::scale_fill_manual(
            name = "Domain", values = colours
        ) +
        ggplot2::geom_vline(xintercept = c(0.25, 0.5, 0.75), linetype = "dashed")

    adjusted_plot <- adjust_barplot(
        plot = plot, n_samples = length(unique(report[, COLNAME_STD_SAMPLE])),
        include_sample_names = include_sample_names, orientation = orientation, 
        filename = filename
    ) 

    # Decide what to do with plot based on user-defined options.
    handlePlot(
        plot = adjusted_plot[[1]], prefix = prefix, return_plot = return_plot, 
        filename = adjusted_plot[[2]], outdir = outdir, fig_width = adjusted_plot[[3]], 
        fig_height = adjusted_plot[[4]]
    )
}

#' PLOT PROPORTION OF MINIMISERS AND SIGNIFICANCE PER SAMPLE
#' 
#' This function takes a Kraken2 report, either in standard or MPA-style format, and creates
#' a bar plot showing the number of classified reads per domain.
#' 
#' @param report Kraken2 report, either in standard or MPA-style format.
#' @param include_sample_names Whether sample names should be displayed.
#' @param orientation Whether plot should be horizontally or vertically oriented.
#' @param include_eukaryotes Whether to include eukaryotes in the plot.
#' @param return_plot Whether plot should be returned.
#' @param outdir Output directory where the plot should be saved.
#' @param prefix Prefix to be added to output plot name.
#' @return 
#' 
plotMinimisers_dotplot <- function(
    std_reports, domain, fig_width, fig_height, return_plot, outdir, prefix = ""
) {

    # Assign NA to outdir in case it has not been provided by the user.
    if (missing(outdir)) outdir <- NA

    df <- prepare_for_plotMinimisers(
        std_reports = std_reports,
        domain = domain
    )
    df[["sample"]] <- df[, COLNAME_STD_SAMPLE]
    df[["taxon"]] <- df[, COLNAME_STD_TAXON]
    df[["rank"]] <- df[, COLNAME_STD_RANK]
    df[["ratio_clade"]] <- df[, COLNAME_STD_RATIO_CLADE]
    df[["padj"]] <- df[, COLNAME_STD_PADJ]
    df[["significance"]] <- df[, COLNAME_STD_SIGNIF]
    df[["n_frag_clade"]] <- df[, COLNAME_STD_N_FRAG_CLADE]

    plot <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = sample, y = taxon, fill = ratio_clade, color = significance, size = log10(n_frag_clade))
    ) +
        ggplot2::geom_point(shape = 21, stroke = 1.25) +
        ggplot2::facet_grid(
            rows = ggplot2::vars(rank), 
            scales = "free_y",
            space = "free_y"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
            axis.text.y = ggplot2::element_text(size = 12),
            legend.text = ggplot2::element_text(size = 15),
            legend.text.align = 0,
            legend.title = ggplot2::element_text(size = 15),
            axis.title = ggplot2::element_text(size = 15),
            plot.title = ggplot2::element_text(size = 16.5, face = "bold", hjust = 0.5),
            strip.text.y = ggplot2::element_text(size = 15, face = "bold")
        ) +
        ggplot2::scale_fill_gradientn(
            colors = colorRampPalette((RColorBrewer::brewer.pal(9, "Reds")))(100),
            name = "Ratio between unique minimisers found in sample\nand total clade-level minimisers in database",
            limits = c(0,1)
        ) +
        ggplot2::scale_color_manual(
            values = c("snow3", "black"),
            name = "Significance",
            labels = c("Non-significant", expression("Adjusted p-value"<="0.05"))
        ) +
        ggplot2::scale_size_continuous(name = expression("log"[10]~"(clade-level reads)")) +
        ggplot2::xlab("Sample") +
        ggplot2::ylab("Taxon") 

    handlePlot(
        plot = plot, prefix = prefix, 
        return_plot = return_plot, 
        filename = "MINIMISERS.pdf", 
        outdir = outdir, 
        fig_width = fig_width, 
        fig_height = fig_height
    )

}