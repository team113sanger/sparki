

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
    report, domain, fig_width, fig_height, return_plot, outdir, prefix = ""
) {

    # Assign NA to outdir in case it has not been provided by the user.
    if (missing(outdir)) outdir <- NA

    df <- prepare_for_plotMinimisers(report, domain)

    if ("Significant" %in% df[[COLNAME_STD_SIGNIF]] && "Non-significant" %in% df[[COLNAME_STD_SIGNIF]]) {
        ggplot_colours <- c("snow3", "black")
        ggplot_labels <- c("Non-significant", expression("Adjusted p-value"<="0.05"))
    } else if ("Significant" %in% df[[COLNAME_STD_SIGNIF]] && !("Non-significant" %in% df[[COLNAME_STD_SIGNIF]])) {
        ggplot_colours <- "black"
        ggplot_labels <- expression("Adjusted p-value"<="0.05")
    } else if (!("Significant" %in% df[[COLNAME_STD_SIGNIF]]) && "Non-significant" %in% df[[COLNAME_STD_SIGNIF]]) {
        ggplot_colours <- "snow3"
        ggplot_labels <- "Non-significant"
    }

    n_all_samples <- report[[COLNAME_STD_SAMPLE]] |> unique() |> length()
    n_subset_samples <- df[[COLNAME_STD_SAMPLE]] |> unique() |> length()

    plot <- ggplot2::ggplot(
        df,
        ggplot2::aes(
            x = get(COLNAME_STD_SAMPLE),
            y = get(COLNAME_STD_TAXON),
            fill = get(COLNAME_STD_RATIO_CLADE),
            color = get(COLNAME_STD_SIGNIF),
            size = get(COLNAME_STD_LOG_N_FRAG_CLADE)
        )
    ) +
        ggplot2::geom_point(shape = 21, stroke = 1.25) +
        ggplot2::facet_grid(
            rows = ggplot2::vars(get(COLNAME_STD_RANK)),
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
            values = ggplot_colours,
            name = "Significance",
            labels = ggplot_labels
        ) +
        ggplot2::scale_size_continuous(name = expression("log"[10]~"(clade-level reads)")) +
        ggplot2::xlab("Sample") +
        ggplot2::ylab("Taxon") +
        ggplot2::ggtitle(
            paste0("Showing ", n_subset_samples, " (out of ",
            n_all_samples,") with results for domain: ", domain)
        )

    handlePlot(
        plot = plot, prefix = prefix,
        return_plot = return_plot,
        filename = paste0(domain, "_MINIMISERS.pdf"),
        outdir = outdir,
        fig_width = fig_width,
        fig_height = fig_height
    )

}

plotTaxon_minimisers <- function(report, taxon, domain) {

    df <- prepare_for_plotMinimisers(report, domain)
    df <- df |> dplyr::filter(!!as.name(COLNAME_STD_TAXON) == taxon)

    ggplot2::ggplot(df, ggplot2::aes(x = "", y = get(COLNAME_STD_RATIO_CLADE))) +
        ggplot2::geom_violin(scale = "width", color = "black", fill = "palegreen1") +
        ggplot2::geom_jitter(ggplot2::aes(size = get(COLNAME_STD_LOG_N_FRAG_CLADE)), color = "black") +
        ggplot2::theme_classic() +
        ggplot2::theme(
            # x-axis
            axis.text.x = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_blank(),
            # y-axis
            axis.text.y = ggplot2::element_text(size = 12),
            axis.title.y = ggplot2::element_text(size = 14, angle = 90),
            # legend
            legend.position = "none"
        ) +
        ggplot2::ylab("Proportion of\nclade-level minimisers")
}
