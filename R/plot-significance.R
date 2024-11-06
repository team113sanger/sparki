
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


plotSignificanceSummary <- function(report, return_plot, outdir, prefix = "") {

    summary <- getSignificanceSummary(report)

    plot <- ggplot2::ggplot(
        summary, 
        ggplot2::aes(
            x = get(COLNAME_SIGNIF_SUMMARY_SAMPLE), 
            y = get(COLNAME_SIGNIF_SUMMARY_N_TAXA), 
            fill = get(COLNAME_SIGNIF_SUMMARY_SIGNIF)
        ) 
    ) +
        ggplot2::geom_bar(position = "fill", stat = "identity") +
        ggplot2::facet_wrap(
            ~get(COLNAME_SIGNIF_SUMMARY_RANK),
            ncol = 1
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            # x-axis
            axis.text.x = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_text(size = 16),
            axis.ticks.x = ggplot2::element_blank(),
            # y-axis
            axis.text.y = ggplot2::element_text(size = 15),
            axis.title.y = ggplot2::element_text(size = 16),
            # legend
            legend.text = ggplot2::element_text(size = 15),
            legend.title = ggplot2::element_text(size = 16),
            legend.justification = "top"
        ) +
        ggplot2::xlab("\nSample\n") +
        ggplot2::ylab("Proportion of taxa\n") +
        ggplot2::scale_fill_manual(name = "Significance", values = c("snow3", "indianred1"))

    # Decide what to do with plot based on user-defined options.
    handlePlot(
        plot = plot, prefix = prefix, return_plot = return_plot, 
        filename = "SignificanceSummary.pdf", outdir = outdir, fig_width = 10, 
        fig_height = 5
    )
}
