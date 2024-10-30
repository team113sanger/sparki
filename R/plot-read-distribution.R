
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

    report[[COLNAME_STD_RANK]] = factor(
        report[[COLNAME_STD_RANK]], 
        levels = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    )

    return(report)
}

plotDistribution_histogram <- function(report, return_plot, outdir, prefix = "") {

    if (is_mpa(report)) stop(paste0("This function does not support MPA-style reports."))

    report <- prepare_for_plotDistribution(report)

    plot <- ggplot2::ggplot(
        report, ggplot2::aes(x = get(COLNAME_STD_N_FRAG_CLADE), fill = get(COLNAME_STD_RANK))
    ) +
            ggplot2::geom_histogram() +
            ggplot2::facet_wrap(~get(COLNAME_STD_RANK), nrow = 2, scale = "free") +
            ggplot2::theme_classic() +
            ggplot2::theme(
                # x-axis
                axis.text.x = ggplot2::element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
                axis.title.x = ggplot2::element_text(size = 14),
                # y-axis
                axis.text.y = ggplot2::element_text(size = 12),
                axis.title.y = ggplot2::element_text(size = 14),
                # legend
                legend.position = "none"
            ) +
            ggplot2::xlab("\n# Clade-level fragments") +
            ggplot2::ylab("Frequency\n")

    # Decide what to do with plot based on user-defined options.
    handlePlot(
        plot = plot, prefix = prefix, return_plot = return_plot, filename = "distribution_histogram.pdf",
        outdir = outdir, fig_width = 8, fig_height = 5
    )    
}

plotDistribution_violin <- function(report, return_plot, outdir, prefix = "") {

    if (is_mpa(report)) stop(paste0("This function does not support MPA-style reports."))

    report <- prepare_for_plotDistribution(report)

    plot <- ggplot2::ggplot(
        report, 
        ggplot2::aes(
            x = get(COLNAME_STD_RANK),
            y = get(COLNAME_STD_N_FRAG_CLADE),
            fill = get(COLNAME_STD_RANK)
        )
    ) +
            ggplot2::geom_violin(scale = "width", color = "black") +
            ggplot2::theme_classic() +
            ggplot2::theme(
                # x-axis
                axis.text.x = ggplot2::element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
                axis.title.x = ggplot2::element_text(size = 14),
                # y-axis
                axis.text.y = ggplot2::element_text(size = 12),
                axis.title.y = ggplot2::element_text(size = 14),
                # legend
                legend.position = "none"
            ) +
            ggplot2::xlab("Rank") +
            ggplot2::ylab("# Clade-level\nfragments\n")

    # Decide what to do with plot based on user-defined options.
    handlePlot(
        plot = plot, prefix = prefix, return_plot = return_plot, filename = "distribution_violin.pdf",
        outdir = outdir, fig_width = 8, fig_height = 5
    )    
}
