prepare_for_plotDomainReads <- function(report, include_eukaryotes) {
  if (is_mpa(report)) stop(paste0("This function does not support MPA-style reports."))

  domains <- "Viruses|Archaea|Bacteria"
  if (include_eukaryotes == TRUE) domains <- paste0("Eukaryota|", domains)

  report <- report |>
    # Select relevant columns.
    dplyr::select(
      !!COLNAME_STD_SAMPLE,
      !!COLNAME_STD_TAXON,
      !!COLNAME_STD_N_FRAG_CLADE
    ) |>
    # Select relevant rows.
    dplyr::filter(grepl(!!as.name(COLNAME_STD_TAXON), pattern = domains)) |>
    # Rename columns.
    dplyr::rename(
      !!COLNAME_DOMAIN_READS_SAMPLE := !!COLNAME_STD_SAMPLE,
      !!COLNAME_DOMAIN_READS_TAXON := !!COLNAME_STD_TAXON,
      !!COLNAME_DOMAIN_READS_N_FRAG := !!COLNAME_STD_N_FRAG_CLADE
    ) |>
    #  Create column with log-transformed number of clade-level fragments.
    dplyr::mutate(
      !!COLNAME_DOMAIN_READS_LOG_N_FRAG := log10(!!as.name(COLNAME_DOMAIN_READS_N_FRAG))
    )

  return(report)
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
    prefix = "", fig_width, fig_height) {
  if (!(include_eukaryotes %in% c(TRUE, FALSE))) {
    stop(paste0(
      "The value of include_eukaryotes should be either TRUE or FALSE, but ",
      include_eukaryotes, " was provided. Please review your input."
    ))
  }

  x_lab <- ifelse(
    include_eukaryotes == TRUE,
    "\nProportion of classified reads\n(all domains)",
    "\nProportion of classified reads\n(non-eukaryotes only)"
  )
  filename <- ifelse(
    include_eukaryotes == TRUE, "nReadsDomains_violin_with_eukaryotes.pdf",
    "nReadsDomains_violin_without_eukaryotes.pdf"
  )
  plot_height <- ifelse(include_eukaryotes == TRUE, 5, 3.5)

  # Prepare data for plotting.
  report <- prepare_for_plotDomainReads(report, include_eukaryotes)

  #  Create violin plot.
  plot <- ggplot2::ggplot(
    report,
    ggplot2::aes(
      x = get(COLNAME_DOMAIN_READS_TAXON),
      y = get(COLNAME_DOMAIN_READS_LOG_N_FRAG),
      fill = get(COLNAME_DOMAIN_READS_TAXON)
    )
  ) +
    ggplot2::geom_violin(scale = "width") +
    ggplot2::geom_jitter(size = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      #  x-axis
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = 15),
      # y-axis
      axis.text.y = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_text(size = 12),
      # legend
      legend.position = "none"
    ) +
    ggplot2::ylab(expression("log"[10] ~ "(# classified reads)")) +
    ggplot2::scale_fill_manual(values = c("gold1", "royalblue", "snow2", "indianred2")) +
    ggplot2::facet_wrap(~ get(COLNAME_DOMAIN_READS_TAXON), scales = "free")

  #  Decide what to do with plot based on user-defined options.
  handlePlot(
    plot = plot, prefix = prefix, return_plot = return_plot, filename = filename,
    outdir = outdir, fig_width = 5, fig_height = plot_height
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
    include_eukaryotes = FALSE, return_plot, outdir, prefix = "") {
  #  Assign NA to outdir in case it has not been provided by the user.
  if (missing(outdir)) outdir <- NA

  if (include_eukaryotes) {
    x_lab <- "\nProportion of classified reads\n(all domains)"
    filename <- "nReadsDomains_barplot_with_eukaryotes"
    colours <- c("gold1", "royalblue", "snow4", "indianred2")
  } else {
    x_lab <- "\nProportion of classified reads\n(non-eukaryotes only)"
    filename <- "nReadsDomains_barplot_without_eukaryotes"
    colours <- c("gold1", "royalblue", "indianred2")
  }

  report <- prepare_for_plotDomainReads(report, include_eukaryotes)

  plot <- ggplot2::ggplot(
    report,
    ggplot2::aes(
      x = get(COLNAME_DOMAIN_READS_N_FRAG),
      y = get(COLNAME_DOMAIN_READS_SAMPLE),
      fill = get(COLNAME_DOMAIN_READS_TAXON)
    )
  ) +
    ggplot2::geom_bar(position = "fill", stat = "identity") +
    ggplot2::theme_void() +
    ggplot2::theme(
      #  x-axis
      axis.text.x = ggplot2::element_text(size = 10),
      axis.title.x = ggplot2::element_text(size = 10),
      # y-axis
      axis.text.y = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_text(size = 10, angle = 90),
      # legend
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 12),
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

  #  Decide what to do with plot based on user-defined options.
  handlePlot(
    plot = adjusted_plot[[1]], prefix = prefix, return_plot = return_plot,
    filename = adjusted_plot[[2]], outdir = outdir, fig_width = adjusted_plot[[3]],
    fig_height = adjusted_plot[[4]]
  )
}
