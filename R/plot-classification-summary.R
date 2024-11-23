getClassificationSummary <- function(report) {
  if (is_mpa(report)) {
    stop(paste0(
      "This function is not applicable to MPA-style reports. The purpose of this function ",
      "is to assess the number of classified and unclassified reads for each sample, but the ",
      "MPA-style reports do not contain information on unclassified reads."
    ))
  }

  summary <- report |>
    #  Select columns of interest.
    dplyr::select(!!COLNAME_STD_SAMPLE, !!COLNAME_STD_TAXON, !!COLNAME_STD_N_FRAG_CLADE) |>
    # Select rows of interest.
    dplyr::filter(!!as.name(COLNAME_STD_TAXON) %in% c("unclassified", "root")) |>
    #  Rename column.
    dplyr::rename(
      !!COLNAME_CLASSIF_SUMMARY_SAMPLE := !!COLNAME_STD_SAMPLE,
      !!COLNAME_CLASSIF_SUMMARY_READ_TYPE := !!COLNAME_STD_TAXON,
      !!COLNAME_CLASSIF_SUMMARY_N_FRAG := !!COLNAME_STD_N_FRAG_CLADE
    ) |>
    #  Rename elements in column.
    dplyr::mutate(!!COLNAME_CLASSIF_SUMMARY_READ_TYPE := dplyr::case_when(
      !!as.name(COLNAME_CLASSIF_SUMMARY_READ_TYPE) == "unclassified" ~ "Unclassified",
      !!as.name(COLNAME_CLASSIF_SUMMARY_READ_TYPE) == "root" ~ "Classified"
    )) |>
    #  Create column with log-transformed number of clade-level fragments.
    dplyr::mutate(
      !!COLNAME_CLASSIF_SUMMARY_LOG_N_FRAG := log10(!!as.name(COLNAME_CLASSIF_SUMMARY_N_FRAG))
    ) |>
    # Reorder columns.
    dplyr::relocate(
      !!COLNAME_CLASSIF_SUMMARY_LOG_N_FRAG,
      .after = !!COLNAME_CLASSIF_SUMMARY_N_FRAG
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
    #  Select relevant rows.
    dplyr::filter(!!as.name(COLNAME_STD_TAXON) == "root") |>
    # Create column with proportion of classified reads.
    dplyr::mutate(
      !!COLNAME_PROP_CLASSIFIED := (!!as.name(COLNAME_STD_N_FRAG_CLADE) / !!as.name(COLNAME_STD_SAMPLE_SIZE))
    ) |>
    #  Remove columns.
    dplyr::select(!c(!!COLNAME_STD_SAMPLE_SIZE, !!COLNAME_STD_TAXON, !!COLNAME_STD_N_FRAG_CLADE)) |>
    # Rename columns.
    dplyr::rename(!!COLNAME_PROP_SAMPLE := !!COLNAME_STD_SAMPLE)

  return(proportion)
}

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

  #  Assign NA to outdir in case it has not been provided by the user.
  if (missing(outdir)) outdir <- NA

  # Prepare data for plotting.
  summary <- getClassificationSummary(report)

  #  Create violin plot.
  plot <- ggplot2::ggplot(
    summary,
    ggplot2::aes(
      x = get(COLNAME_CLASSIF_SUMMARY_READ_TYPE),
      y = get(COLNAME_CLASSIF_SUMMARY_LOG_N_FRAG)
    )
  ) +
    ggplot2::geom_violin(scale = "width", fill = "white", color = "black") +
    ggplot2::geom_line(ggplot2::aes(group = get(COLNAME_CLASSIF_SUMMARY_SAMPLE)), alpha = 0.25) +
    ggplot2::geom_point(ggplot2::aes(color = get(COLNAME_CLASSIF_SUMMARY_READ_TYPE)), alpha = 0.5) +
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
    ggplot2::ylab(expression("log"[10] ~ "(# reads)")) +
    ggplot2::scale_color_manual(values = c("indianred2", "royalblue"))

  #  Decide what to do with plot based on user-defined options.
  handlePlot(
    plot = plot, prefix = prefix, return_plot = return_plot, filename = paste0(PLOT_CLASSIF_SUMMARY_VIOLIN, ".pdf"),
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
    return_plot, outdir, prefix = "") {
  # Stop if report provided is in MPA-style format.
  if (is_mpa(report)) stop("This function is not applicable to MPA-style reports.")

  #  Assign NA to outdir in case it has not been provided by the user.
  if (missing(outdir)) outdir <- NA

  # Prepare data for plotting.
  summary <- getClassificationSummary(report)

  #  Create bar plot.
  plot <- ggplot2::ggplot(
    summary,
    ggplot2::aes(
      x = get(COLNAME_CLASSIF_SUMMARY_N_FRAG),
      y = get(COLNAME_CLASSIF_SUMMARY_SAMPLE),
      fill = get(COLNAME_CLASSIF_SUMMARY_READ_TYPE)
    )
  ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      #  x-axis
      axis.text.x = ggplot2::element_text(size = 12, vjust = 0.5),
      axis.title.x = ggplot2::element_text(size = 14),
      #  y-axis
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
    filename = PLOT_CLASSIF_SUMMARY_BARPLOT
  )

  #  Decide what to do with plot based on user-defined options.
  handlePlot(
    plot = adjusted_plot[[1]], prefix = prefix, return_plot = return_plot,
    filename = adjusted_plot[[2]], outdir = outdir, fig_width = adjusted_plot[[3]],
    fig_height = adjusted_plot[[4]]
  )
}

plotClassificationProportion <- function(report, return_plot, outdir, prefix = "") {
  # Stop if report provided is in MPA-style format.
  if (is_mpa(report)) stop("This function is not applicable to MPA-style reports.")

  #  Assign NA to outdir in case it has not been provided by the user.
  if (missing(outdir)) outdir <- NA

  # Prepare data for plotting.
  summary <- getClassificationProportion(report)

  #  Create violin plot.
  plot <- ggplot2::ggplot(
    summary,
    ggplot2::aes(x = "", y = get(COLNAME_PROP_CLASSIFIED))
  ) +
    ggplot2::geom_violin(scale = "width", fill = "palegreen1", color = "black") +
    ggplot2::geom_jitter(color = "black") +
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
    ggplot2::ylab("Proportion of\nclassified reads")

  #  Decide what to do with plot based on user-defined options.
  handlePlot(
    plot = plot, prefix = prefix, return_plot = return_plot, filename = paste0(PLOT_CLASSIF_PROPORTION, ".pdf"),
    outdir = outdir, fig_width = 3, fig_height = 4
  )
}
