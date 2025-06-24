#' Save plot to output file
#'
#' @param plot A ggplot2 object.
#' @param filename The file name for the plot.
#' @param outdir Output directory where the plot should be saved.
#' @param fig_width The width of the output PDF file.
#' @param fig_height The height of the output PDF file.
#'
exportPlot <- function(plot, filename, outdir, fig_width, fig_height) {
  if (missing(fig_width) && missing(fig_height)) {
    fig_width <- 10
    fig_height <- 10
  }

  ggpubr::ggexport(
    plot,
    filename = paste0(outdir, filename),
    width = fig_width,
    height = fig_height,
    verbose = FALSE
  )
}

#' Add prefix to output file name
#'
#' @param filename The file name for the plot.
#' @param prefix The prefix that should be added to the file name.
#'
#' @return An updated version of the file name provided.
#'
handlePrefix <- function(filename, prefix) {
  if (prefix != "" && substr(prefix, nchar(prefix), nchar(prefix)) != "_") {
    prefix <- paste0(prefix, "_")
  }

  filename <- paste0(prefix, filename)

  return(filename)
}

#' Handle plot based on the user-defined arguments
#'
#' @param plot A ggplot2 object.
#' @param filename The file name for the plot.
#' @param prefix The prefix that should be added to the file name.
#' @param outdir Output directory where the plot should be saved.
#' @param fig_width The width of the output PDF file.
#' @param fig_height The height of the output PDF file.
#' @param return_plot Whether plot should be returned.
#'
#' @return Depending on the arguments provided, the plot will be returned.
#'
handlePlot <- function(plot, return_plot = FALSE, filename, prefix = "", outdir = NA, fig_width, fig_height) {
  #  If outdir is not NA, it means it has been provided by the user and therefore
  #  the plot should be saved to an output directory.
  if (!(is.na(outdir))) {
    #  Since the plot will be saved to an output directory, we need to deal with
    #  the file name and ensure a prefix is added in case that was specified by
    #  the user.
    filename <- handlePrefix(filename = filename, prefix = prefix)

    # Now we need to export the plot.
    exportPlot(
      plot = plot, filename = filename, outdir = outdir,
      fig_width = fig_width, fig_height = fig_height
    )
    #  If outdir is NA, it means the user has not provided an output directory and
    #  therefore the plot will not be saved.
  } else if (is.na(outdir)) {
    if (return_plot == TRUE) {
      return(plot)
    } else {
      print(plot)
    }
  }
}

#' Process a ggplot2 bar plot to define its orientation and dimensions
#'
#' @param plot A ggplot2 object.
#' @param n_samples The number of samples that will be displayed in the plot.
#' @param orientation The orientation of the plot ("horizontal"/"h" or "vertical"/"v").
#' @param factor A factor to adjust the dimensions of the plot.
#'
#' @return An updated version of the plot, with the appropriate orientation and
#'  dimensions.
#'
process_barplot_orientation <- function(plot, n_samples, orientation, factor) {
  #  proportion (x-axis) vs samples (y-axis)
  if (orientation %in% c("vertical", "v")) {
    pdf_width <- 5
    pdf_height <- determine_pdf_height(
      n_elements = n_samples,
      base_size = 3,
      factor = factor
    )

  #  samples (x-axis) vs proportion (y-axis)
  } else if (orientation %in% c("horizontal", "h")) {
    plot <- plot + ggplot2::coord_flip()

    pdf_height <- 3.5
    pdf_width <- determine_pdf_width(
      n_elements = n_samples,
      base_size = 4,
      factor = factor
    )
  } else {
    stop(
      "The orientation provided (", orientation, ") is not valid. ",
      "Please choose from the following options: 'horizontal' or 'h' ",
      "for horizontal plots; 'vertical' or 'v' for vertical plots."
    )
  }

  return(list(plot, pdf_width, pdf_height))
}

#' Process a ggplot2 bar plot to define its axis tick names
#'
#' @param plot A ggplot2 object.
#' @param orientation The orientation of the plot ("horizontal"/"h" or "vertical"/"v").
#' @param include_sample_names Whether sample names should be included in the plot.
#'
#' @return An updated version of the plot, with the appropriate tick names.
#'
process_barplot_ticknames <- function(plot, orientation, include_sample_names) {
  #  proportion (x-axis) vs samples (y-axis)
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

#' Adjust ggplot2 bar plot.
#'
#' @param plot A ggplot2 object.
#' @param n_samples The number of samples that will be displayed in the plot.
#' @param include_sample_names Whether sample names should be included in the plot.
#' @param orientation The orientation of the plot ("horizontal"/"h" or "vertical"/"v").
#' @param filename The file name for the plot.
#'
#' @return An updated version of the plot, with the appropriate tick names.
#'
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
    factor = pdf_factor
  )

  plot <- processed_plot[[1]]
  pdf_width <- processed_plot[[2]]
  pdf_height <- processed_plot[[3]]

  return(list(plot, filename, pdf_width, pdf_height))
}

#' DETERMINE WIDTH FOR PDF FILE
#'
#' @param n_elements Number of elements on the x-axis.
#' @param base_size Base plot width when there is only one element on the x-axis.
#' @param factor Value to help adjust the width.
#'
#' @return Plot width for PDF file.
#'
determine_pdf_width <- function(n_elements, base_size = 1, factor = 1) {
  #  Increment width if there are 2 elements or more.
  for (i in seq_len(n_elements)) {
    pdf_width <- base_size + (0.25 * factor)
  }

  return(pdf_width)
}

#' DETERMINE HEIGHT FOR PDF FILE
#'
#' @param n_elements Number of elements on the y-axis.
#' @param base_size Base plot height when there is only one element on the y-axis.
#' @param factor Value to help adjust the height.
#'
#' @return Plot height for PDF file.
#'
determine_pdf_height <- function(n_elements, base_size = 2, factor = 1) {
  #  Increment height if there are 2 elements or more.
  for (i in seq_len(n_elements)) {
    pdf_height <- base_size + (1 * factor)
  }

  return(pdf_height)
}

#' Parse a symbol-delimited list to get individual elements
#'
#' This function takes a symbol-delimited list of elements and splits it up into individual elements.
#' Symbols can be anything, e.g. ",", ";", "/", "//", "@", "." etc.
#'
#' @param del_list A list with symbol-delimited values; the symbol can be a comma, for example.
#' @param delimiter The symbol that delimits the elements in the list.
#'
#' @return A vector with individual elements.
#'
parse_delimited_list <- function(del_list, delimiter) {
  #  Split comma-separated list.
  elements <- unlist(stringr::str_split(del_list, delimiter))

  return(elements)
}

is_mpa <- function(report) {
  ifelse(
    COLNAME_MPA_TAXON_LEAF %in% colnames(report),
    return(TRUE),
    return(FALSE)
  )
}
