% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-classification-summary.R
\name{plotClassificationSummary_barplot}
\alias{plotClassificationSummary_barplot}
\title{Plot the numbers of classified/unclassified reads per sample (bar plot)}
\usage{
plotClassificationSummary_barplot(
  report,
  include_sample_names,
  orientation = "vertical",
  return_plot,
  outdir,
  prefix
)
}
\arguments{
\item{report}{Standard Kraken2 report.}

\item{include_sample_names}{Whether to include sample names in the plot.}

\item{orientation}{Whether plot should be horizontally or vertically oriented.}

\item{return_plot}{Whether plot should be returned.}

\item{outdir}{Output directory where the plot should be saved.}

\item{prefix}{Prefix to be added to output plot name.}
}
\value{
If return_plot = TRUE, returns a bar plot with the classification summary.
}
\description{
This function takes a report and creates a bar plot showing the number of classified and
unclassified reads per sample. Note that reports that were processed with subsetReports()
should not be used with this function!
}
