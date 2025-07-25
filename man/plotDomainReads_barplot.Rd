% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-domain-reads.R
\name{plotDomainReads_barplot}
\alias{plotDomainReads_barplot}
\title{Plot numbers of classified reads per domain (bar plot)}
\usage{
plotDomainReads_barplot(
  report,
  include_sample_names = FALSE,
  orientation = "vertical",
  include_eukaryotes = FALSE,
  return_plot,
  outdir,
  prefix
)
}
\arguments{
\item{report}{A report generated by mergeReports() or loaded with load_STDreports().}

\item{include_sample_names}{Whether to include sample names in the plot.}

\item{orientation}{Whether plot should be horizontally or vertically oriented.}

\item{include_eukaryotes}{Whether to include eukaryotes in the plot.}

\item{return_plot}{Whether plot should be returned.}

\item{outdir}{Output directory where the plot should be saved.}

\item{prefix}{Prefix to be added to output plot name.}
}
\value{
If return_plot = TRUE, returns a bar plot with the reads classified under each domain.
}
\description{
This function takes a report generated by mergeReports() or loaded with load_STDreports()
and creates a bar plot showing the number of classified reads per domain per sample.
}
