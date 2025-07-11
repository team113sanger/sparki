% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-read-distribution.R
\name{plotDistribution_violin}
\alias{plotDistribution_violin}
\title{Plot the distribution of reads classified under each domain (violin plot)}
\usage{
plotDistribution_violin(report, return_plot, outdir, prefix)
}
\arguments{
\item{report}{A report generated by mergeReports().}

\item{return_plot}{Whether plot should be returned.}

\item{outdir}{Output directory where the plot should be saved.}

\item{prefix}{Prefix to be added to output file name.}
}
\value{
If return_plot = TRUE, returns the violin plot.
}
\description{
Plot the distribution of reads classified under each domain (violin plot)
}
