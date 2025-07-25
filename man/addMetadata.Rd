% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reports-add-info.R
\name{addMetadata}
\alias{addMetadata}
\title{Add metadata columns to report}
\usage{
addMetadata(report, metadata, metadata_sample_col, metadata_columns)
}
\arguments{
\item{report}{A report generated by load_STDreports(), load_MPAreports() or mergeReports().}

\item{metadata}{A metadata table loaded with loadMetadata().}

\item{metadata_sample_col}{The name of the column in the metadata table that contains sample IDs (these sample
IDs must match the sample IDs present in the report).}

\item{metadata_columns}{The names of the columns in the metadata table that should be added to the report.}
}
\value{
An updated version of the report, containing the metadata columns specified.
}
\description{
This function takes as input a report (loaded with load_STDreports() or load_MPAreports(), or generated by
mergeReports()) and a metadata table (loaded with loadMetadata()) and adds specified columns from the latter
to the former.
}
