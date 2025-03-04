#' Subset report in order to filter out a given species.
#'
#' @export
subsetReports <- function(report, species) {
  species <- check_species_in_report(report, species)

  # Determine the genus and family of an user-defined species.
  genus <- get_name_at_rank(report, species, "G")
  family <- get_name_at_rank(report, species, "F")

  # Create vector with species, genus and family to be filtered out.
  taxa_to_remove <- c(species, genus, family)

  logger::log_info(paste0("Filtering out ", species, ", ", genus, ", and ", family, "..."))

  report <- report |>
    #  Select rows corresponding to species, genus and family.
    dplyr::filter(!!as.name(COLNAME_STD_RANK) %in% c("F", "G", "S")) |>
    #  Filter out the user-defined organism.
    dplyr::filter(!(!!as.name(COLNAME_STD_TAXON) %in% taxa_to_remove))

  return(report)
}

get_name_at_rank <- function(report, species, rank) {
  # Get the species position (row number) in the report (first appearance only).
  species_position <- which(report[, COLNAME_STD_TAXON] == species)[1]

  #  Find the corresponding rank position (row number) in the report (first appearance only).
  i <- 1
  while (i > 0) {
    #  If the desired rank has been found, break the loop.
    if (report[, COLNAME_STD_RANK, drop = TRUE][(species_position - i)] == rank) {
      rank_position <- (species_position - i)
      break
    } else {
      i <- i + 1
    }
  }

  subset_report <- report |> dplyr::slice(rank_position)

  return(subset_report[[COLNAME_STD_TAXON]])
}
