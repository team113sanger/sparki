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

check_species_in_report <- function(report, species) {
  # Here the species name will be specified by the user, so it might be that it will not
  # be present in the report exactly as provided by the user.

  # If the species name specified by the user is not present in the report as is...
  if (species %in% report[[COLNAME_STD_TAXON]]) {

    logger::log_debug(paste0(
      "Species ", species, " was found in the report at rows ", which(report[[COLNAME_STD_TAXON]] == species), "."
    ))

    return(species)

  } else if (!(species %in% report[[COLNAME_STD_TAXON]])) {

    logger::log_debug(paste0("Species ", glue::single_quote(species), " was not found in the report."))

    # Check if the species name specified by the user has an underscore...
    if (grepl("_", species)) {

      logger::log_debug(paste0("Species ", glue::single_quote(species), " is delimited with an underscore."))

      # Then replace the underscore with a white space and look for the species name
      # in the report...
      species_alt <- gsub("_", " ", species)

      logger::log_debug(paste0(
        "The white space-delimited version of ", glue::single_quote(species), " is ", glue::single_quote(species_alt),
        "; checking if it is present in the report..."
      ))

      if (species_alt %in% report[[COLNAME_STD_TAXON]]) {

        logger::log_debug(paste0(
          "The white space-delimited version of ", glue::single_quote(species), " (", glue::single_quote(species_alt),
          ") was found in the report!"
        ))

        logger::log_warn(paste0(
          glue::single_quote(species), " was not found in the report, but ", glue::single_quote(species_alt), " was."
        ))

        # Return the "corrected" species name ("_" -> " ").
        return(species_alt)
      }

    # Alternatively, check if the species name specified by the user has a white space...
    } else if (grepl("[[:space:]]", species)) {

      logger::log_debug(paste0("Species ", glue::single_quote(species), " is delimited with a white space."))

      # Then replace the white space with an underscore and look for the species name
      # in the report...
      species_alt <- gsub(" ", "_", species)

    logger::log_debug(paste0(
        "The underscore-delimited version of ", glue::single_quote(species), " is ", glue::single_quote(species_alt),
        "; checking if it is present in the report..."
      ))

      if (species_alt %in% report[[COLNAME_STD_TAXON]]) {

        logger::log_debug(paste0(
          "The underscore-delimited version of ", glue::single_quote(species), " (", glue::single_quote(species_alt),
          ") was found in the report!"
        ))

        logger::log_warn(paste0(
          glue::single_quote(species), " was not found in the report, but ", glue::single_quote(species_alt), " was."
        ))

        # Return the "corrected" species name (" " -> "_").
        return(species_alt)
      }

    # If none of the above works, raise an error as the species name specified was not
    # found at all in the report (bearing in mind that the delimiters supported here
    # are "_" [underscore] and " " [white space]).
    } else {
      stop(paste0(
        "The species ", glue::single_quote(species), " is not present in the dataframe. Please review your input!"
      ))
    }
  }
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
