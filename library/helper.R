########################################
## HELPER FUNCTIONS FOR READING FILES ##
#######################################################################################################

#' LOAD INFORMATION FROM KRAKEN2'S REFERENCE DATABASE
#' 
#' This function takes a path to an inspect.txt file inside a Kraken2 reference database
#' and reads the table. The first 7 lines are skipped, as they correspond to header lines.
#' Column names are added before the reference database dataframe is returned.
#' 
#' @param reference_path Path to an inspect.txt file inside a Kraken2 reference database.
#' @return 
#' 
load_reference <- function(reference_path) {

    # Read Kraken2 reference file (inspect.txt).
    ref <- read.table(
        reference_path,
        row.names = NULL, 
        header = FALSE,
        sep = "\t",
        skip = 7 # Skip header lines.
    )
    # Add column names.
    colnames(ref) <- c(
        "percent_minimisers_clade", 
        "number_minimisers_clade",
        "number_minimisers_taxon",
        "rank", 
        "taxid", 
        "scientific_name"
    )

    return(ref)

}

#' LOAD METADATA TABLE
#' 
#' This function takes a path to a metadata file and reads the table.
#' and reads the table. The first 7 lines are skipped, as they correspond to header lines.
#' Column names are added before the reference database dataframe is returned.
#' 
#' @param mdata_path Path to a metadata table.
#' @return 
#' 
load_metadata <- function(mdata_path) {

    # Read metadata file.
    mdata <- read.csv(
        mdata_path,
        row.names = 1,
        header = TRUE
    )

    return(mdata)
}

#' PREPARE DATAFRAME WITH CONCATENATED MPA-STYLE REPORTS 
#' 
#' This function takes a path to a metadata file and reads the table.
#' and reads the table. The first 7 lines are skipped, as they correspond to header lines.
#' Column names are added before the reference database dataframe is returned.
#' 
#' @param mdata_path Path to a metadata table.
#' @return 
#' 
load_mpa <- function(mpa_dir, verbose = TRUE) {

    # Get paths to MPA-style reports in a specified directory.
    mpa_file_names <- list.files(mpa_dir, pattern = "mpa$")

    # Check if directory really has any MPA-style reports...
    if (length(mpa_file_names) == 0) {
        stop(paste0("No MPA-style reports were found at ", mpa_dir, ". Please review your input."))
    
    # If it does, read MPA-style reports.
    } else {

        # Create empty dataframe to store MPA-style reports.
        merged_mpa <- data.frame(matrix(nrow = 0, ncol = 2))

        # Iterate over files in specified MPA directory (ideally each file should be a different sample)...
        for (file_name in mpa_file_names) {

            # Get sample name.
            sample_name <- gsub("\\..*", "", basename(file_name))

            if (verbose == TRUE) message(paste0("Reading MPA-style report for sample ", sample_name))

            # Read MPA-style report.
            mpa <- read.csv(file = paste0(mpa_dir, "/", file_name), header = FALSE, sep = "\t")

            # Add column names.
            colnames(mpa) <- c(
                "classification", 
                "number_fragments_clade"
            )

            # Add sample name to MPA dataframe.
            mpa$sample <- rep(sample_name, times = nrow(mpa))

            # Reorder columns.
            mpa <- mpa[, c(
                "sample", 
                "classification", 
                "number_fragments_clade"
            )]

            # Add sample results to dataframe where all results will be stored.
            merged_mpa <- rbind(merged_mpa, mpa)
        } 

        return(merged_mpa) 
    }
}

load_reports <- function(reports_dir, verbose = TRUE) {

    # Get paths to reports in a specified directory.
    reports_file_names <- list.files(reports_dir, pattern = "kraken$")

    # Check if directory really has any reports...
    if (length(reports_file_names) == 0) {
        stop(paste0("No reports were found at ", reports_dir, ". Please review your input."))
    
    # If it does, read reports.
    } else {

        # Create empty dataframe to store reports.
        merged_reports <- data.frame(matrix(nrow = 0, ncol = 2))

        # Iterate over files in specified reports directory (ideally each file should be a different sample)...
        for (file_name in reports_file_names) {

            # Get sample name.
            sample_name <- gsub("\\..*", "", basename(file_name))               
            
            if (verbose == TRUE) message(paste0("Reading standard report for sample ", sample_name))

            # Read report.
            report <- read.csv(file = paste0(reports_dir, "/", file_name), header = FALSE, sep = "\t")
            
            # Add column names.
            colnames(report) <- c(
                "percent_fragments_clade", 
                "number_fragments_clade", 
                "number_fragments_taxon", 
                "number_minimisers", 
                "number_distinct_minimisers",
                "rank", 
                "taxid", 
                "scientific_name"
            )

            # Add sample name to report dataframe.
            report$sample <- rep(sample_name, times = nrow(report))

            # Reorder columns.
            report <- report[, c(
                "sample",
                "percent_fragments_clade", 
                "number_fragments_clade", 
                "number_fragments_taxon", 
                "number_minimisers", 
                "number_distinct_minimisers",
                "rank", 
                "taxid", 
                "scientific_name"
            )]

            # Remove extra spaces that are present before the scientific names.
            report$scientific_name <- gsub("^[[:space:]]\\s*(.*?)", "", report$scientific_name, perl = TRUE)

            # Add sample results to dataframe where all results will be stored.
            merged_reports <- rbind(merged_reports, report)
        } 

        return(merged_reports) 
    }
}

############################################
## HELPER FUNCTIONS FOR MPA-STYLE REPORTS ##
#######################################################################################################

mpa_determine_last_rank <- function(merged_mpa) {

    merged_mpa$last_rank_in_line <- sapply(seq_len(nrow(merged_mpa)), function(x) {

        # Get scientific name (in this case it will have a prefix indicating the taxonomic rank).
        last_rank <- tail(unlist(stringr::str_split(merged_mpa$classification[x], "\\|")), n = 1)

        # Identify which taxonomic rank the scientific name belongs to.
        if(grepl("s__", last_rank)) return("species")
        if(grepl("g__", last_rank)) return("genus")
        if(grepl("f__", last_rank)) return("family")
        if(grepl("o__", last_rank)) return("order")
        if(grepl("c__", last_rank)) return("class")
        if(grepl("p__", last_rank)) return("phylum")
        if(grepl("k__", last_rank)) return("kingdom")
        if(grepl("d__", last_rank)) return("domain")
    })

    return(merged_mpa)
}

extract_domain <- function(mpa_line, last_rank = TRUE) {

    if (last_rank) { 
        return(stringr::str_match(mpa_line, "^d__(.+)")[2])
    } else {
        return(stringr::str_match(mpa_line, "^d__\\s*(.*?)\\s*[|](.*?)")[2])
    }
}

extract_family <- function(mpa_line, last_rank = TRUE) {
    
    if (last_rank) { 
        return(stringr::str_match(mpa_line, "f__(.+)")[2])
    } else {
        return(stringr::str_match(mpa_line, "(.*?)[|]\\s*f__(.*?)\\s*[|](.*?)")[3])
    }
}

extract_genus <- function(mpa_line, last_rank = TRUE) {
    
    if (last_rank) { 
        return(stringr::str_match(mpa_line, "g__(.+)")[2])
    } else {
        return(stringr::str_match(mpa_line, "(.*?)[|]\\s*g__(.*?)\\s*[|](.*?)")[3])
    }
}

extract_species <- function(mpa_line) {
    
    return(stringr::str_match(mpa_line, "s__(.+)")[2])
}

mpa_add_rank_columns <- function(merged_mpa) {
    
    merged_mpa$domain <- sapply(seq_len(nrow(merged_mpa)), function (x) {

        if (merged_mpa$last_rank_in_line[x] == "domain") {
            return(extract_domain(merged_mpa$classification[x], last_rank = TRUE))
        } else {
            return(extract_domain(merged_mpa$classification[x], last_rank = FALSE))
        }
    })

    merged_mpa$family <- sapply(seq_len(nrow(merged_mpa)), function (x) {

        if (merged_mpa$last_rank_in_line[x] == "family") {
            return(extract_family(merged_mpa$classification[x], last_rank = TRUE))
        } else if (merged_mpa$last_rank_in_line[x] %in% c("kingdom", "phylum", "class", "order")) {
            return(NA)
        } else {
            return(extract_family(merged_mpa$classification[x], last_rank = FALSE))
        }
    })

    merged_mpa$genus <- sapply(seq_len(nrow(merged_mpa)), function (x) {

        if (merged_mpa$last_rank_in_line[x] == "genus") {
            return(extract_genus(merged_mpa$classification[x], last_rank = TRUE))
        } else if (merged_mpa$last_rank_in_line[x] %in% c("kingdom", "phylum", "class", "order")) {
            return(NA)
        } else {
            return(extract_genus(merged_mpa$classification[x], last_rank = FALSE))
        }
    })

    merged_mpa$species <- sapply(seq_len(nrow(merged_mpa)), function (x) {

        if (merged_mpa$last_rank_in_line[x] == "species") {
            return(extract_species(merged_mpa$classification[x]))
        } else if (merged_mpa$last_rank_in_line[x] %in% c("kingdom", "phylum", "class", "order")) {
            return(NA)
        } else {
            return(extract_species(merged_mpa$classification[x]))
        }
    })

    return(merged_mpa)
}

mpa_get_n_reads_in_samples <- function(merged_mpa) {

    n_reads_in_samples <- data.frame(matrix(nrow = 0, ncol = 4))

    # Iterate over samples in merged MPA dataframe.
    for (sample in unique(merged_mpa$sample)) {

        # Subset MPA dataframe to keep only a given sample.
        merged_mpa_sample <- merged_mpa[merged_mpa$sample == sample,]

        # Put sample-level results of all domains (4 domains) or non-eukaryote domains (3 domains) in 2 different dataframes.
        merged_mpa_sample_4_domains <- merged_mpa_sample[grep("d__Eukaryota$|d__Viruses$|d__Archaea$|d__Bacteria$", merged_mpa_sample$classification),]
        merged_mpa_sample_3_domains <- merged_mpa_sample[grep("d__Viruses$|d__Archaea$|d__Bacteria$", merged_mpa_sample$classification),]
        
        # Sum all sample-level reads classified under all domains.
        sum_4_domains <- sum(merged_mpa_sample_4_domains$number_fragments_clade)
        n_reads_in_samples <- rbind(n_reads_in_samples, c(sample, "all", sum_4_domains))
        
        # Sum all sample-level reads classified under non-eukaryote domains.
        sum_3_domains <- sum(merged_mpa_sample_3_domains$number_fragments_clade)
        n_reads_in_samples <- rbind(n_reads_in_samples, c(sample, "subset", sum_3_domains))
    }

    colnames(n_reads_in_samples) <- c("sample", "domains_considered", "number_reads_clade")
    n_reads_in_samples$number_reads_clade <- as.numeric(n_reads_in_samples$number_reads_clade)

    return(n_reads_in_samples)
}

mpa_get_n_ranks_in_samples <- function(merged_mpa) {

    n_ranks_in_samples <- data.frame(matrix(nrow = 0, ncol = 4))

    # Iterate over samples...
    for (sample in unique(merged_mpa$sample)) {

        # Iterate over domains...
        for (domain in unique(merged_mpa$domain)) {

            # Get number of species, genuses and families identified in a given sample.
            n_species <- length(unique(merged_mpa$species[merged_mpa$sample == sample & merged_mpa$domain == domain]))
            n_genus <- length(unique(merged_mpa$genus[merged_mpa$sample == sample & merged_mpa$domain == domain]))
            n_family <- length(unique(merged_mpa$family[merged_mpa$sample == sample & merged_mpa$domain == domain]))

            # Add numbers to dataframe.
            n_ranks_in_samples <- rbind(n_ranks_in_samples, c(sample, domain, "species", n_species))
            n_ranks_in_samples <- rbind(n_ranks_in_samples, c(sample, domain, "genus", n_genus))
            n_ranks_in_samples <- rbind(n_ranks_in_samples, c(sample, domain, "family", n_family))
        }
    }
    colnames(n_ranks_in_samples) <- c("sample", "domain", "rank", "value")

    return(n_ranks_in_samples)
}

#subset_mpa <- function(merged_mpa, include_human = TRUE) {

    # Subsetting family-, genus- and species-level results
#    merged_mpa <- merged_mpa[merged_mpa$last_rank_in_line %in% c("family", "genus", "species"),]

#    if (include_human) {
#        return(merged_mpa)
#    } else {
#        merged_mpa <- merged_mpa[!(merged_mpa$family %in% c("Hominidae", "Homo", "Homo sapiens")), ]
#        return(merged_mpa)
#    }
#}

##################################
## HELPER FUNCTIONS FOR REPORTS ##
#######################################################################################################

get_class_unclass_numbers <- function(merged_reports) {

    class_unclass_df <- data.frame(matrix(nrow = 0, ncol = 4))

    for(sample in unique(merged_reports$sample)) {

        unclass_value <- as.numeric(merged_reports$number_fragments_clade[merged_reports$scientific_name == "unclassified" & merged_reports$sample == sample])
        class_value <- as.numeric(merged_reports$number_fragments_clade[merged_reports$scientific_name == "root" & merged_reports$sample == sample])

        class_unclass_df <- rbind(class_unclass_df, c(sample, "unclassified", unclass_value))
        class_unclass_df <- rbind(class_unclass_df, c(sample, "classified", class_value))
    }

    colnames(class_unclass_df) <- c("sample", "type", "value")
    class_unclass_df$value <- as.numeric(class_unclass_df$value)

    return(class_unclass_df)
}

subset_reports <- function(merged_reports, include_human = TRUE) {

    # Subsetting family-, genus- and species-level results
    merged_reports <- merged_reports[merged_reports$rank %in% c("F", "G", "S"),]

    if (include_human) {
        return(merged_reports)
    } else {
        merged_reports <- merged_reports[!(merged_reports$scientific_name %in% c("Hominidae", "Homo", "Homo sapiens")),]
        return(merged_reports)
    }
}

report_add_domains <- function(merged_reports, merged_mpa) {

    merged_reports$domain <- sapply(seq_len(nrow(merged_reports)), function(x) {

        domain <- NA

        if (merged_reports$rank[x] == "F") {
            domain <- unique(merged_mpa$domain[which(merged_mpa$family == merged_reports$scientific_name[x])])
        } else if (merged_reports$rank[x] == "G") {
            domain <- unique(merged_mpa$domain[which(merged_mpa$genus == merged_reports$scientific_name[x])])
        } else if (merged_reports$rank[x] == "S") {
            domain <- unique(merged_mpa$domain[which(merged_mpa$species == merged_reports$scientific_name[x])])
        }

        return(domain)
    })

    return(merged_reports)
}





#' DETERMINE WIDTH FOR PDF FILE
#' 
#' This function takes a number of elements (e.g. number of samples) and determines the width
#' that should be used for a plot to be saved in PDF format. This function is tailored for the
#' plots generated by plot_mahalanobis_with_metadata() and plot_mahalanobis_with_driver_genes(),
#' so should be used with caution for any other ends.
#' 
#' @param n_elements Number of elements on the x-axis.
#' @param factor Value to help adjust the width.
#' @return Plot width for PDF file.
#' @export
determine_pdf_width <- function(n_elements, factor = 1) {

    # Basic plot width when there is only one element on the x-axis.
    pdf_width <- 1

    # Increment width if there are 2 elements or more.
    for (i in 1:n_elements) {
        pdf_width <- pdf_width + (0.25 * factor)
    }

    return(pdf_width)
}

#' DETERMINE HEIGHT FOR PDF FILE
#' 
#' This function takes a number of elements (e.g. number of genes) and determines the height
#' that should be used for a plot to be saved in PDF format. This function is tailored for the
#' plots generated by plot_mahalanobis_with_metadata() and plot_mahalanobis_with_driver_genes(),
#' so should be used with caution for any other ends.
#' 
#' @param n_elements Number of elements on the y-axis.
#' @param factor Value to help adjust the height.
#' @return Plot height for PDF file.
#' @export
determine_pdf_height <- function(n_elements, factor = 1) {

    # Basic plot height when there is only one element on the y-axis.
    pdf_height <- 2

    # Increment height if there are 2 elements or more.
    for (i in 1:n_elements) {
        pdf_height <- pdf_height + (1 * factor)
    }

    return(pdf_height)
}
