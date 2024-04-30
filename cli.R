library(optparse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(cowplot)
library(factoextra, lib.loc = "/nfs/users/nfs_j/jb62/r_packages/")
library(ggpubr, lib.loc = "/nfs/users/nfs_j/jb62/r_packages/")

my_library <- "/lustre/scratch126/casm/team113da/users/jb62/projects/single-cell-and-spatial-profiling-of-skin-tumours/R/similarity_assessment/"
source(paste0(my_library, "main.R"))
source(paste0(my_library, "exceptions.R"))
source(paste0(my_library, "library/helper.R"))
source(paste0(my_library, "library/utilities.R"))
source(paste0(my_library, "library/plotting.R"))

###################
# Parse arguments #
###################

# Options that overlap between Mahalanobis analysis and
analysis_option_list_overlap <- list(
    make_option(
        c("-i", "--infile"),
        dest = "infile",
        action = "store",
        default = NA,
        type = "character",
        help = "Path to gene expression matrix."
    ),
    make_option(
        c("-o", "--outdir"),
        dest = "outdir",
        action = "store",
        default = NA,
        type = "character",
        help = "Path to output directory."
    ),
    make_option(
        c("-p", "--prefix"),
        dest = "prefix",
        action = "store",
        default = NA,
        type = "character",
        help = "Prefix of output files."
    ), 
    make_option(
        c("-m", "--metadata"),
        dest = "metadata",
        action = "store",
        default = NA,
        type = "character",
        help = "Path to metadata file. [Optional]"
    ),
    make_option(
        c("-s", "--sample-column"),
        dest = "samples",
        action = "store",
        default = NA,
        type = "character",
        help = "Column in metadata file that contains sample IDs. [Optional]"
    ),
    make_option(
        c("-c", "--cat-column"),
        dest = "categories",
        action = "store",
        default = NA,
        type = "character",
        help = "Column in metadata file that contains categories of interest. [Optional]"
    )
)

# Options that are Mahalanobis-specific.
analysis_option_list_mahalanobis <- list(
    make_option(
        c("-d", "--mode"),
        dest = "mode",
        action = "store",
        default = NA,
        type = "character",
        help = "Mode that should be used to assess the Mahalanobis 
                distances (either 'PCA' or 'genes'). [Default = PCA]"
    ),
    make_option(
        c("-n", "--n-pcs"),
        dest = "npcs",
        action = "store",
        default = NA,
        type = "numeric",
        help = "If mode is PCA, users can define the number of principal 
                components (PCs) that should be taken into account in the 
                Mahalanobis distance assessment."
    ),
    make_option(
        c("-a", "--columns-to-plot"),
        dest = "cols_to_plot",
        action = "store",
        default = NA,
        type = "character",
        help = "Comma-separated list of metadata columns whose info should be 
                included in the plots. [Optional]"
    ),
    make_option(
        c("-v", "--drivers-column"),
        dest = "driver_genes",
        action = "store",
        default = NA,
        type = "character",
        help = "Column in metadata file that contains driver genes delimited by a 
                symbol (e.g. comma, slash etc). [Optional]"
    ),
    make_option(
        c("-l", "--drivers-delimiter"),
        dest = "driver_genes_del",
        action = "store",
        default = NA,
        type = "character",
        help = "If --driver-column was provided, the symbol delimiting the driver 
                genes must also be provided. [Optional]"
    )
)

# Options that are PCA-specific.
analysis_option_list_pca <- list()

####################
# Helper functions #
####################

# 
read_prog_name_env_var <- function() {
    # Read the environment variable that contains the name of the program.
    prog_name <- Sys.getenv("PROG_NAME")
    if (is.null(prog_name)) {
        stop("Environment variable PROG_NAME not set.")
    }
    return(prog_name)
}

# Parse input arguments.
parse_options <- function(option_list, description = "") {
    # We need to fetch the program name from the environment variables.
    prog_name <- read_prog_name_env_var()

    # CLI arguments need some pre-processing before they can be parsed.
    raw_arguments <- commandArgs(trailingOnly = TRUE) # From the CLI
    raw_arguments <- raw_arguments[-1] # We only want to handle the arguments that come after the function name which is the first element.

    # Parse the arguments.
    parser <- OptionParser(
        option_list = option_list,
        add_help_option = TRUE,
        prog = prog_name,
        description = description
    )
    arguments <- parse_args(parser, args = raw_arguments)
    return(arguments)
}

# Assess Mahalanobis distances.
run_mahalanobis_analysis <- function(arguments) {
    mahalanobis_analysis(
        gex_data = arguments$options$infile,
        outdir = arguments$options$outdir,
        prefix = arguments$options$prefix,
        mode = arguments$options$mode,
        npcs = arguments$options$npcs,
        metadata = arguments$options$metadata,
        sample_column = arguments$options$samples,
        cat_column = arguments$options$categories,
        columns_for_plot = arguments$options$cols_to_plot,
        drivers_column = arguments$options$driver_genes,
        drivers_del = arguments$options$driver_genes_del
    )
}

# Compute PCA.
run_pca_analysis <- function(arguments) {
    pca_analysis(
        gex_data = arguments$options$infile,
        outdir = arguments$options$outdir,
        prefix = arguments$options$prefix,
        metadata = arguments$options$metadata,
        sample_column = arguments$options$samples,
        cat_column = arguments$options$categories
    )
}

###################
# CLI entrypoints #
###################

cli_example_function <- function() {
    # Minimal function that parses the option --msg and prints it to stderr
    option_list <- list(
        make_option(
            c("-m", "--msg"),
            dest = "msg",
            action = "store",
            default = "Default message.",
            type = "character",
            help = "Message to print."
        )
    )
    description <- "Example function that prints a message to stdout and stderr."
    parsed_arguments <- parse_options(option_list, description)
    msg <- parsed_arguments$msg
    cat("Stdout message:", msg, "\n")
    message("Stderr message: ", msg)
}

cli_pca_analysis <- function() {

    analysis_option_list <- c(analysis_option_list_overlap, analysis_option_list_pca)
    description <- "Perform principal component analysis (PCA) on a gene expression matrix."

    parsed_arguments <- parse_options(analysis_option_list, description)
    run_pca_analysis(parsed_arguments)
}

cli_mahalanobis_analysis <- function() {

    analysis_option_list <- c(analysis_option_list_overlap, analysis_option_list_mahalanobis)
    description <- "Assess the Mahalanobis distance of each sample in a cohort based on their gene expression patterns."

    parsed_arguments <- parse_options(analysis_option_list, description)
    run_mahalanobis_analysis(parsed_arguments)
}

cli_run_multiple_analyses <- function() {

    analysis_option_list <- c(analysis_option_list_overlap, analysis_option_list_mahalanobis, analysis_option_list_pca)
    description <- "Carry out a thorough similarity assessment on samples of a given cohort based on their gene expression patterns.\nThe similarity assessment involves (1) Mahalanobis distance calculation and (2) principal component analysis (PCA)."

    parsed_arguments <- parse_options(analysis_option_list, description)
    run_pca_analysis(parsed_arguments)
    run_mahalanobis_analysis(parsed_arguments)
}

######################
# Run CLI Entrypoint #
######################

args <- commandArgs(trailingOnly = TRUE)[1]
func_name <- args[1]

# Empty list is to satisfy the function signature, the actual args are passed via command line.
do.call(func_name, list()) 

