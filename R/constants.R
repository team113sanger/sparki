######################
## HELPER CONSTANTS ##
#####################################################################

#### Reference database ####
COLNAME_REF_DB_PCT_FRAG_CLADE <- "pct_fragments_clade"
COLNAME_REF_DB_MINIMISERS_CLADE <- "n_minimisers_clade"
COLNAME_REF_DB_MINIMISERS_TAXON <- "n_minimisers_taxon"
COLNAME_REF_DB_RANK <- "rank"
COLNAME_REF_DB_NCBI_ID <- "ncbi_id"
COLNAME_REF_DB_TAXON <- "taxon"

#### Standard report ####
COLNAME_STD_SAMPLE <- "sample"
COLNAME_STD_PCT_FRAG_CLADE <- "pct_fragments_clade"
COLNAME_STD_N_FRAG_CLADE <- "n_fragments_clade"
COLNAME_STD_N_FRAG_TAXON <- "n_fragments_taxon"
COLNAME_STD_MINIMISERS <- "n_minimisers"
COLNAME_STD_UNIQ_MINIMISERS <- "n_distinct_minimisers"
COLNAME_STD_RANK <- "rank"
COLNAME_STD_NCBI_ID <- "ncbi_id"
COLNAME_STD_TAXON <- "taxon"
COLNAME_STD_TOTAL_READS <- "n_total_reads_in_sample"
COLNAME_STD_DB_MINIMISERS_CLADE <- "db_number_minimisers_clade"
COLNAME_STD_DB_MINIMISERS_TAXON <- "db_number_minimisers_taxon"

#### MPA-style report #### 
COLNAME_MPA_SAMPLE <- "sample"
COLNAME_MPA_N_FRAG_CLADE <- "n_fragments_clade"
COLNAME_MPA_TAXON <- "taxon_hierarchy"
COLNAME_MPA_RANK <- "rank"
COLNAME_MPA_NCBI_ID <- "ncbi_id"
COLNAME_MPA_TOTAL_READS <- "n_total_reads_in_sample"
COLNAME_MPA_DB_MINIMISERS_CLADE <- "db_number_minimisers_clade"
COLNAME_MPA_DB_MINIMISERS_TAXON <- "db_number_minimisers_taxon"

#### Auxiliary dataframes ####

# Number of domain reads per sample.
COLNAME_DOMAIN_READS_SAMPLE <- "sample"
COLNAME_DOMAIN_READS_SCOPE <- "domains_considered"
COLNAME_DOMAIN_READS_N_READS <- "n_clade_reads"

# Number of classified/unclassified reads per sample.
COLNAME_CLASSIF_SUMMARY_SAMPLE <- "sample"
COLNAME_CLASSIF_SUMMARY_READ_TYPE <- "type"
COLNAME_CLASSIF_SUMMARY_N_READS <- "n_reads"