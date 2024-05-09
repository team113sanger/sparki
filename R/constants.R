######################
## HELPER CONSTANTS ##
#######################################################################################################

#### Reference database ####
COLNAME_PCT_FRAG_CLADE_REF_DB <- "pct_fragments_clade"
COLNAME_MINIMISERS_CLADE_REF_DB <- "n_minimisers_clade"
COLNAME_MINIMISERS_TAXON_REF_DB <- "n_minimisers_taxon"
COLNAME_RANK_REF_DB <- "rank"
COLNAME_NCBI_ID_REF_DB <- "ncbi_id"
COLNAME_TAXON_REF_DB <- "taxon"

#### Standard report ####
COLNAME_SAMPLE_STD <- "sample"
COLNAME_PCT_FRAG_CLADE_STD <- "pct_fragments_clade"
COLNAME_N_FRAG_CLADE_STD <- "n_fragments_clade"
COLNAME_N_FRAG_TAXON_STD <- "n_fragments_taxon"
COLNAME_MINIMISERS_STD <- "n_minimisers"
COLNAME_UNIQ_MINIMISERS_STD <- "n_distinct_minimisers"
COLNAME_RANK_STD <- "rank"
COLNAME_NCBI_ID_STD <- "ncbi_id"
COLNAME_TAXON_STD <- "taxon"
COLNAME_DB_MINIMISERS_CLADE_STD <- "db_number_minimisers_clade"
COLNAME_DB_MINIMISERS_TAXON_STD <- "db_number_minimisers_taxon"

#### MPA-style report #### 
COLNAME_SAMPLE_MPA <- "sample"
COLNAME_N_FRAG_CLADE_MPA <- "n_fragments_clade"
COLNAME_TAXON_MPA <- "taxon_hierarchy"
COLNAME_RANK_MPA <- "rank"
COLNAME_NCBI_ID_MPA <- "ncbi_id"
COLNAME_DB_MINIMISERS_CLADE_MPA <- "db_number_minimisers_clade"
COLNAME_DB_MINIMISERS_TAXON_MPA <- "db_number_minimisers_taxon"

#### Auxiliary dataframes ####

# Number of domain reads per sample.
COLNAME_SAMPLE_DOMAIN_READS <- "sample"
COLNAME_SCOPE_DOMAIN_READS <- "domains_considered"
COLNAME_N_READS_DOMAIN_READS <- "n_clade_reads"

# Number of classified/unclassified reads per sample.
COLNAME_SAMPLE_CLASSIF_SUMMARY <- "sample"
COLNAME_READ_TYPE_CLASSIF_SUMMARY <- "type"
COLNAME_N_READS_CLASSIF_SUMMARY <- "n_reads"