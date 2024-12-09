test_that("addRank() handles the taxon hierarchies appropriately", {
    logger::log_threshold(logger::OFF)
    # Define a few taxon hierarchies to be tested and the expected
    # outputs of addRank(), then carry out tests.

    # Domain.
    hierarchy_domain <- "d__Domain"
    expected_domain_out <- "D"
    actual_domain_out <- SPARKI:::addRank(hierarchy_domain)
    expect_identical(actual_domain_out, expected_domain_out)

    # Kingdom.
    hierarchy_kingdom <- "d__Domain|k__Kingdom"
    expected_kingdom_out <- "K"
    actual_kingdom_out <- SPARKI:::addRank(hierarchy_kingdom)
    expect_identical(actual_kingdom_out, expected_kingdom_out)

    # Phylum.
    hierarchy_phylum <- "d__Domain|k__Kingdom|p__Phylum"
    expected_phylum_out <- "P"
    actual_phylum_out <- SPARKI:::addRank(hierarchy_phylum)
    expect_identical(actual_phylum_out, expected_phylum_out)

    # Class.
    hierarchy_class <- "d__Domain|k__Kingdom|p__Phylum|c__Class"
    expected_class_out <- "C"
    actual_class_out <- SPARKI:::addRank(hierarchy_class)
    expect_identical(actual_class_out, expected_class_out)

    # Order.
    hierarchy_order <- "d__Domain|k__Kingdom|p__Phylum|c__Class|o__Order"
    expected_order_out <- "O"
    actual_order_out <- SPARKI:::addRank(hierarchy_order)
    expect_identical(actual_order_out, expected_order_out)

    # Family.
    hierarchy_family <- "d__Domain|k__Kingdom|p__Phylum|c__Class|o__Order|f__Family"
    expected_family_out <- "F"
    actual_family_out <- SPARKI:::addRank(hierarchy_family)
    expect_identical(actual_family_out, expected_family_out)

    # Genus.
    hierarchy_genus <- "d__Domain|k__Kingdom|p__Phylum|c__Class|o__Order|f__Family|g__Genus"
    expected_genus_out <- "G"
    actual_genus_out <- SPARKI:::addRank(hierarchy_genus)
    expect_identical(actual_genus_out, expected_genus_out)

    # Species.
    hierarchy_species <- "d__Domain|k__Kingdom|p__Phylum|c__Class|o__Order|f__Family|g__Genus|s__Species"
    expected_species_out <- "S"
    actual_species_out <- SPARKI:::addRank(hierarchy_species)
    expect_identical(actual_species_out, expected_species_out)
})
