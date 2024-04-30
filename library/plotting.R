plot_classified_vs_unclassified_summary <- function(merged_reports, return_plot = FALSE, outdir) {

    class_unclass_df <- get_class_unclass_numbers(merged_reports)

    plot <- ggplot2::ggplot(
            class_unclass_df, 
            ggplot2::aes(x = type, y = value, group = sample)
        ) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 12, vjust = 0.5),
            axis.text.y = ggplot2::element_text(size = 12),
            axis.title.x = ggplot2::element_text(size = 14),
            axis.title.y = ggplot2::element_text(size = 14, angle = 90),
            legend.text = ggplot2::element_text(size = 12),
            legend.title = ggplot2::element_text(size = 14),
            legend.position = "left",
            legend.justification = "top",
            plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = "bold"),
            plot.margin = ggplot2::unit(c(0, 0, 0, 0), "pt")
        ) +
        ggplot2::xlab("Type\n") +
        ggplot2::ylab("# Reads\n")

}

plot_classified_vs_unclassified_per_sample <- function(merged_reports, return_plot = FALSE, outdir) {

    class_unclass_df <- get_class_unclass_numbers(merged_reports)

    plot <- ggplot2::ggplot(
            class_unclass_df, 
            ggplot2::aes(x = value, y = sample, fill = type)
        ) +
        ggplot2::geom_bar(position = "fill", stat = "identity") +
        ggplot2::theme_void() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 12, vjust = 0.5),
            axis.text.y = ggplot2::element_text(size = 12),
            axis.title.x = ggplot2::element_text(size = 14),
            axis.title.y = ggplot2::element_text(size = 14, angle = 90),
            legend.text = ggplot2::element_text(size = 12),
            legend.title = ggplot2::element_text(size = 14),
            legend.position = "left",
            legend.justification = "top",
            plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = "bold"),
            plot.margin = ggplot2::unit(c(0, 0, 0, 0), "pt")
        ) +
        ggplot2::xlab("Proportion of reads\n") +
        ggplot2::ylab("Sample\n") +
        ggplot2::scale_fill_brewer(
            palette = "Accent",
            name = "Type",
            labels = c("Classified", "Unclassified")
        )

    pdf_height = determine_pdf_height(n_elements = length(unique(merged_mpa$sample)), factor = 0.25)

    if (return_plot == FALSE) {

        if (!(missing(outdir))) {
            ggpubr::ggexport(
                plot, 
                filename = paste0(outdir, "/proportion_classified_unclassified.pdf"),
                width = 10,
                height = pdf_height
            )
        } else {
            stop("Plot is to be saved to a file but no output directory has been provided.")
        }
    } else if (return_plot == TRUE && missing(outdir)) {
        return(plot)
    }
}

plot_reads_all_domains <- function(merged_mpa, return_plot = FALSE, outdir) {

    domains <- c("Archaea", "Bacteria", "Eukaryota", "Viruses")
    names(domains) <- c("d__Archaea", "d__Bacteria", "d__Eukaryota", "d__Viruses")

    merged_mpa <- merged_mpa[grep("d__Eukaryota$|d__Viruses$|d__Archaea$|d__Bacteria$", merged_mpa$classification), ]

    plot <- ggplot2::ggplot(
            merged_mpa, ggplot2::aes(x = classification, y = log10(number_fragments_clade), fill = classification)
        ) +
        ggplot2::geom_violin(scale = "width") +
        ggplot2::geom_jitter(size = 0.5) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(size = 15),
            axis.title = ggplot2::element_text(size = 16),
            legend.position = "none",
            strip.text.x = ggplot2::element_text(size = 15),
            axis.ticks.x = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(size = 16.5, face = "bold", hjust = 0.5)
        ) +
        ggplot2::xlab(NULL) +
        ggplot2::ylab(expression("log"[10]~"(number of reads classified)")) +
        ggplot2::scale_fill_manual(values = c("gold1", "royalblue", "snow2", "indianred2")) +
        ggplot2::facet_wrap(
            ~classification, 
            scales = "free",
            labeller = ggplot2::labeller(classification = domains)
        )

    if (return_plot == FALSE) {

        if (!(missing(outdir))) {
            ggpubr::ggexport(plot, filename = paste0(outdir, "/reads_all_domains.pdf"))
        } else {
            stop("Plot is to be saved to a file but no output directory has been provided.")
        }
    } else if (return_plot == TRUE && missing(outdir)) {
        return(plot)
    }
}

plot_reads_non_eukaryote_domains <- function(n_reads_in_samples_df, merged_mpa, return_plot = FALSE, outdir) {

    merged_mpa <- merged_mpa[grep("d__Viruses$|d__Archaea$|d__Bacteria$", merged_mpa$classification), ]

    subplot1 <- ggplot2::ggplot(
            merged_mpa, 
            ggplot2::aes(x = number_fragments_clade, y = sample, fill = classification)
        ) +
        ggplot2::geom_bar(position = "fill", stat = "identity") +
        ggplot2::theme_void() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 12),
            axis.text.y = ggplot2::element_text(size = 12),
            axis.title.x = ggplot2::element_text(size = 14),
            axis.title.y = ggplot2::element_text(size = 14, angle = 90),
            legend.text = ggplot2::element_text(size = 14),
            legend.title = ggplot2::element_text(size = 16),
            legend.position = "left",
            legend.justification = "top",
            plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = "bold"),
            plot.margin = ggplot2::unit(c(0, 0, 0, 0), "pt")
        ) +
        ggplot2::xlab("Proportion of reads") +
        ggplot2::ylab("Sample\n") +
        ggplot2::scale_fill_brewer(
            palette = "Pastel1",
            name = "Domain",
            labels = c("Archaea", "Bacteria", "Viruses")
        ) 

    subplot2 <- ggplot2::ggplot(
            n_reads_in_samples_df[n_reads_in_samples_df$domains_considered == "subset", ],
            ggplot2::aes(x = log10(number_reads_clade), y = sample)
        ) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_void() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 12),
            axis.title.x = ggplot2::element_text(size = 14),
            plot.margin = ggplot2::unit(c(0, 0, 0, 0), "pt")
        ) +
        ggplot2::xlab(expression("log"[10]~"(total non-eukaryote reads)")) 

    plot <- cowplot::plot_grid(
        plotlist = list(subplot1, subplot2),
        ncol = 2,
        align = "hv",
        axis = "tb",
        rel_widths = c(1.15,1)
    ) 

    pdf_height = determine_pdf_height(n_elements = length(unique(merged_mpa$sample)), factor = 0.25)

    if (return_plot == FALSE) {

        if (!(missing(outdir))) {
            ggpubr::ggexport(
                plot, 
                filename = paste0(outdir, "/reads_non_eukaryote_domains.pdf"),
                width = 10,
                height = pdf_height
            )
        } else {
            stop("Plot is to be saved to a file but no output directory has been provided.")
        }
    } else if (return_plot == TRUE && missing(outdir)) {
        return(plot)
    }
}
