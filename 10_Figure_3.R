################################################################################
# FIGURE 3: MANHATTAN PLOTS FOR SIX PESTICIDES
#
# Run 09_enrichment_and_CTD.R first. This script uses its SITableS17.csv
################################################################################

rm(list = ls())
Sys.setenv(LANG = "en")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(vegan)
  library(qvalue)
})

################################################################################
# 1. SETTINGS AND INPUTS
################################################################################

chemical_list <- c(
  "Cypermethrin", "Dimethomorph", "Folpet",
  "Malathion", "Pyrimethanil", "Thiram"
)

# These must match the final RDA used in the manuscript.
conditioning_variables <- c(
  "agricultural_percentage", "human_population", 
  "latitude", "longitude"
)

file <- list(
  genotype = "wb_matrix_no_na_282.RData",
  metadata = "wb_final_metadata_4km.RData",
  snp_manifest = "forR_PorcineSNP60.csv",
  map_names = "names_SNPs.txt",
  final_lg = "all_tox_wb_after_sens_duplicate_4km.RData",
  ctd_overlap = file.path("TableS17_CTD_LG_overlap.csv")
)

output_dir <- "revision_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

missing_files <- unlist(file)[!file.exists(unlist(file))]
if (length(missing_files)) {
  stop("Required files are missing:\n  ", paste(missing_files, collapse = "\n  "))
}

wb_matrix <- readRDS(file$genotype)
metadata <- as.data.frame(readRDS(file$metadata))
lg <- as.data.table(readRDS(file$final_lg))
# Read SI Table S12
ctd_overlap <- fread(
  file.path(output_dir, "TableS17_CTD_LG_overlap.csv")
)

# Expand the genes column so each gene has its own row
ctd_overlap_long <- as.data.table(
  tidyr::separate_rows(
    ctd_overlap,
    genes,
    sep = "\\|"
  )
)

# Now standardise the names
ctd_overlap_long[, ChemName := tolower(trimws(as.character(ChemName)))]
ctd_overlap_long[, genes := toupper(trimws(as.character(genes)))]

ctd_overlap_long <- ctd_overlap_long[
  !is.na(genes) & genes != ""
]

manifest <- fread(file$snp_manifest)
map <- fread(file$map_names, header = FALSE)
data.table::setnames(map, c("chrom_old", "snp", "x", "pos_old"))
map[, locus := as.character(seq_len(.N))]

required_metadata <- c("id", chemical_list, conditioning_variables)
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata)) {
  stop("Metadata columns are missing: ", paste(missing_metadata, collapse = ", "))
}
if (!identical(rownames(wb_matrix), as.character(metadata$id))) {
  stop("Genotype-matrix row order does not exactly match metadata$id")
}

################################################################################
# 2. RUN SIX RDAs AND RETAIN ALL TESTED SNPs
################################################################################

run_one_rda <- function(chemical) {
  message("Running RDA for ", chemical)
  toxic_column <- metadata[[chemical]]
  rda_formula <- as.formula(
    paste0(
      "wb_matrix ~ toxic_column + Condition(",
      paste(conditioning_variables, collapse = " + "), ")"
    ),
    env = environment()
  )
  set.seed(123)
  model <- vegan::rda(
    rda_formula, data = metadata, scale = TRUE, na.action = na.omit
  )

  loadings <- as.data.frame(vegan::scores(model, choices = 1, display = "species"))
  loadings$locus <- rownames(loadings)
  loading_col <- setdiff(names(loadings), "locus")[1]
  loadings$scaled_loading <- as.numeric(scale(loadings[[loading_col]]))
  loadings$pvalue <- 2 * (1 - pnorm(abs(loadings$scaled_loading)))
  loadings$qvalue <- qvalue::qvalue(loadings$pvalue)$qvalues

  mapped <- merge(
    map[, .(snp, locus)],
    manifest[, .(snp = Name, chromosome = Chr, position = MapInfo)],
    by = "snp"
  )
  mapped <- merge(
    mapped,
    as.data.table(loadings)[, .(locus = as.character(locus), qvalue)],
    by = "locus"
  )
  mapped[, tox_name := chemical]
  mapped[, .(tox_name, snp, chromosome, position, qvalue)]
}

all_tested <- rbindlist(lapply(chemical_list, run_one_rda), fill = TRUE)
all_tested[, `:=`(
  chromosome = as.numeric(as.character(chromosome)),
  position = as.numeric(as.character(position)),
  qvalue = as.numeric(as.character(qvalue))
)]
all_tested <- all_tested[
  chromosome %between% c(1, 18) & !is.na(position) & !is.na(qvalue)
]

################################################################################
# 3. APPLY THE FINAL SENSITIVITY FILTER TO INITIALLY SIGNIFICANT SNPs
################################################################################

# All non-significant points remain as the grey background. Among initially
# significant SNPs, retain only those present in the final LG results.
final_significant <- unique(lg[
  tox_name %in% chemical_list, .(tox_name, snp)
])
initial_significant <- unique(all_tested[qvalue < 0.05, .(tox_name, snp)])

initial_significant[, initially_significant := TRUE]
final_significant[, retained_final := TRUE]
manhattan <- merge(all_tested, initial_significant,
                   by = c("tox_name", "snp"), all.x = TRUE)
manhattan <- merge(manhattan, final_significant,
                   by = c("tox_name", "snp"), all.x = TRUE)
manhattan <- manhattan[is.na(initially_significant) | retained_final == TRUE]
manhattan[, c("initially_significant", "retained_final") := NULL]

################################################################################
# 4. CTD HIGHLIGHT CATEGORIES
################################################################################

stopifnot(all(c("ChemName", "Class", "term", "Gene_Source", "genes") %in%
                names(ctd_overlap)))
ctd_long <- as.data.table(separate_rows(ctd_overlap, genes, sep = "[|]"))
ctd_long[, `:=`(
  ChemName = tolower(trimws(ChemName)),
  genes = toupper(trimws(genes))
)]
ctd_long <- ctd_long[!is.na(genes) & genes != ""]

ctd_overlap_long[, ChemName := tolower(trimws(as.character(ChemName)))]
ctd_overlap_long[, genes := toupper(trimws(as.character(genes)))]

# Category 1: LG gene contributing to a term shared with human CTD for the same
# pesticide. Category 2: exact same pesticide + term + gene in both datasets.
pathway_gene_lookup <- unique(ctd_long[
  Gene_Source == "LG_genes", .(ChemName, gene_name = genes)
])
lg_pathway_genes <- unique(ctd_long[
  Gene_Source == "LG_genes", .(ChemName, Class, term, gene_name = genes)
])
human_pathway_genes <- unique(ctd_long[
  Gene_Source == "CTD_genes", .(ChemName, Class, term, gene_name = genes)
])
exact_gene_pathway <- merge(
  lg_pathway_genes, human_pathway_genes,
  by = c("ChemName", "Class", "term", "gene_name")
)
exact_gene_lookup <- unique(exact_gene_pathway[, .(ChemName, gene_name)])

lg_gene_annotations <- unique(lg[
  tox_name %in% chemical_list & !is.na(gene_name),
  .(
    tox_name,
    ChemName = tolower(tox_name),
    snp,
    gene_name = toupper(gene_name)
  )
])
lg_gene_annotations[, Highlight_SNP := "X"]
lg_gene_annotations[pathway_gene_lookup,
  on = .(ChemName, gene_name),
  Highlight_SNP := "Homo sapiens gene in overlapping pathway"
]
lg_gene_annotations[exact_gene_lookup,
  on = .(ChemName, gene_name),
  Highlight_SNP := "Homo sapiens gene & pathway"
]
lg_gene_annotations[, priority := fcase(
  Highlight_SNP == "Homo sapiens gene & pathway", 2L,
  Highlight_SNP == "Homo sapiens gene in overlapping pathway", 1L,
  default = 0L
)]

lg_gene_annotations[, ChemName := tolower(trimws(as.character(ChemName)))]
lg_gene_annotations[, gene_name := toupper(trimws(as.character(gene_name)))]

snp_highlights <- lg_gene_annotations[, {
  p <- max(priority)
  list(
    Highlight_SNP = fcase(
      p == 2L, "Homo sapiens gene & pathway",
      p == 1L, "Homo sapiens gene in overlapping pathway",
      default = "X"
    ),
    gene_label = if (p > 0L) {
      paste(sort(unique(gene_name[priority == p])), collapse = ", ")
    } else ""
  )
}, by = .(tox_name, snp)]

manhattan <- merge(manhattan, snp_highlights,
                   by = c("tox_name", "snp"), all.x = TRUE)
manhattan[is.na(Highlight_SNP), Highlight_SNP := "X"]
manhattan[is.na(gene_label), gene_label := ""]

################################################################################
# 5. TRUE GENOMIC X-AXIS AND LABEL LOCATIONS
################################################################################

chromosome_info <- manhattan[, .(
  chromosome_length = max(position, na.rm = TRUE)
), by = chromosome][order(chromosome)]
chromosome_info[, offset := c(0, head(cumsum(chromosome_length), -1))]
chromosome_info[, chromosome_center := offset + chromosome_length / 2]

manhattan <- merge(
  manhattan, chromosome_info[, .(chromosome, offset)],
  by = "chromosome", all.x = TRUE
)
manhattan[, genomic_position := position + offset]

label_data <- manhattan[
  Highlight_SNP != "X" & gene_label != "" & !is.na(qvalue)
][order(qvalue), .SD[1], by = .(tox_name, gene_label)]

manhattan[, tox_name := factor(tox_name, levels = chemical_list)]
label_data[, tox_name := factor(tox_name, levels = chemical_list)]

################################################################################
# 6. FIGURE 3 AND SUPPORTING TABLES
################################################################################

p <- ggplot() +
  geom_point(
    data = manhattan[Highlight_SNP == "X"],
    aes(genomic_position, -log10(qvalue)),
    colour = "grey70", size = 0.55, alpha = 0.65
  ) +
  geom_point(
    data = manhattan[Highlight_SNP != "X"],
    aes(genomic_position, -log10(qvalue), colour = Highlight_SNP),
    size = 1.5
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.8) +
  geom_label_repel(data = label_data,aes(x = genomic_position,y = -log10(qvalue),
      label = gene_label,colour = Highlight_SNP),
    size = 2.3,show.legend = FALSE,
    
    # Repelling
    direction = "both",
    force = 2,
    force_pull = 0.2,
    box.padding = 0.5,
    point.padding = 0.3,
    
    # Allow labels to move further away
    max.overlaps = Inf,
    max.time = 5,
    max.iter = 100000,
    
    # Always draw connector lines
    min.segment.length = 0,
    
    # Give labels more vertical room
    ylim = c(-Inf, Inf)
  )+
  scale_colour_manual(
    values = c(
      "Homo sapiens gene in overlapping pathway" = "darkblue",
      "Homo sapiens gene & pathway" = "red"
    ),
    breaks = c(
      "Homo sapiens gene in overlapping pathway",
      "Homo sapiens gene & pathway"
    ),
    name = NULL
  )  +
  scale_x_continuous(
    breaks = chromosome_info$chromosome_center,
    labels = chromosome_info$chromosome,
    expand = c(0.005, 0.005)
  ) +
  facet_wrap(~tox_name, ncol = 2, scales="free_y") +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.30))
  )+
  labs(x = "Chromosome", y = expression(-log[10]("q-value"))) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(size = 11),
    axis.text.x = element_text(size = 7),
    legend.position = "bottom"
  )

ggsave(file.path(output_dir, "Figure3_six_chemical_manhattan.jpeg"), p,
       width = 9, height = 11, units = "in", dpi = 300)
ggsave(file.path(output_dir, "Figure3_six_chemical_manhattan.pdf"), p,
       width = 9, height = 11, units = "in")

highlighted_genes <- unique(manhattan[
  Highlight_SNP != "X", .(tox_name, gene_label, Highlight_SNP)
])
figure_summary <- manhattan[, .(
  N_SNPs_plotted = uniqueN(snp),
  N_significant = uniqueN(snp[qvalue < 0.05]),
  N_highlighted = uniqueN(snp[Highlight_SNP != "X"])
), by = tox_name]

fwrite(manhattan, file.path(output_dir, "Figure3_manhattan_data.csv"))
fwrite(highlighted_genes, file.path(output_dir, "Figure3_highlighted_CTD_genes.csv"))
fwrite(figure_summary, file.path(output_dir, "Figure3_summary.csv"))
fwrite(exact_gene_pathway,
       file.path(output_dir, "Figure3_exact_gene_pathway_matches.csv"))

cat("\nFigure 3 complete.\n")
print(figure_summary)
cat("\nExact CTD/LG pesticide-gene-pathway matches:\n")
print(exact_gene_pathway[order(ChemName, gene_name)])



