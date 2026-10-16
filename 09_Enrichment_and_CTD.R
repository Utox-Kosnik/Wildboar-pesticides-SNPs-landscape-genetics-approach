################################################################################
# Main outputs
#  Figure 2, Figure S3, Figure S4, Figure S5
#  Tables S12-S17
################################################################################

rm(list = ls())
Sys.setenv(LANG = "en")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(ggtext)
  library(igraph)
})

################################################################################
# 1. FILES AND ANALYSIS SETTINGS
################################################################################

file <- list(
  lg_3km = "all_tox_wb_after_sens_duplicate_3km.RData",
  lg_4km = "all_tox_wb_after_sens_duplicate_4km.RData",
  lg_5km = "all_tox_wb_after_sens_duplicate_5km.RData",
  chemical_class = "ChemicalClass_pesticides.csv",
  hypergeom_function = "getHypergeomAssociation.R",
  pig_kegg = "250710_ssc_gene_kegg.csv",
  human_kegg = "250710_hsa_gene_kegg.csv",
  pig_annotation_dir = "Sus_scrofa",
  human_annotation_dir = "Homo_Sapiens",
  ctd_dir = "250630-CTD_Specific_Data",
  input_gene_map = "all_genes_represented_in_40700_inputdata.RData"
)

output_dir <- "revision_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

MIN_TERM_SIZE <- 5L
MAX_TERM_SIZE <- 1000L
MIN_OVERLAP <- 3L
Q_CUTOFF <- 0.05

NEURAL_OR_TRANSPORT_REGEX <-
  "nervous|neuro|neuron|axon|synap|dendrit|potentiation|transmembrane|channel|cation|ion|transport"
NEURAL_REGEX <- "nervous|neuro|neuron|axon|synap|dendrit|potentiation|spinocerebellar"

MF_DETOX_REGEX <-
  "detox|xenobiotic|oxidoreduct|transferase|hydrolase|lyase|catalytic|antioxid|glutathione|transmembrane|transport|channel"

required_files <- c(
  file$lg_4km, file$chemical_class, file$hypergeom_function,
  file$pig_kegg, file$human_kegg, file$input_gene_map,
  file.path(file$pig_annotation_dir, "OFFICIAL_GENE_SYMBOL2GOTERM_BP_ALL.txt"),
  file.path(file$pig_annotation_dir, "OFFICIAL_GENE_SYMBOL2GOTERM_CC_ALL.txt"),
  file.path(file$pig_annotation_dir, "OFFICIAL_GENE_SYMBOL2GOTERM_MF_ALL.txt"),
  file.path(file$pig_annotation_dir, "OFFICIAL_GENE_SYMBOL2REACTOME_PATHWAY.txt"),
  file.path(file$human_annotation_dir, "OFFICIAL_GENE_SYMBOL2GOTERM_BP_ALL.txt"),
  file.path(file$human_annotation_dir, "OFFICIAL_GENE_SYMBOL2GOTERM_CC_ALL.txt"),
  file.path(file$human_annotation_dir, "OFFICIAL_GENE_SYMBOL2GOTERM_MF_ALL.txt"),
  file.path(file$human_annotation_dir, "OFFICIAL_GENE_SYMBOL2REACTOME_PATHWAY.txt"),
  file.path(file$ctd_dir, "CTD_chem_gene_ixns.tsv.gz"),
  file.path(file$ctd_dir, "CTD_diseases.tsv.gz"),
  file.path(file$ctd_dir, "CTD_diseases_pathways.tsv.gz"),
  file.path(file$ctd_dir, "CTD_Phenotype-Disease_molecular_function_associations.tsv.gz"),
  file.path(file$ctd_dir, "CTD_Phenotype-Disease_biological_process_associations.tsv.gz"),
  file.path(file$ctd_dir, "CTD_Phenotype-Disease_cellular_component_associations.tsv.gz")
)

missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required)) {
  stop(
    "Required input files are missing:\n  ",
    paste(missing_required, collapse = "\n  ")
  )
}

################################################################################
# 2. HELPERS
################################################################################

`%notin%` <- function(x, table) !(x %in% table)

standardise_genes <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x == ""] <- NA_character_
  x
}

normalise_pesticide_class <- function(x) {
  z <- tolower(trimws(as.character(x)))
  fcase(
    grepl("insect", z), "Insecticides",
    grepl("fung", z), "Fungicides",
    grepl("herb", z), "Herbicides",
    default = "Other"
  )
}

find_column_name <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit)) hit[1] else NA_character_
}

process_go_file <- function(path) {
  ans <- fread(path, header = FALSE, sep = "\t", col.names = c("gene", "go_info")) %>%
    tidyr::separate(go_info, into = c("term_code", "term"), sep = "~", extra = "merge") %>%
    transmute(
      term = tolower(term),
      term_code,
      gene = standardise_genes(gene)
    ) %>%
    distinct()
  as.data.table(ans)
}

process_kegg_file <- function(path) {
  ans <- fread(path)[, .(
    term = pathway,
    term_code = pathway_code,
    gene = standardise_genes(gene)
  )]
  unique(ans)
}

filter_term_size <- function(x) {
  x[, term_n := .N, by = term]
  x[term_n >= MIN_TERM_SIZE & term_n <= MAX_TERM_SIZE][, term_n := NULL]
}

load_annotation_sets <- function(species = c("pig", "human")) {
  species <- match.arg(species)
  annotation_dir <- if (species == "pig") file$pig_annotation_dir else file$human_annotation_dir
  kegg_file <- if (species == "pig") file$pig_kegg else file$human_kegg

  sets <- list(
    BP = process_go_file(file.path(annotation_dir, "OFFICIAL_GENE_SYMBOL2GOTERM_BP_ALL.txt")),
    CC = process_go_file(file.path(annotation_dir, "OFFICIAL_GENE_SYMBOL2GOTERM_CC_ALL.txt")),
    MF = process_go_file(file.path(annotation_dir, "OFFICIAL_GENE_SYMBOL2GOTERM_MF_ALL.txt")),
    KEGG = process_kegg_file(kegg_file),
    Reactome = process_go_file(file.path(annotation_dir, "OFFICIAL_GENE_SYMBOL2REACTOME_PATHWAY.txt"))
  )

  background <- unique(unlist(lapply(sets, function(x) x$gene)))
  sets <- lapply(sets, filter_term_size)
  list(sets = sets, background = background)
}

source(file$hypergeom_function)
if (exists("dg_getHypergeomAssociation_v6")) {
  getHypergeomAssociation <- dg_getHypergeomAssociation_v6
}
if (!exists("getHypergeomAssociation")) {
  stop("getHypergeomAssociation() was not found after sourcing ", file$hypergeom_function)
}

run_enrichment <- function(annotation, input, group_col, background_n) {
  result <- rbindlist(
    lapply(names(annotation), function(annotation_class) {
      z <- getHypergeomAssociation(
        annotation[[annotation_class]], input,
        "gene", "gene_name", "fdr", background_n, TRUE
      )
      z$Class <- annotation_class
      as.data.table(z)
    }),
    fill = TRUE
  )

  group_gene_col <- paste0("genes_", group_col)
  if (!all(c("genes_term", group_gene_col) %in% names(result))) {
    stop("Unexpected output from getHypergeomAssociation(); missing genes_term or ", group_gene_col)
  }

  result[, genes := Map(intersect, genes_term, get(group_gene_col))]
  result <- result[num_gene_overlap >= MIN_OVERLAP]
  result[]
}

collapse_gene_list <- function(x) {
  vapply(x, function(z) paste(sort(unique(na.omit(z))), collapse = "|"), character(1))
}

add_manual_pesticide_classes <- function(x) {
  x[ChemicalName == "Forchlorfenuron", Pesticide_Type := "Other, Plants"]
  x[ChemicalName == "Bromadiolone", Pesticide_Type := "Other, Mammals"]
  x[ChemicalName == "Clodinafop", Pesticide_Type := "Herbicides"]
  x[ChemicalName == "Fosetyl", Pesticide_Type := "Fungicides"]
  x[ChemicalName == "Maleic hydrazide", Pesticide_Type := "Other, Plants"]
  x[ChemicalName == "Fluazifop-P", Pesticide_Type := "Herbicides"]
  x[, Pesticide_Class := normalise_pesticide_class(Pesticide_Type)]
  x
}

make_enrichment_panel <- function(dat, title) {
  z <- copy(dat)
  z[, term := factor(term, levels = rev(unique(term[order(p_value, -num_gene_overlap)])))]
  ggplot(z, aes(num_gene_overlap, term, fill = p_value)) +
    geom_col() +
    scale_fill_gradient(
      name = "p value", high = "#f54842", low = "#425af5",
      trans = "reverse",
      guide = guide_colorbar(direction = "horizontal", barwidth = 10, barheight = 2)) +
    scale_x_continuous(breaks = scales::breaks_width(5))+
    labs(x = NULL, y = NULL, title = title) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0),
      legend.position = "none",
      axis.text.y = element_text(size = 14)
    )
}

################################################################################
# 3. FINAL 4-KM LANDSCAPE-GENETIC DATA AND PESTICIDE CLASSES
################################################################################

lg <- as.data.table(readRDS(file$lg_4km))
stopifnot(all(c("snp", "tox_name", "gene_name") %in% names(lg)))
lg[, gene_name := standardise_genes(gene_name)]
lg[, ChemName_modified := tolower(trimws(tox_name))]

lg_type_col <- find_column_name(
  lg, c("Pesticide type", "Pesticide.type", "Pesticide_Type")
)
if (is.na(lg_type_col)) {
  lg[, Pesticide_Type_raw := NA_character_]
} else {
  lg[, Pesticide_Type_raw := as.character(get(lg_type_col))]
}

chemclass <- fread(file$chemical_class)
chemclass[, ChemName_modified := tolower(trimws(ChemName_modified))]
chemclass <- unique(chemclass)

lookup_type_col <- find_column_name(
  chemclass, c("Pesticide type", "Pesticide.type", "Pesticide_Type")
)
if (is.na(lookup_type_col)) {
  chemclass[, Pesticide_Type_lookup := NA_character_]
} else {
  chemclass[, Pesticide_Type_lookup := as.character(get(lookup_type_col))]
}

lg <- merge(lg, chemclass, by = "ChemName_modified", all.x = TRUE)
lg[, Pesticide_Type_raw := fcoalesce(Pesticide_Type_raw, Pesticide_Type_lookup)]
lg[, Pesticide_Class := normalise_pesticide_class(Pesticide_Type_raw)]

# Manual classifications also used in the earlier analysis.
lg[tox_name == "Bromadiolone", Pesticide_Class := "Other"]
lg[tox_name == "Dichlorprop-P", Pesticide_Class := "Herbicides"]
lg[tox_name == "Glufosinate", Pesticide_Class := "Herbicides"]
lg[tox_name == "Maleic hydrazide", Pesticide_Class := "Other"]
lg[tox_name == "Mecoprop-P", Pesticide_Class := "Herbicides"]
lg[tox_name == "Metaldehyde", Pesticide_Class := "Other"]
lg[tox_name == "Propamocarb", Pesticide_Class := "Fungicides"]


pesticide_class_counts <- unique(lg[, .(tox_name, Pesticide_Class)])[, .(
  N_pesticides = .N
), by = Pesticide_Class][order(Pesticide_Class)]
fwrite(pesticide_class_counts, file.path(output_dir, "Pesticide_class_counts.csv"))

gene_counts_per_pesticide <- unique(lg[!is.na(gene_name), .(tox_name, gene_name)])[, .(
  N_genes = uniqueN(gene_name)
), by = tox_name][order(-N_genes)]
fwrite(gene_counts_per_pesticide, file.path(output_dir, "Genes_per_pesticide.csv"))

################################################################################
# 4. SNP SHARING AMONG PESTICIDE CLASSES
################################################################################

snp_class <- unique(lg[, .(snp, Pesticide_Class)])
snp_class[, N_classes := uniqueN(Pesticide_Class), by = snp]

class_specific_snp_summary <- snp_class[, .(
  N_class_snps = uniqueN(snp),
  N_unique_to_class = uniqueN(snp[N_classes == 1L])
), by = Pesticide_Class]
class_specific_snp_summary[, Percent_of_class_snps_unique :=
  100 * N_unique_to_class / N_class_snps]
class_specific_snp_summary[, Percent_of_all_snps_unique_to_class :=
  100 * N_unique_to_class / uniqueN(snp_class$snp)]
setorder(class_specific_snp_summary, Pesticide_Class)
fwrite(class_specific_snp_summary, file.path(output_dir, "SNP_sharing_between_pesticide_classes.csv"))

################################################################################
# 5. SNP / SNP-PESTICIDE OVERLAP AMONG 3-, 4-, AND 5-KM RADII
################################################################################

pig <- load_annotation_sets("pig")

radius_files <- c(`3 km` = file$lg_3km, `4 km` = file$lg_4km, `5 km` = file$lg_5km)
if (all(file.exists(radius_files))) {
  radius_data <- rbindlist(lapply(names(radius_files), function(radius) {
    z <- as.data.table(readRDS(radius_files[[radius]]))
    stopifnot(all(c("snp", "tox_name") %in% names(z)))
    unique(z[, .(radius, snp, tox_name)])
  }))

  make_membership <- function(dat, key_cols, label) {
    z <- unique(dat[, c(key_cols, "radius"), with = FALSE])
    z[, present := 1L]
    radius_formula <- as.formula(
      paste(paste(key_cols, collapse = " + "), "~ radius")
    )
    z <- dcast(z, radius_formula, value.var = "present", fill = 0L)
    radius_cols <- intersect(c("3 km", "4 km", "5 km"), names(z))
    z[, combination := apply(.SD > 0, 1, function(v) {
      paste(radius_cols[v], collapse = " + ")
    }), .SDcols = radius_cols]
    z[, .(N = .N), by = combination][, Unit := label][]
  }

  radius_overlap <- rbindlist(list(
    make_membership(radius_data, "snp", "Unique SNPs"),
    make_membership(radius_data, c("snp", "tox_name"), "SNP-pesticide associations")
  ))
  radius_overlap[, combination := factor(
    combination,
    levels = c("3 km", "4 km", "5 km", "3 km + 4 km", "3 km + 5 km",
               "4 km + 5 km", "3 km + 4 km + 5 km")
  )]
  setorder(radius_overlap, Unit, combination)
  fwrite(radius_overlap, file.path(output_dir, "Radius_overlap_counts.csv"))

  p_radius <- ggplot(radius_overlap, aes(combination, N, fill = Unit)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = N), vjust = -0.25, size = 3) +
    facet_wrap(~Unit, scales = "free_y", ncol = 1) +
    labs(x = "Radius membership", y = "Count") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(file.path(output_dir, "FigureS3A_radius_overlap.pdf"), p_radius,
         width = 8, height = 7)
  
  ##########################################################################
  # ENRICHMENT OF SNP-PESTICIDE ASSOCIATIONS FOUND AT ALL THREE RADII
  ##########################################################################
  
  # Retain exact SNP-pesticide combinations that occur at 3, 4 and 5 km.
  common_radius_pairs <- radius_data[,.(N_radii = uniqueN(radius)),
    by = .(snp,tox_name)][
    N_radii == length(radius_files),
    .(snp,tox_name)]
  
  fwrite(
    common_radius_pairs,file.path(output_dir,"Common_3km_4km_5km_SNP_pesticide_associations.csv"))
  
  # Retrieve the corresponding gene annotations from the final 4-km
  # landscape-genetic dataset.
  common_radius_gene_map <- unique(
    merge(
      common_radius_pairs,
      lg[!is.na(gene_name), .(snp,tox_name,gene_name)],by = c("snp","tox_name")))
  
  fwrite(
    common_radius_gene_map,file.path(output_dir,"Common_3km_4km_5km_SNP_pesticide_gene_map.csv"))
  
  # Create one unique list of genes for the enrichment analysis.
  common_radius_gene_input <- unique(common_radius_gene_map[,.(input = "Common_3km_4km_5km",gene_name)])
  
  if (nrow(common_radius_gene_input) > 0) {
    
    # Run the same enrichment analysis used for the overall Figure 2.
    common_radius_enrichment <- run_enrichment(
      annotation = pig$sets,
      input = common_radius_gene_input,
      group_col = "input",
      background_n = length(pig$background)
    )[
      ,
      .(Class,term,qvalue,mvalue,nvalue,kvalue,p_value,num_gene_overlap,genes)]
    
    # Export the complete enrichment results.
    common_radius_enrichment_out <- copy(
      common_radius_enrichment)
    
    common_radius_enrichment_out[,genes := collapse_gene_list(genes)]
    
    fwrite(common_radius_enrichment_out,file.path(output_dir,"SITableS13_common_3km_4km_5km_enrichment.csv"))
    
    # Count enriched terms within each annotation class.
    common_radius_enrichment_counts <-
      common_radius_enrichment[,.(N_enriched_terms = uniqueN(term) ),by = Class]
    
    fwrite(common_radius_enrichment_counts,file.path(output_dir,"Common_3km_4km_5km_enrichment_counts.csv" ))
    
    # Summary numbers for reporting.
    common_radius_summary <- data.table(
      N_common_associations = nrow(common_radius_pairs),
      N_unique_snps = uniqueN(common_radius_pairs$snp),
      N_pesticides = uniqueN(common_radius_pairs$tox_name),
      N_intragenic_associations = nrow(common_radius_gene_map),
      N_unique_genes = uniqueN(common_radius_gene_map$gene_name),
      N_enriched_terms = uniqueN(
        common_radius_enrichment[,.(Class,term)]))
    
    fwrite(common_radius_summary,file.path(output_dir,"Common_3km_4km_5km_enrichment_summary.csv"))
    
    ########################################################################
    # ENRICHMENT FIGURE
    ########################################################################
    
    if (nrow(common_radius_enrichment)) {
      common_radius_figure_data <- copy(common_radius_enrichment)
      common_radius_figure_data[, Plot_group := fcase(
        Class == "BP", "BP",
        Class == "CC", "CC",
        Class == "MF", "MF",
        Class %in% c("KEGG", "Reactome"), "KEGG + Reactome"
      )]
      
      common_radius_figure_data[, Plot_group := factor(
        Plot_group, levels = c("BP", "CC", "MF", "KEGG + Reactome")
      )]
      
      # Order the most significant terms at the top of each panel
      setorder(common_radius_figure_data, Plot_group, p_value, -num_gene_overlap)
      common_radius_figure_data[, term_panel := paste(term, Class, sep = "|||")]
      common_radius_figure_data[, term_panel := factor(
        term_panel, levels = rev(unique(term_panel))
      )]
      
      # Round the maximum x-axis value upwards to the nearest 5
      common_x_max <- max(5,5 * ceiling(max(common_radius_figure_data$num_gene_overlap, na.rm = TRUE) / 5))
      panel_labels <- c(BP = "B P",CC = "C C",MF = "M F", `KEGG + Reactome` = "KEGG /\nREACTOME")
      common_radius_figure <- ggplot(
        common_radius_figure_data,
        aes(x = num_gene_overlap, y = term_panel, fill = p_value)
      ) +
        geom_col() +
        scale_x_continuous(
          limits = c(0, common_x_max),
          breaks = seq(0, common_x_max, by = 5),
          expand = c(0, 0)
        ) +
        scale_y_discrete(labels = function(x) sub("\\|\\|\\|.*$", "", x)) +
        scale_fill_gradient(
          name = "p value",
          high = "#f54842",
          low = "#425af5",
          trans = "reverse",
          guide = guide_colorbar(
            direction = "vertical",
            barwidth = 0.7,
            barheight = 8
          )
        ) +
        facet_grid(
          Plot_group ~ .,
          scales = "free_y",
          space = "free_y",
          switch = "y",
          labeller = labeller(
            Plot_group = as_labeller(panel_labels)
          )
        ) +
        labs(x = NULL, y = NULL) +
        theme_bw() +
        theme(
          legend.position = "right",
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 10),
          panel.spacing.y = grid::unit(0.08, "lines"),
          strip.background = element_rect(fill = "grey85", colour = "black"),
          strip.text.y.left = element_text(
            angle = 0,
            face = "bold",
            size = 12,
            margin = margin(6, 8, 6, 8)
          )
        )
      
      ggsave(
        file.path(output_dir, "FigureS3B_common_3km_4km_5km_enrichment.pdf"),
        common_radius_figure,
        width = 10,
        height = 5
      )
      
    } else {
      
      warning(
        paste(
          "No significant enrichment terms were found for SNP-pesticide",
          "associations shared by all three radii; no enrichment plot was made."
        )
      )
    }
    
  } else {
    
    warning(
      paste(
        "No intragenic SNP-pesticide associations were shared by all three",
        "radii; no common-radius enrichment analysis was run."
      )
    )
  }
} else {
  warning(
    "Radius-overlap outputs were skipped. Missing: ",
    paste(radius_files[!file.exists(radius_files)], collapse = ", ")
  )
}

################################################################################
# 6. PIG ANNOTATIONS AND ENRICHMENT
################################################################################

all_genes <- unique(lg[!is.na(gene_name), .(input = "AllBoar", gene_name)])
overall_enrichment_raw <- run_enrichment(
  pig$sets, all_genes, "input", length(pig$background)
)
overall_enrichment <- overall_enrichment_raw[, .(
  Class, term, qvalue, mvalue, nvalue, kvalue, p_value,
  num_gene_overlap, genes
)]

per_pesticide_input <- unique(lg[!is.na(gene_name), .(tox_name, gene_name)])
per_pesticide_raw <- run_enrichment(
  pig$sets, per_pesticide_input, "tox_name", length(pig$background)
)
per_pesticide_raw[, ChemName_modified := tolower(tox_name)]
chemical_info <- unique(lg[, .(
  ChemName_modified, CAS, Pesticide_Type_raw, Chem.group, MOA
)])
per_pesticide <- merge(
  per_pesticide_raw, chemical_info,
  by = "ChemName_modified", all.x = TRUE
)
per_pesticide <- per_pesticide[, .(
  Class, term, ChemicalName = tox_name, ChemName_modified,
  CASRN = CAS, Pesticide_Type = Pesticide_Type_raw, Chem.group, MOA,
  qvalue, mvalue, nvalue, kvalue, p_value, num_gene_overlap, genes
)]
per_pesticide <- add_manual_pesticide_classes(per_pesticide)

class_input <- unique(lg[!is.na(gene_name), .(Pesticide_Class, gene_name)])
class_enrichment_raw <- run_enrichment(
  pig$sets, class_input, "Pesticide_Class", length(pig$background)
)
class_enrichment <- class_enrichment_raw[, .(
  Class, term, Pesticide_Class, qvalue, mvalue, nvalue, kvalue,
  p_value, num_gene_overlap, genes
)]

class_gene_counts <- class_input[, .(N_genes = uniqueN(gene_name)), by = Pesticide_Class]
fwrite(class_gene_counts, file.path(output_dir, "Genes_per_pesticide_class.csv"))

################################################################################
# 7. OVERALL ENRICHMENT COUNTS, KEYWORD COUNTS, TOP TERMS, FIGURE 2, TABLE S5
################################################################################

overall_counts <- overall_enrichment[, .(N_terms = uniqueN(term)), by = Class]
overall_counts <- rbind(
  overall_counts,
  data.table(Class = "Total", N_terms = uniqueN(overall_enrichment[, .(Class, term)]))
)
fwrite(overall_counts, file.path(output_dir, "Overall_enrichment_counts.csv"))

overall_term_flags <- unique(overall_enrichment[, .(Class, term)])
overall_term_flags[, neural_or_transport :=
  grepl(NEURAL_OR_TRANSPORT_REGEX, term, ignore.case = TRUE)]
overall_term_flags[, neural := grepl(NEURAL_REGEX, term, ignore.case = TRUE)]
overall_term_flags[, mf_detox :=
  Class == "MF" & grepl(MF_DETOX_REGEX, term, ignore.case = TRUE)]

keyword_summary <- rbindlist(list(
  overall_term_flags[Class %in% c("BP", "CC","MF","KEGG","Reactome"), .(
    Category = "Neuronal or membrane-transport term",
    N_matched = sum(neural_or_transport),
    N_total = .N
  ), by = Class]
))

keyword_summary[, Percent := 100 * N_matched / N_total]
fwrite(keyword_summary, file.path(output_dir, "Overall_keyword_summary.csv"))
fwrite(
  overall_term_flags[neural_or_transport | mf_detox | neural],
  file.path(output_dir, "Overall_keyword_matched_terms.csv")
)

# 'Top' is defined reproducibly as the largest gene overlap, with p value used
# as the tie-breaker. This is preferable to the ambiguous phrase 'in terms of
# significance and gene overlap'.
top_terms <- copy(overall_enrichment)
top_terms[, Plot_group := fcase(Class == "BP", "BP",Class == "CC", "CC",Class == "MF", "MF",default = "KEGG + Reactome")]
setorder(top_terms, Plot_group, -num_gene_overlap, p_value)
top_terms <- top_terms[, head(.SD, 3), by = Plot_group]
fwrite(top_terms, file.path(output_dir, "Top_three_enriched_terms.csv"))

figure2_data <- copy(overall_enrichment)
figure2_data[, Plot_group := fcase(Class == "BP", "BP",Class == "CC", "CC",Class == "MF", "MF",default = "KEGG + Reactome")]

p_bp <- make_enrichment_panel(figure2_data[Plot_group == "BP"], "Biological Processes")
p_cc <- make_enrichment_panel(figure2_data[Plot_group == "CC"], "Cellular Components")
p_mf <- make_enrichment_panel(figure2_data[Plot_group == "MF"], "Molecular Functions")
p_path <- make_enrichment_panel(figure2_data[Plot_group == "KEGG + Reactome"],
                                "KEGG / Reactome Pathways")

legend_plot <- ggplot(figure2_data, aes(num_gene_overlap, term, fill = p_value)) +
  geom_col() +
  scale_fill_gradient(
    name = "p value", high = "#f54842", low = "#425af5", trans = "reverse",
    guide = guide_colorbar(direction = "horizontal", barwidth = 15, barheight = 1)
  ) +
  theme_void() + theme(legend.position = "bottom",
                       legend.text = element_text(size = 14),
                       legend.title = element_text(size = 14))

figure2 <- plot_grid(
  plot_grid(p_bp, p_cc, p_mf, p_path, ncol = 2,
            labels = c("A", "B", "C", "D"), label_size = 20),
  get_legend(legend_plot),
  ncol = 1, rel_heights = c(1, 0.06)
)
ggsave(file.path(output_dir, "Figure2_overall_enrichment.pdf"), figure2,
       width = 18, height = 12)

si7 <- copy(overall_enrichment)
si7[, genes := collapse_gene_list(genes)]
fwrite(si7[, .(Class, term, qvalue, p_value, num_gene_overlap, genes)],
       file.path(output_dir, "SITableS12.csv"))

################################################################################
# 8. PER-PESTICIDE ENRICHMENT SUMMARY AND TABLE S6
################################################################################

per_pesticide_term_counts <- unique(per_pesticide[, .(
  ChemicalName, Pesticide_Class, Class, term
)])[, .(N_enriched_terms = .N), by = .(ChemicalName, Pesticide_Class)][
  order(-N_enriched_terms)
]
fwrite(per_pesticide_term_counts,
       file.path(output_dir, "Enriched_terms_per_pesticide.csv"))

si9 <- copy(per_pesticide)
si9[, genes := collapse_gene_list(genes)]
fwrite(si9[, .(
  ChemicalName, CASRN, Pesticide_Class, Chem.group, MOA, Class, term,
  qvalue, p_value, num_gene_overlap, genes
)], file.path(output_dir, "SITableS14.csv"))

################################################################################
# 9. CTD DISEASE CLASSIFICATION FOR ENRICHED TERMS
################################################################################

ctd_disease <- fread(file.path(file$ctd_dir, "CTD_diseases.tsv.gz"))[
  , .(DiseaseName = `#`, DiseaseID = V2, SlimMappings = V9)
]
ctd_disease <- as.data.table(separate_rows(ctd_disease, SlimMappings, sep = "[|]"))

read_ctd_term_disease <- function(path, class) {
  x <- fread(path)
  if (class %in% c("KEGG", "Reactome")) {
    x <- x[, .(DiseaseName = `#`, DiseaseID = V2, term = V3, termID = V4)]
    x[, Class := fifelse(grepl("KEGG", termID), "KEGG", "Reactome")]
  } else {
    x <- x[, .(term = `#`, termID = V2, DiseaseName = V3, DiseaseID = V4)]
    x[, Class := class]
  }
  unique(x)
}

ctd_term_disease <- rbindlist(list(
  read_ctd_term_disease(file.path(file$ctd_dir, "CTD_diseases_pathways.tsv.gz"), "KEGG"),
  read_ctd_term_disease(file.path(file$ctd_dir,
    "CTD_Phenotype-Disease_molecular_function_associations.tsv.gz"), "MF"),
  read_ctd_term_disease(file.path(file$ctd_dir,
    "CTD_Phenotype-Disease_biological_process_associations.tsv.gz"), "BP"),
  read_ctd_term_disease(file.path(file$ctd_dir,
    "CTD_Phenotype-Disease_cellular_component_associations.tsv.gz"), "CC")
), fill = TRUE)

# Optional CTD node downloads used in the original script. If present, include
# them; if absent, report which enriched terms remain unclassified.
node_files <- data.table(
  filename = c(
    "CTD_GO_0098916_diseases_20250722122559.tsv",
    "CTD_GO_0045229_diseases_20250722122657.tsv",
    "CTD_GO_0009894_diseases_20250722122828.tsv",
    "CTD_GO_0031344_diseases_20250722122954.tsv",
    "CTD_GO_0051588_diseases_20250722123031.tsv",
    "CTD_GO_0033043_diseases_20250722123104.tsv",
    "CTD_GO_0099177_diseases_20250722123243.tsv",
    "CTD_GO_0062023_diseases_20250722123332.tsv",
    "CTD_GO_0099501_diseases_20250722123415.tsv",
    "CTD_GO_0030312_diseases_20250722123435.tsv",
    "CTD_GO_0016835_diseases_20250722123638.tsv",
    "CTD_GO_0140098_diseases_20250722123701.tsv",
    "CTD_GO_0015318_diseases_20250722123730.tsv",
    "CTD_GO_0015399_diseases_20250722123757.tsv"
  ),
  term = c(
    "anterograde trans-synaptic signaling",
    "external encapsulating structure organization",
    "regulation of catabolic process",
    "regulation of cell projection organization",
    "regulation of neurotransmitter transport",
    "regulation of organelle organization",
    "regulation of trans-synaptic signaling",
    "collagen-containing extracellular matrix",
    "plasma membrane signaling receptor complex",
    "external encapsulating structure",
    "carbon-oxygen lyase activity",
    "catalytic activity, acting on RNA",
    "inorganic molecular entity transmembrane transporter activity",
    "primary active transmembrane transporter activity"
  ),
  termID = c(
    "GO:0098916", "GO:0045229", "GO:0009894", "GO:0031344",
    "GO:0051588", "GO:0033043", "GO:0099177", "GO:0062023",
    "GO:0099501", "GO:0030312", "GO:0016835", "GO:0140098",
    "GO:0015318", "GO:0015399"
  ),
  Class = c(rep("BP", 7), rep("CC", 3), rep("MF", 4))
)

available_nodes <- node_files[file.exists(file.path(file$ctd_dir, filename))]
if (nrow(available_nodes)) {
  node_terms <- rbindlist(lapply(seq_len(nrow(available_nodes)), function(i) {
    meta <- available_nodes[i]
    x <- fread(file.path(file$ctd_dir, meta$filename))[, 1:4]
    setnames(x, c("term", "termID", "DiseaseName", "DiseaseID"))
    x[, `:=`(term = meta$term, termID = meta$termID, Class = meta$Class)]
    x
  }), fill = TRUE)
  ctd_term_disease <- rbind(ctd_term_disease, node_terms, fill = TRUE)
}

ctd_term_class <- merge(
  unique(ctd_term_disease), ctd_disease,
  by = c("DiseaseName", "DiseaseID"), allow.cartesian = TRUE
)[, .(N_diseases = uniqueN(DiseaseID)), by = .(Class, term, termID, SlimMappings)]
ctd_term_class[, max_n := max(N_diseases), by = .(Class, term, termID)]
ctd_term_class <- unique(ctd_term_class[N_diseases == max_n,
  .(Class, term, termID, SlimMappings, N_diseases)])
ctd_term_class[, term_key := tolower(trimws(term))]

classify_enrichment_terms <- function(x) {
  
  x <- copy(x)
  x[, term_key := tolower(trimws(term))]
  
  z <- merge(
    x,ctd_term_class[, .(Class,term_key,termID,SlimMappings,N_diseases)],
    by = c("Class", "term_key"),all.x = TRUE,allow.cartesian = TRUE)
  
  # Identify terms with multiple equally supported MEDIC categories
  ambiguous <- z[,uniqueN(na.omit(SlimMappings)),by = .(Class, term)][V1 > 1L]
  
  # Ambiguous terms -> Other
  z[ambiguous,on = .(Class, term),SlimMappings := "Other"]
  
  # Unclassified terms -> Other
  z[is.na(SlimMappings),SlimMappings := "Other"]
  
  # Collapse disease classes into the three comparison groups
  z[, CompareClass := fifelse(SlimMappings %in% c("Nervous system disease","Cancer"),
    SlimMappings,"Other")]
  
  dedup_cols <- names(z)[!vapply(z, is.list, logical(1))]
  z[!duplicated(z[, ..dedup_cols])]
}

per_pesticide_classified <- classify_enrichment_terms(per_pesticide)
class_enrichment_classified <- classify_enrichment_terms(class_enrichment)


chemical_disease_summary <- unique(per_pesticide_classified[, .(
  ChemicalName, Pesticide_Class, Class, term, CompareClass
)])[, .(
  N_pesticides = uniqueN(ChemicalName),
  N_enriched_terms = uniqueN(paste(Class, term))
), by = .(Pesticide_Class, CompareClass)][order(Pesticide_Class, CompareClass)]
fwrite(chemical_disease_summary,
       file.path(output_dir, "Table_Pesticide_class_disease_summary_per_chemical.csv")) #ENRICHMENT BASED ON Each pesticide separately

class_neural_summary <- unique(class_enrichment_classified[, .(
  Pesticide_Class, Class, term, CompareClass
)])[, .(
  N_total_terms = .N,
  N_neural_terms = sum(CompareClass == "Nervous system disease")
), by = Pesticide_Class]
class_neural_summary[, Percent_neural := 100 * N_neural_terms / N_total_terms]
setorder(class_neural_summary, Pesticide_Class)
fwrite(class_neural_summary,
       file.path(output_dir, "Table_Pesticide_class_level_neural_enrichment.csv")) #ENRICHMENT BASED ON Genes pooled within each pesticide class

class_plot_data <- unique(class_enrichment_classified[, .(
  Pesticide_Class, Class, term, num_gene_overlap, CompareClass
)])
class_plot_data[, Plot_group := fcase(
  Class == "BP", "Biological Processes",
  Class == "CC", "Cellular Components",
  Class == "MF", "Molecular Functions",
  default = "KEGG / Reactome Pathways"
)]
class_plot_data[, term_order := sum(num_gene_overlap), by = .(Plot_group, term)]
class_plot_data[, term := factor(term, levels = unique(term[order(term_order)]))]

class_cols <- c(
  Herbicides = "green3", Fungicides = "salmon",
  Insecticides = "skyblue2", Other = "grey70"
)
p_class <- ggplot(class_plot_data,
                  aes(num_gene_overlap, term, fill = Pesticide_Class)) +
  geom_col() +
  facet_wrap(~Plot_group, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = class_cols, drop = FALSE) +
  labs(x = "Gene count", y = NULL, fill = NULL) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 14),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )
ggsave(file.path(output_dir, "FigureS4_class_enrichment.pdf"), p_class,
       width = 18, height = 12)

class_table <- copy(class_enrichment_classified)
class_table[, genes := collapse_gene_list(genes)]
fwrite(class_table[, .(
  Pesticide_Class, Class, term, qvalue, p_value, num_gene_overlap,
  CompareClass, genes
)], file.path(output_dir, "SITableS15_class_enrichment.csv"))

################################################################################
# 10. CTD CHEMICAL-GENE DATA AND SUS SCROFA GENE-LEVEL VALIDATION
################################################################################

ctd <- fread(file.path(file$ctd_dir, "CTD_chem_gene_ixns.tsv.gz"))
ctd <- ctd[, .(
  N_publications = uniqueN(V11)
), by = .(
  ChemicalName = `#`, CASRN = V3, GeneSymbol = standardise_genes(V4),
  GeneID = V5, Species = V7
)]
ctd[, ChemName_modified := tolower(ChemicalName)]

# Match by standardised name first and CASRN second, then restrict to chemicals
# actually represented in the final LG dataset.
lg_chemicals <- unique(lg[, .(ChemName_modified, tox_name)])
ctd_by_name <- ctd[ChemName_modified %in% lg_chemicals$ChemName_modified]
ctd_unmatched <- ctd[ChemName_modified %notin% lg_chemicals$ChemName_modified]
cas_lookup <- unique(chemclass[
  ChemName_modified %in% lg_chemicals$ChemName_modified & !is.na(CAS),
  .(CASRN = CAS, ChemName_modified)
])
ctd_by_cas <- merge(
  ctd_unmatched, cas_lookup, by = "CASRN",
  suffixes = c("_CTD", ""), allow.cartesian = TRUE
)
ctd_by_cas[, ChemName_modified_CTD := NULL]
ctd_lg <- unique(rbind(ctd_by_name, ctd_by_cas, fill = TRUE))
ctd_lg <- ctd_lg[ChemName_modified %in% lg_chemicals$ChemName_modified]

ctd_pig <- ctd_lg[Species == "Sus scrofa"]
pig_ctd_pesticides <- sort(unique(ctd_pig$ChemName_modified))

lg_pig_genes <- unique(lg[
  ChemName_modified %in% pig_ctd_pesticides & !is.na(gene_name), gene_name
])
ctd_pig_genes <- unique(na.omit(ctd_pig$GeneSymbol))

input_gene_map <- as.data.table(readRDS(file$input_gene_map))
input_genes <- unique(na.omit(standardise_genes(input_gene_map$gene_name)))

pig_validation_summary <- data.table(
  Metric = c(
    "CTD pesticides represented in Sus scrofa",
    "Unique LG genes for those pesticides",
    "Unique Sus scrofa CTD genes",
    "Sus scrofa CTD genes present in initial SNP input",
    "Overlapping unique genes",
    "Exact pesticide-gene overlaps"
  ),
  Value = c(
    length(pig_ctd_pesticides), length(lg_pig_genes), length(ctd_pig_genes),
    length(intersect(ctd_pig_genes, input_genes)),
    length(intersect(lg_pig_genes, ctd_pig_genes)), NA_integer_))

lg_pairs <- unique(lg[
  ChemName_modified %in% pig_ctd_pesticides & !is.na(gene_name),
  .(ChemName_modified, gene_name)])

ctd_pig_pairs <- unique(ctd_pig[, .(ChemName_modified, gene_name = GeneSymbol)])
pig_direct_overlap <- merge(lg_pairs, ctd_pig_pairs,
                            by = c("ChemName_modified", "gene_name"))
pig_validation_summary[Metric == "Exact pesticide-gene overlaps", Value := nrow(pig_direct_overlap)]
fwrite(pig_validation_summary,
       file.path(output_dir, "Sus_scrofa_CTD_validation_summary.csv"))

################################################################################
# 11. HUMAN CTD ENRICHMENT, OVERLAPPING TERMS, TABLES S16-S17, FIGURE S5
################################################################################

human <- load_annotation_sets("human")
ctd_human <- ctd_lg[Species == "Homo sapiens"]
human_ctd_input <- unique(ctd_human[, .(
  ChemName_modified, gene_name = GeneSymbol)])

human_ctd_enrichment_raw <- run_enrichment(
  human$sets, human_ctd_input, "ChemName_modified", length(human$background))

human_ctd_enrichment <- human_ctd_enrichment_raw[, .(
  ChemName_modified, Class, term, p_value, qvalue,
  num_gene_overlap, genes)]

# Work in long gene format for exact gene-within-pathway comparisons.
ctd_enrich_long <- copy(human_ctd_enrichment)
ctd_enrich_long[, genes := collapse_gene_list(genes)]
ctd_enrich_long <- as.data.table(separate_rows(ctd_enrich_long, genes, sep = "[|]"))

lg_enrich_long <- copy(per_pesticide_classified)
lg_enrich_long[, genes := collapse_gene_list(genes)]
lg_enrich_long <- as.data.table(separate_rows(lg_enrich_long, genes, sep = "[|]"))

term_overlap <- merge(ctd_enrich_long, lg_enrich_long,by = c("ChemName_modified", "Class", "term"), 
                      suffixes = c("_CTD", "_LG"), allow.cartesian = TRUE)

overlap_terms <- unique(term_overlap[, .(
  ChemName_modified, Class, term,
  ChemicalName, CASRN, Pesticide_Class, CompareClass)])

exact_gene_pathway <- unique(term_overlap[
  genes_CTD == genes_LG,
  .(ChemName_modified, Class, term, gene_name = genes_CTD)])

ctd_term_genes <- unique(term_overlap[, .(
  ChemName = ChemName_modified, Class, term,
  Gene_Source = "CTD_genes", genes = genes_CTD)])

lg_term_genes <- unique(term_overlap[, .(
  ChemName = ChemName_modified, Class, term,
  Gene_Source = "LG_genes", genes = genes_LG)])

si8 <- rbind(ctd_term_genes, lg_term_genes)
si8 <- si8[, .(genes = paste(sort(unique(genes)), collapse = "|")),
           by = .(ChemName, Class, term, Gene_Source)]
fwrite(si8, file.path(output_dir, "TableS17_CTD_LG_overlap.csv"))

# Direct pesticide-gene overlap with CTD across all available species. Aggregate
# species so the number of associations is the number of unique pesticide-gene
# pairs rather than duplicated species rows.
lg_all_pairs <- unique(lg[!is.na(gene_name), .(ChemName_modified, ChemicalName = tox_name, gene_name)])
ctd_all_pairs <- unique(ctd_lg[, .(ChemName_modified, gene_name = GeneSymbol, Species)])
direct_overlap_all_species <- merge(lg_all_pairs, ctd_all_pairs,by = c("ChemName_modified", "gene_name"))[, .(
  CTD_Species = paste(sort(unique(Species)), collapse = "|")), by = .(ChemicalName, gene_name)]
fwrite(direct_overlap_all_species, file.path(output_dir, "SITableS16.csv"))

# Genes contributing to overlapping terms: distinguish genes represented in
# the human CTD enrichment but linked to another pesticide from genes absent
# from the corresponding CTD starting data.
lg_overlap_genes <- unique(lg_term_genes$genes)
ctd_enriched_genes_any_pesticide <- unique(ctd_enrich_long$genes)
ctd_starting_genes <- unique(ctd_human$GeneSymbol)

# Require the supporting CTD enrichment to involve a different pesticide.
lg_gene_pesticide <- unique(lg_term_genes[, .(LG_pesticide = ChemName, genes)])
ctd_gene_pesticide <- unique(ctd_enrich_long[, .(CTD_pesticide = ChemName_modified, genes)])
different_pesticide_support <- merge(lg_gene_pesticide, ctd_gene_pesticide, 
                                     by = "genes", allow.cartesian = TRUE)[LG_pesticide != CTD_pesticide]

lg_gene_support_summary <- data.table(
  Metric = c(
    "LG genes contributing to overlapping terms",
    "Represented in human CTD enrichment for another pesticide",
    "Absent from corresponding human CTD starting data"),
  Value = c(length(lg_overlap_genes),uniqueN(different_pesticide_support$genes),
              length(setdiff(lg_overlap_genes, ctd_starting_genes))))
fwrite(lg_gene_support_summary,
       file.path(output_dir, "LG_gene_support_within_overlapping_terms.csv"))

human_validation_summary <- data.table(
  Metric = c(
    "LG pesticides with human CTD chemical-gene interactions",
    "Unique direct pesticide-gene overlaps across CTD species",
    "Pesticides sharing at least one enriched term with human CTD",
    "Unique enriched terms shared for the same pesticide",
    "Exact pesticide-gene-pathway overlaps"
  ),
  Value = c(
    uniqueN(ctd_human$ChemName_modified),
    nrow(direct_overlap_all_species),
    uniqueN(overlap_terms$ChemName_modified),
    uniqueN(overlap_terms[, .(Class, term)]),
    nrow(exact_gene_pathway)
  )
)
fwrite(human_validation_summary,
       file.path(output_dir, "Human_CTD_validation_summary.csv"))
fwrite(exact_gene_pathway,
       file.path(output_dir, "Exact_CTD_LG_gene_pathway_overlaps.csv"))

# Network for Figure S5. Term squares are coloured by pesticide. Terms shared
# by several pesticides receive a separate "Multiple pesticides" colour.
network_edges <- rbind(
  unique(lg_term_genes[,.(from = term,to = genes,Gene_Source,ChemName)]),
  unique(ctd_term_genes[,.(from = term,to = genes,Gene_Source,ChemName)]))

network_edges[,ChemName := tolower(trimws(ChemName))]

# Determine whether each term belongs to one pesticide or several pesticides.
term_chemical <- unique(network_edges[,.(term = from,ChemName)])
term_chemical[ ,N_chemicals := uniqueN(ChemName),by = term]
term_chemical[,Display_chemical := fifelse(N_chemicals > 1L, "multiple",ChemName)]
term_chemical <- unique(term_chemical[,.(term,Display_chemical)])

# Pesticide colours for GO/pathway squares.
chemical_colors <- c("cypermethrin" = "goldenrod1","dimethomorph" = "green3","folpet"       = "pink",
  "malathion"    = "skyblue3","pyrimethanil" = "cyan3","thiram"       = "salmon","multiple"     = "mediumpurple3")

term_chemical[,term_color := unname(chemical_colors[Display_chemical])]

# Fallback colour for unexpected pesticide names.
term_chemical[is.na(term_color),term_color := "grey70"]

if (nrow(network_edges) > 0) {
  
  ##############################################################################
  # Build network
  ##############################################################################
  
  net <- graph_from_data_frame(network_edges,directed = FALSE)
  
  ##############################################################################
  # Identify term and gene nodes
  ##############################################################################
  
  term_names <- unique(network_edges$from)
  lg_gene_names <- unique(network_edges[Gene_Source == "LG_genes",to])
  ctd_gene_names <- unique(network_edges[Gene_Source == "CTD_genes",to])
  
  # Genes represented somewhere in both datasets.
  both_source_genes <- intersect(lg_gene_names,ctd_gene_names)
  
  # Genes matching the same pesticide and GO/pathway term.
  exact_genes <- unique(exact_gene_pathway$gene_name)
  term_nodes <- V(net)$name %in% term_names
  lg_nodes <- V(net)$name %in% lg_gene_names
  ctd_nodes <- V(net)$name %in% ctd_gene_names
  both_source_nodes <- V(net)$name %in% both_source_genes
  exact_nodes <- V(net)$name %in% exact_genes
  
  ##############################################################################
  # Assign node shapes
  ##############################################################################
  
  # GO/pathway terms are squares; genes are circles.
  V(net)$shape <- ifelse(term_nodes,"square","circle")
  
  ##############################################################################
  # Assign node colours
  ##############################################################################
  
  V(net)$color <- "grey70"
  
  # Colour GO/pathway squares by pesticide.
  term_match <- match(V(net)$name[term_nodes],term_chemical$term)
  V(net)$color[term_nodes] <-term_chemical$term_color[term_match]
  
  # Landscape-genetic genes.
  V(net)$color[lg_nodes] <- "orange3"
  
  # Human CTD genes.
  V(net)$color[ctd_nodes] <- adjustcolor("seagreen2",alpha.f = 0.45)
  
  # Genes present in both datasets, but not an exact match.
  V(net)$color[both_source_nodes & !exact_nodes] <- "lightskyblue"
  
  # Exact pesticide–gene–pathway matches.
  V(net)$color[exact_nodes] <- "yellow"
  V(net)$frame.color <- "black"
  V(net)$label <- ""
  
  ##############################################################################
  # Create and save Figure S5
  ##############################################################################
  
  set.seed(123)
  pdf(file.path(output_dir,"FigureS5_CTD_network.pdf"),width = 12,height = 10)
  par(mar = c(1, 1, 1, 1),xpd = NA)
  
  plot(net,layout = layout_with_fr(net),vertex.size = 5,vertex.label = NA,
    edge.color = adjustcolor("grey60",alpha.f = 0.5 ),edge.width = 0.5)
  
  ##############################################################################
  # Pesticide legend
  ##############################################################################
  
  legend("topleft",
    legend = c("Cypermethrin","Dimethomorph","Folpet","Malathion","Pyrimethanil","Thiram","Multiple pesticides"),
    pch = 22,pt.bg = chemical_colors[
      c("cypermethrin","dimethomorph","folpet","malathion","pyrimethanil","thiram","multiple")],
    col = "black",pt.cex = 1.5,cex = 0.8,bty = "n",title = "GO/pathway terms")
  
  ##############################################################################
  # Gene legend
  ##############################################################################
  
  legend("bottomleft",legend = c("Landscape-genetic gene","Human CTD gene",
                                 "Gene in both datasets, but not exact","Exact pesticide-gene-pathway match"),
    pch = 21,pt.bg = c("orange3",adjustcolor("seagreen2",alpha.f = 0.45 ),"lightskyblue","yellow"),
    col = "black",pt.cex = 1.5,cex = 0.8,bty = "n",title = "Genes")
  dev.off()
}
