# Coral virome auxiliary metabolic gene (AMG) heatmap.
# Usage: Rscript scripts/run_amg_heatmap.R [input_directory] [output_directory]
# Defaults: module inputs/ and module root. Output subdirectories: results/, figures/.
# Inputs: final 24,083-genome DNA TPM table, sample metadata and AMG annotations.
#
# Analysis: join positive TPM records to AMG annotations; sum within sample and
# gene/category; average observed positive records within ecological groups; fill
# absent group/gene cells with zero; display log1p(mean TPM). Gene and group
# ordering, colors, ecological grouping and Project labels follow the supplied
# final analysis. Project labels annotate groups without changing aggregation.
# No clustering is applied. A TSV matrix and separate row/column annotations
# provide the exact values and ordering used in the final heatmap.
args <- commandArgs(trailingOnly = TRUE)
script_flag <- grep('^--file=', commandArgs(), value = TRUE)
module_root <- if (length(script_flag)) dirname(dirname(normalizePath(sub('^--file=', '', script_flag[[1]]), winslash = '/', mustWork = TRUE))) else getwd()
input_dir <- normalizePath(if (length(args) >= 1) args[[1]] else file.path(module_root, 'inputs'), winslash = '/', mustWork = TRUE)
output_root <- if (length(args) >= 2) args[[2]] else module_root
results_dir <- file.path(output_root, 'results')
figures_dir <- file.path(output_root, 'figures')
for (folder in c(results_dir, figures_dir)) dir.create(folder, recursive = TRUE, showWarnings = FALSE)
results_dir <- normalizePath(results_dir, winslash = '/', mustWork = TRUE)
figures_dir <- normalizePath(figures_dir, winslash = '/', mustWork = TRUE)
setwd(results_dir)

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(gridExtra)

# Reproducible seed for the supplied analysis.
set.seed(123)

# ==============================================================================
# 1. Read the final input tables.
# ==============================================================================
print(">>> Reading input tables...")
tpm_df <- read_tsv(file.path(input_dir, "dna_tpm.tsv"), show_col_types = FALSE)
meta_df <- read_tsv(file.path(input_dir, "sample_metadata.tsv"), show_col_types = FALSE)
amg_info <- read_tsv(file.path(input_dir, "amg_annotations.tsv"), show_col_types = FALSE)
stopifnot(nrow(tpm_df) == 24083L, !anyDuplicated(tpm_df$Genome))

# ==============================================================================
# 2. Harmonize sample and status metadata.
# ==============================================================================
print(">>> Harmonizing metadata...")

if(!"Sample" %in% colnames(meta_df) && "SampleID" %in% colnames(meta_df)) {
  meta_df <- meta_df %>% rename(Sample = SampleID)
}

allowed_tissues <- c("Coral fragment", "Coral mucus", "Coral skeleton", "Coral tissue")

clean_status_string <- function(x) {
  x <- gsub("−", "-", x)
  return(x)
}

meta_clean <- meta_df %>%
  filter(Sample %in% colnames(tpm_df)) %>%
  mutate(
    Status_Clean = clean_status_string(Status),
    Sample_Type_Raw = if("Sample_Type" %in% colnames(.)) Sample_Type else "Unknown",
    Sample_Type = ifelse(Sample_Type_Raw %in% allowed_tissues, Sample_Type_Raw, "Other"),
    Sample_Type = factor(Sample_Type, levels = c(allowed_tissues, "Other")),
    
    # Preserve missing ocean metadata as Unknown.
    Ocean_Basin = if("Ocean_Basin" %in% colnames(.)) Ocean_Basin else "Unknown",
    Ocean_Basin = ifelse(is.na(Ocean_Basin) | Ocean_Basin == "", "Unknown", Ocean_Basin),
    
    Block_Group = case_when(
      grepl("DSS.*enrichment", Status_Clean, ignore.case = TRUE) ~ "Other",
      grepl("Treated|Recovery", Status_Clean, ignore.case = TRUE) ~ "Other",
      grepl("Unknown", Status_Clean, ignore.case = TRUE) ~ "Other",
      grepl("Healthy|Normal", Status_Clean, ignore.case = TRUE) ~ "Healthy",
      Status_Clean == "Bleached" ~ "Bleached",
      grepl("SCTLD|White Plague|Black Band|Dark Spot|Pathogen|Bleaching lesion", Status_Clean, ignore.case = TRUE) ~ "Disease",
      grepl("Disease", Status_Clean, ignore.case = TRUE) ~ "Disease",
      grepl("Stress", Status_Clean, ignore.case = TRUE) ~ "Stressed",
      TRUE ~ "Other"
    ),
    
    Detail_Status = case_when(
      grepl("Healthy|Normal", Status_Clean, ignore.case = TRUE) ~ "Healthy/Normal",
      Status_Clean == "Bleached" ~ "Bleached",
      grepl("Bleaching lesion", Status_Clean, ignore.case = TRUE) ~ "Disease (Bleaching lesion)",
      grepl("Pathogen", Status_Clean, ignore.case = TRUE) ~ "Disease (Pathogen-induced disease)",
      grepl("Dark Spot", Status_Clean, ignore.case = TRUE) ~ "Disease(Dark Spot Syndrome)",
      grepl("SCTLD", Status_Clean, ignore.case = TRUE) ~ "Disease(SCTLD)",
      grepl("Heat stress", Status_Clean, ignore.case = TRUE) ~ "Stressed(Heat stress)",
      grepl("Acidification", Status_Clean, ignore.case = TRUE) ~ "Stressed(Acidification stress)",
      grepl("Combined", Status_Clean, ignore.case = TRUE) ~ "Stressed(Combined)",
      grepl("Treated|Recovery", Status_Clean, ignore.case = TRUE) ~ "Treated/Recovery",
      TRUE ~ "Other"
    ),
    Species = Coral_Host_Species
  )

block_levels <- c("Healthy", "Bleached", "Disease", "Stressed", "Other")
meta_clean$Block_Group <- factor(meta_clean$Block_Group, levels = block_levels)
all_species <- sort(unique(meta_clean$Species))
all_oceans <- sort(unique(meta_clean$Ocean_Basin))

# ==============================================================================
# 3. Assign positive genome TPM to AMG annotation records and sum by sample.
# ==============================================================================
print(">>> Calculating annotated AMG abundance...")

amg_clean <- amg_info %>%
  select(Viral_genome, Level1_Name, Level2_Name, KO_Short) %>%
  mutate(Gene_Name = gsub(",.*", "", KO_Short)) %>% 
  filter(!is.na(Gene_Name) & Gene_Name != "-" & Gene_Name != "") %>%
  mutate(
    Level1_Name = ifelse(is.na(Level1_Name) | Level1_Name == "", "Unclassified", Level1_Name),
    Level2_Name = ifelse(is.na(Level2_Name) | Level2_Name == "", "Unclassified", Level2_Name)
  )

# A genome can contribute more than one AMG annotation; retain annotation rows
# as supplied, so repeated gene annotations contribute separately to Total_TPM.
tpm_long <- tpm_df %>%
  pivot_longer(-Genome, names_to = "Sample", values_to = "TPM") %>%
  filter(TPM > 0)

amg_abundance <- tpm_long %>%
  inner_join(amg_clean, by = c("Genome" = "Viral_genome")) %>%
  group_by(Sample, Level1_Name, Level2_Name, Gene_Name) %>%
  summarise(Total_TPM = sum(TPM), .groups = "drop")

# ==============================================================================
# 4. Order genes by functional category and total annotated abundance.
# ==============================================================================
print(">>> Ordering genes...")

gene_stats <- amg_abundance %>%
  group_by(Level1_Name, Level2_Name, Gene_Name) %>%
  summarise(Gene_Total_TPM = sum(Total_TPM), .groups = "drop")

l2_stats <- gene_stats %>%
  group_by(Level1_Name, Level2_Name) %>%
  summarise(L2_Sum_TPM = sum(Gene_Total_TPM), .groups = "drop")

l1_stats <- gene_stats %>%
  group_by(Level1_Name) %>%
  summarise(L1_Gene_Count = n_distinct(Gene_Name), .groups = "drop") %>%
  mutate(L1_Sort_Key = ifelse(grepl("Other", Level1_Name, ignore.case = TRUE), -1, L1_Gene_Count)) %>%
  arrange(desc(L1_Sort_Key))

master_gene_list <- gene_stats %>%
  left_join(l2_stats, by = c("Level1_Name", "Level2_Name")) %>%
  left_join(l1_stats, by = "Level1_Name") %>%
  mutate(Level1_F = factor(Level1_Name, levels = l1_stats$Level1_Name)) %>%
  arrange(Level1_F, desc(L2_Sum_TPM), desc(Gene_Total_TPM))

fixed_genes <- master_gene_list$Gene_Name
ordered_l1_levels <- l1_stats$Level1_Name[!is.na(l1_stats$Level1_Name)]
ordered_l2_levels <- unique(master_gene_list$Level2_Name)

ann_genes <- data.frame(
  Level1 = factor(master_gene_list$Level1_Name, levels = ordered_l1_levels),
  Level2 = factor(master_gene_list$Level2_Name, levels = ordered_l2_levels)
)
rownames(ann_genes) <- fixed_genes

# ==============================================================================
# 5. Construct the grouped abundance matrix and order ecological groups.
# ==============================================================================
print(">>> Constructing the grouped matrix...")

plot_data_full <- amg_abundance %>%
  inner_join(meta_clean, by = "Sample")

# Group by coral species, basin, broad/detail status and sample type.
# Only observed positive-TPM sample/gene records contribute to this mean; absent
# group/gene combinations are filled with zero after aggregation.
heatmap_data_full <- plot_data_full %>%
  group_by(Species, Ocean_Basin, Block_Group, Detail_Status, Sample_Type, Gene_Name) %>% 
  summarise(Mean_TPM = mean(Total_TPM, na.rm=TRUE), .groups="drop") %>%
  mutate(Row_ID = paste(Species, Ocean_Basin, Block_Group, Detail_Status, Sample_Type, sep = "||"))

# Convert group/gene records to a wide matrix.
mat_raw <- heatmap_data_full %>%
  select(Row_ID, Gene_Name, Mean_TPM) %>%
  pivot_wider(names_from = Gene_Name, values_from = Mean_TPM, values_fill = 0) %>%
  column_to_rownames("Row_ID")

mat_fixed <- matrix(0, nrow = nrow(mat_raw), ncol = length(fixed_genes))
colnames(mat_fixed) <- fixed_genes
rownames(mat_fixed) <- rownames(mat_raw)
common_genes <- intersect(colnames(mat_raw), fixed_genes)
mat_fixed[, common_genes] <- as.matrix(mat_raw[, common_genes])

# Natural-log display transformation: ln(mean TPM + 1).
mat_log <- log1p(mat_fixed) 

# Order groups by species representation, tissue, ocean and status.
row_info_raw <- data.frame(Row_ID = rownames(mat_fixed)) %>%
  separate(Row_ID, into = c("Species", "Ocean", "Block", "Detail", "Tissue"), sep = "\\|\\|", remove = FALSE)

species_row_counts <- row_info_raw %>% count(Species) %>% rename(Total_Conditions = n)

row_info <- row_info_raw %>%
  left_join(species_row_counts, by = "Species") %>%
  mutate(
    Block = factor(Block, levels = block_levels),
    Tissue = factor(Tissue, levels = c(allowed_tissues, "Other")),
    Ocean = factor(Ocean, levels = sort(unique(Ocean)))
  ) %>%
  arrange(desc(Total_Conditions), Species, Tissue, Ocean, Block, Detail)

# Apply the final group ordering.
mat_ordered_samples <- mat_log[row_info$Row_ID, , drop=FALSE]

# ==============================================================================
# 6. Transpose: genes become rows and ecological groups become columns.
# ==============================================================================
print(">>> Transposing the matrix...")
mat_final <- t(mat_ordered_samples) 

# ==============================================================================
# 7. Define functional and ecological annotations and colors.
# ==============================================================================

ann_col_final <- data.frame(
  Species = row_info$Species,
  Tissue = row_info$Tissue,
  Ocean = row_info$Ocean,   
  Block = row_info$Block, 
  Detail = row_info$Detail
)
rownames(ann_col_final) <- rownames(mat_ordered_samples)

ann_row_final <- ann_genes

# Annotation palettes used in the final heatmap.

# 1. Ocean colors used in the final analysis
ocean_colors <- c(
  "Mediterranean Sea" = "#E41A1C", # Red
  "Pacific Ocean"     = "#377EB8", # Blue
  "Red Sea"           = "#4DAF4A", # Green
  "Arabian Gulf"      = "#984EA3", # Purple
  "Atlantic Ocean"    = "#FF7F00", # Orange
  "Caribbean Sea"     = "#FFFF33", # Yellow
  "Indian Ocean"      = "#A65628", # Brown
  "Unknown"           = "#999999"  # Grey for others
)

# 2. Tissue Colors
tissue_colors <- c(
  "Coral tissue"   = "#D95F02", 
  "Coral mucus"    = "#1B9E77", 
  "Coral skeleton" = "#8D40CA", 
  "Coral fragment" = "#E72989", 
  "Other"          = "#999999"
)

# 3. Block Colors
block_colors <- c(
  "Healthy"  = "#7FC87E", 
  "Bleached" = "#FCC085", 
  "Disease"  = "#E21A1C", 
  "Stressed" = "#3C15D3", 
  "Other"    = "#999999"
)

# 4. Detail Colors
detail_colors <- c(
  "Healthy/Normal"                     = "#7EC87E", 
  "Bleached"                           = "#FCC085", 
  "Disease (Bleaching lesion)"         = "#FC8D3B", 
  "Disease(SCTLD)"                     = "#800025", 
  "Disease(Dark Spot Syndrome)"        = "#BC0025", 
  "Disease (Pathogen-induced disease)" = "#E26E7C", 
  "Stressed(Heat stress)"              = "#756BB0", 
  "Stressed(Acidification stress)"     = "#9E58F5", 
  "Stressed(Combined)"                 = "black",
  "Treated/Recovery"                   = "#D9D9D9", 
  "Other"                              = "#D9D9D9"
)

l1_colors <- brewer.pal(min(8, length(ordered_l1_levels)), "Dark2")
if(length(ordered_l1_levels) > 8) l1_colors <- colorRampPalette(l1_colors)(length(ordered_l1_levels))
names(l1_colors) <- ordered_l1_levels

l2_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(ordered_l2_levels))
names(l2_colors) <- ordered_l2_levels

species_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(all_species))
names(species_colors) <- all_species

ann_colors <- list(
  Species = species_colors,
  Tissue = tissue_colors,
  Ocean = ocean_colors, # Fixed ocean palette.
  Block = block_colors, 
  Detail = detail_colors, 
  Level1 = l1_colors, 
  Level2 = l2_colors
)

heatmap_palette <- colorRampPalette(c("white", "orange", "red", "black"))(100)
global_max_val <- max(mat_final, na.rm = TRUE)
global_breaks <- seq(0, global_max_val, length.out = 100)

# ==============================================================================
# 8. Insert visual gaps between species groups and functional categories.
# ==============================================================================
gaps_col <- which(row_info$Species[-1] != row_info$Species[-nrow(row_info)])
l2_vec <- ann_row_final$Level2
gaps_row <- which(l2_vec[-1] != l2_vec[-length(l2_vec)])

# ==============================================================================
# Project mapping used for display labels; aggregation keys and values are unchanged.
# Use all matching metadata samples to check that each original group has one Project.
project_group_lookup <- meta_clean %>%
  mutate(Row_ID = paste(Species, Ocean_Basin, Block_Group, Detail_Status,
                        Sample_Type, sep = "||")) %>%
  group_by(Row_ID) %>%
  summarise(Project_Count = n_distinct(Project),
            Project = paste(sort(unique(as.character(Project))), collapse = "+"),
            .groups = "drop")
project_group_index <- match(row_info$Row_ID, project_group_lookup$Row_ID)
stopifnot(!anyNA(project_group_index))
project_group_match <- project_group_lookup[project_group_index, , drop = FALSE]
stopifnot(all(project_group_match$Project_Count == 1),
          all(!is.na(project_group_match$Project)),
          all(nzchar(project_group_match$Project)),
          !any(project_group_match$Project == "NA"))
row_info$Project <- project_group_match$Project
display_labels_project <- paste0(row_info$Project, "_", row_info$Species)

# Attach Project to exported source labels while retaining matrix values and keys.
source_data_project <- mat_final
colnames(source_data_project) <- paste0(row_info$Project, "||", colnames(mat_final))
source_group_project <- data.frame(Project = row_info$Project,
                                   ann_col_final,
                                   Display_Label = display_labels_project,
                                   check.names = FALSE)
rownames(source_group_project) <- colnames(source_data_project)
stopifnot(identical(unname(source_data_project), unname(mat_final)))

# ==============================================================================
# 9. Export group and gene annotations and the exact display matrix.
# ==============================================================================
print(">>> Exporting source data...")

# The numerical source matrix is exported below as amg_heatmap_log1p.tsv.
# Matrix rows are genes; matrix columns are ecological groups.
# Compound group keys are explained by the separate group annotation table.

# Export column metadata, including Project and display labels.
# These fields map directly to the annotation bars above the heatmap.
write.csv(source_group_project, file = "group_annotations.csv", quote = FALSE)

# Export row metadata for the functional classification bars.
# The row order matches the exported numerical matrix.
write.csv(ann_row_final, file = "gene_annotations.csv", quote = FALSE)

print(">>> Source data exported.")

# ==============================================================================
# 10. Render and export the final vector PDF.
# ==============================================================================
print(">>> Rendering the final heatmap...")

my_cell_width <- 18  
my_cell_height <- 14 

dynamic_width <- max(25, 8 + (ncol(mat_final) * my_cell_width * 0.045)) 
dynamic_height <- max(20, 5 + (nrow(mat_final) * my_cell_height * 0.045))

pdf(file.path(figures_dir, "amg_heatmap.pdf"), width = dynamic_width, height = dynamic_height)

heatmap_plot <- pheatmap(mat_final,
         color = heatmap_palette, 
         breaks = global_breaks,
         border_color = NA, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE,
         
         annotation_col = ann_col_final, 
         annotation_row = ann_row_final, 
         annotation_colors = ann_colors,
         
         gaps_col = gaps_col, 
         gaps_row = gaps_row,
         
         labels_col = display_labels_project,
         
         show_colnames = TRUE, 
         show_rownames = TRUE,
         
         fontsize = 12,      
         fontsize_row = 10,   
         fontsize_col = 10,   
         
         angle_col = 90,
         
         cellwidth = my_cell_width, 
         cellheight = my_cell_height, 
         
         main = "Global Viral AMGs Landscape (X: Species/Geo, Y: Genes)"
)

dev.off()
print(paste0(">>> Heatmap PDF and source data written."))
# Export the exact log1p matrix displayed in the PDF.
write.table(data.frame(Gene = rownames(source_data_project), source_data_project, check.names = FALSE),
            file = "amg_heatmap_log1p.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
# Save the execution environment alongside the numerical source data.
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

# R-native full-page preview of the same gtable, including all annotation legends.
# Low-resolution screen preview; use the PDF for vector artwork.
ragg::agg_png(file.path(figures_dir, "amg_heatmap_preview.png"), width = dynamic_width,
              height = dynamic_height, units = "in", res = 72,
              background = "white")
grid::grid.newpage()
grid::grid.draw(heatmap_plot$gtable)
dev.off()
writeLines(sprintf("Completed: %d genes x %d groups; %d input genomes x %d samples.",
                   nrow(mat_final), ncol(mat_final), nrow(tpm_df), ncol(tpm_df) - 1))
