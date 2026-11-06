# Coral virome alpha diversity: render the final DNA and RNA source values.
# Usage: Rscript scripts/plot_alpha_diversity.R [results_directory] [figures_directory]
# This plotting entry point does not recalculate diversity or statistical tests.
options(warn = 1)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
script_flag <- grep("^--file=", commandArgs(), value = TRUE)
bundle_root <- if (length(script_flag)) dirname(dirname(normalizePath(sub("^--file=", "", script_flag[[1]]), winslash = "/", mustWork = TRUE))) else getwd()
input_dir <- if (length(args) >= 1) args[[1]] else file.path(bundle_root, "results")
output_dir <- if (length(args) >= 2) args[[2]] else file.path(bundle_root, "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"
# Use the available R raster device for accurate Arial metrics during layout.
cowplot::set_null_device("agg")

status_order <- c(
  "Healthy",
  "Bleached",
  "Diseased(SCTLD)",
  "Diseased(black band disease)",
  "Diseased(white plague disease)",
  "Diseased(Dark Spot Syndrome)",
  "Diseased(Bleaching lesion)",
  "Diseased(Pathogen-induced)",
  "Stressed(Heat stress)",
  "Stressed(Acidification stress)",
  "Stressed(Combined)",
  "Stressed(Cold stress)",
  "Other"
)

status_colors <- c(
  "Healthy" = "#33a02c",
  "Bleached" = "#e31a1c",
  "Diseased(SCTLD)" = "#542788",
  "Diseased(black band disease)" = "#40004b",
  "Diseased(white plague disease)" = "#9970ab",
  "Diseased(Dark Spot Syndrome)" = "#762a83",
  "Diseased(Bleaching lesion)" = "#FD8D3C",
  "Diseased(Pathogen-induced)" = "#8073ac",
  "Stressed(Heat stress)" = "#ff7f00",
  "Stressed(Acidification stress)" = "#fdbf6f",
  "Stressed(Combined)" = "#b15928",
  "Stressed(Cold stress)" = "#1f78b4",
  "Other" = "#999999"
)

source_files <- c(
  DNA = file.path(input_dir, "SourceData_DNA_Alpha_Diversity_ProjectAware.csv"),
  RNA = file.path(input_dir, "SourceData_RNA_Alpha_Diversity_ProjectAware.csv")
)
observed_statuses <- status_order[
  status_order %in% c(
    as.character(read.csv(source_files[["DNA"]], check.names = FALSE, stringsAsFactors = FALSE)$Status_Group),
    as.character(read.csv(source_files[["RNA"]], check.names = FALSE, stringsAsFactors = FALSE)$Status_Group)
  )
]

read_alpha_source <- function(analysis_type) {
  path <- file.path(
    input_dir,
    paste0("SourceData_", analysis_type, "_Alpha_Diversity_ProjectAware.csv")
  )
  if (!file.exists(path)) {
    stop("Source data file not found: ", path)
  }

  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) %>%
    mutate(
      Status_Group = factor(Status_Group, levels = status_order),
      Display_Species = as.character(Display_Species)
    )
}

# Assign explicit x positions within each species. This makes the box centre
# and the point cloud centre identical even when some Status levels are absent.
prepare_plot_data <- function(source_data) {
  display_order <- source_data %>%
    distinct(Ocean_Basin, Clean_Type, Species, Analysis_Project, Display_Species) %>%
    arrange(Ocean_Basin, Clean_Type, Species, Analysis_Project) %>%
    pull(Display_Species) %>%
    unique()

  long_data <- source_data %>%
    pivot_longer(
      cols = c(Shannon, Chao1),
      names_to = "Index",
      values_to = "Value"
    ) %>%
    mutate(
      Display_Species = factor(Display_Species, levels = display_order),
      Status_Group = factor(Status_Group, levels = status_order)
    )

  status_positions <- long_data %>%
    distinct(Display_Species, Status_Group) %>%
    filter(!is.na(Status_Group)) %>%
    arrange(Display_Species, Status_Group) %>%
    group_by(Display_Species) %>%
    mutate(
      n_status = n(),
      offset = if (n() == 1) 0 else seq(-0.38, 0.38, length.out = n())
    ) %>%
    ungroup()

  long_data %>%
    left_join(status_positions, by = c("Display_Species", "Status_Group")) %>%
    mutate(
      x_base = as.numeric(Display_Species),
      x_pos = x_base + offset,
      box_group = interaction(Display_Species, Status_Group, drop = TRUE)
    )
}

make_alpha_plot <- function(source_data, analysis_type) {
  plot_data <- prepare_plot_data(source_data)
  display_order <- levels(plot_data$Display_Species)

  ggplot(
    plot_data,
    aes(x = x_pos, y = Value, fill = Status_Group, group = box_group)
  ) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      alpha = 0.72,
      linewidth = 0.38,
      colour = "black",
      position = "identity"
    ) +
    geom_point(
      position = position_jitter(width = 0.035, height = 0, seed = 123),
      size = 1.8,
      alpha = 0.78,
      shape = 21,
      colour = "black",
      stroke = 0.22
    ) +
    facet_grid(
      Index ~ Ocean_Basin + Clean_Type,
      scales = "free",
      space = "free_x",
      drop = TRUE,
      labeller = labeller(
        Ocean_Basin = label_wrap_gen(12),
        Clean_Type = label_wrap_gen(14)
      )
    ) +
    scale_x_continuous(
      breaks = seq_along(display_order),
      labels = display_order,
      expand = expansion(mult = c(0.08, 0.08))
    ) +
    scale_fill_manual(
      values = status_colors,
      limits = observed_statuses,
      drop = FALSE,
      na.translate = FALSE,
      guide = guide_legend(nrow = 2, byrow = TRUE)
    ) +
    labs(
      x = "Coral taxon",
      y = paste0(analysis_type, " viral alpha diversity"),
      fill = "Status"
    ) +
    theme_classic(base_size = 9, base_family = font_family) +
    theme(
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.42),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.38, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text.x = element_text(
        angle = 55,
        hjust = 1,
        vjust = 1,
        face = "italic",
        size = 7.6,
        colour = "black"
      ),
      axis.text.y = element_text(size = 8, colour = "black"),
      axis.title = element_text(size = 9, colour = "black"),
      strip.text = element_text(size = 7.2, face = "bold", colour = "black"),
      strip.background = element_rect(fill = "grey96", colour = "black", linewidth = 0.35),
      legend.position = "bottom",
      legend.title = element_text(size = 8.5),
      legend.text = element_text(size = 7.5),
      legend.key.size = grid::unit(0.42, "cm"),
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.box.just = "center",
      legend.spacing.x = grid::unit(0.22, "cm"),
      plot.margin = margin(5.5, 7, 5.5, 5.5)
    )
}

dna <- read_alpha_source("DNA")
rna <- read_alpha_source("RNA")

p_dna <- make_alpha_plot(dna, "DNA")
p_rna <- make_alpha_plot(rna, "RNA")

max_species <- max(
  n_distinct(dna$Display_Species),
  n_distinct(rna$Display_Species)
)
fig_width <- max(16, 0.75 * max_species + 7)

save_plot <- function(plot, filename, width, height) {
  # Retain editable vector text in both formats without changing the layout.
  svglite::svglite(file = file.path(output_dir, sub('\\.pdf$', '.svg', filename)),
                   width = width, height = height, system_fonts = list(sans = font_family))
  print(plot)
  grDevices::dev.off()
  grDevices::cairo_pdf(filename = file.path(output_dir, filename), width = width, height = height, family = font_family)
  print(plot)
  grDevices::dev.off()
}

save_plot(p_dna, "Alpha_Diversity_DNA_ProjectAware_Refined_Wide.pdf", fig_width, 8.2)
save_plot(p_rna, "Alpha_Diversity_RNA_ProjectAware_Refined_Wide.pdf", fig_width, 8.2)

# Build the shared legend explicitly from the union of observed Status levels.
# This prevents patchwork from inheriting a legend layout from only one assay.
legend_data <- data.frame(
  Status_Group = factor(observed_statuses, levels = observed_statuses),
  x = 1,
  y = 1
)

legend_plot <- ggplot(
  legend_data,
  aes(x = x, y = y, fill = Status_Group)
) +
  geom_point(shape = 22, size = 3.2, colour = "black", stroke = 0.25) +
  scale_fill_manual(
    values = status_colors,
    limits = observed_statuses,
    drop = FALSE,
    guide = guide_legend(
      title = "Status",
      nrow = 2,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(shape = 22, size = 3.2, alpha = 1)
    )
  ) +
  theme_void(base_size = 9, base_family = font_family) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 7.5),
    legend.key.size = grid::unit(0.42, "cm"),
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.box.just = "center",
    legend.spacing.x = grid::unit(0.22, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(0, 0, 0, 0)
  )

shared_legend <- cowplot::get_legend(legend_plot)

combined <- cowplot::plot_grid(
  p_dna + theme(legend.position = "none"),
  p_rna + theme(legend.position = "none"),
  shared_legend,
  ncol = 1,
  align = "v",
  axis = "tb",
  rel_heights = c(1, 1, 0.20)
)

save_plot(combined, "Alpha_Diversity_ProjectAware_Final_Refined_Wide.pdf", fig_width, 17.1)

# PNG previews are generated by R for visual checking only.
ggsave(
  filename = file.path(output_dir, "Alpha_Diversity_ProjectAware_Final_Refined_Wide_preview.png"),
  plot = combined,
  width = fig_width,
  height = 17.1,
  units = "in",
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

writeLines(
  c(
    "Input: per-assay source tables generated by run_alpha_diversity.R.",
    "Output: editable SVG/PDF figures and PNG screen previews.",
    paste0("Figure width (in): ", format(fig_width, digits = 4)),
    "Point and box positions were explicitly aligned by species and status.",
    "The wide layout uses wider panels and expanded within-species status spacing.",
    "DNA and RNA values and statistics are calculated from their respective count matrices."
  ),
  con = file.path(output_dir, "Alpha_Diversity_ProjectAware_Refined_plot_log.txt")
)

message("Refined alpha-diversity figures written to: ", output_dir)

# The requested vector deliverable is the original-layout Cairo PDF above.
# DNA-only preview for panel-level visual checking.
ggsave(filename = file.path(output_dir, "Alpha_Diversity_DNA_24083_preview.png"), plot = p_dna, width = fig_width, height = 8.2, units = "in", dpi = 300, bg = "white", limitsize = FALSE)
