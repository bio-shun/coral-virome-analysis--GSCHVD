######################################################################################

# Portable command-line contract:
#   Rscript run_lyso.R [input_directory] [output_directory]
# The defaults are the package's inputs/ directory and package root.  All
# output files are written below results/ and figures/ in the output directory.
lyso_cli_args <- commandArgs(trailingOnly = TRUE)
lyso_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
lyso_package_dir <- if (length(lyso_file_arg)) {
  dirname(normalizePath(sub("^--file=", "", lyso_file_arg[1]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
lyso_input_dir <- if (length(lyso_cli_args) >= 1) {
  normalizePath(lyso_cli_args[1], mustWork = TRUE)
} else {
  file.path(lyso_package_dir, "inputs")
}
lyso_output_dir <- if (length(lyso_cli_args) >= 2) {
  normalizePath(lyso_cli_args[2], mustWork = FALSE)
} else {
  lyso_package_dir
}
lyso_results_dir <- file.path(lyso_output_dir, "results")
lyso_figures_dir <- file.path(lyso_output_dir, "figures")
dir.create(lyso_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(lyso_figures_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Script: run_lyso.R — phage lifestyle abundance and nested figure
# 目标: 
#   1. [计算] 严格计算 (Fragment独立 + Ocean维度)
#   2. [配色] 严格使用用户提供的 Ocean/Tissue/Block/Detail 配色方案
#   3. [绘图] Species -> Compartment -> Ocean -> Block -> Detail 嵌套分面
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggh4x) 
library(cowplot)

# ==============================================================================
# PART 1: Strict Data Calculation 
# ==============================================================================
message(">>> [PART 1] 开始执行严格计算 (Fragment + Ocean)...")

# 1. 读取数据
lyso_input_files <- c("24083_metagenome_TPM.tsv", "sample_info_metagenome.tsv", "lyso.tsv")
lyso_missing <- lyso_input_files[!file.exists(file.path(lyso_input_dir, lyso_input_files))]
if (length(lyso_missing)) {
  stop("Missing input file(s): ", paste(lyso_missing, collapse = ", "))
}

tpm_df <- read_tsv(file.path(lyso_input_dir, "24083_metagenome_TPM.tsv"), show_col_types = FALSE)
meta_df <- read_tsv(file.path(lyso_input_dir, "sample_info_metagenome.tsv"), show_col_types = FALSE)
lyso_df <- read_tsv(file.path(lyso_input_dir, "lyso.tsv"), show_col_types = FALSE)

# 2. 清洗 Lifestyle
lyso_clean <- lyso_df %>% 
  select(Genome, Lifestyle) %>%
  mutate(Lifestyle = case_when(
    grepl("Lysogenic", Lifestyle, ignore.case = T) ~ "Lysogenic", 
    TRUE ~ "Unknown" 
  ))

# 3. 清洗元数据
if(!"Sample" %in% colnames(meta_df) && "SampleID" %in% colnames(meta_df)) { 
  meta_df <- meta_df %>% rename(Sample = SampleID) 
}

meta_clean <- meta_df %>% 
  filter(Sample %in% colnames(tpm_df)) %>% 
  filter(!grepl("DSS", Status, ignore.case = TRUE)) %>% # 过滤 DSS
  mutate(
    # 1. Compartment (简写，用于分面标签)
    Compartment = case_when(
      grepl("Mucus", Sample_Type, ignore.case = T) ~ "Mucus", 
      grepl("Tissue", Sample_Type, ignore.case = T) ~ "Tissue",   
      grepl("Fragment", Sample_Type, ignore.case = T) ~ "Fragment", 
      grepl("Skeleton", Sample_Type, ignore.case = T) ~ "Skeleton", 
      TRUE ~ NA_character_
    ),
    
    # 2. Ocean
    Ocean_Basin = if("Ocean_Basin" %in% colnames(.)) Ocean_Basin else "Unknown",
    Ocean_Basin = ifelse(is.na(Ocean_Basin) | Ocean_Basin == "", "Unknown", Ocean_Basin),
    
    # 3. Block Group
    Block_Group = case_when(
      grepl("Healthy", Status, ignore.case = TRUE) ~ "Healthy", 
      grepl("Bleached", Status, ignore.case = TRUE) ~ "Bleached", 
      grepl("Disease|SCTLD|White Plague|Black Band", Status, ignore.case = TRUE) ~ "Disease", 
      grepl("Stress", Status, ignore.case = TRUE) ~ "Stressed", 
      TRUE ~ "Other"
    ),
    
    # 4. Simplified Status
    Simplified_Status = case_when(
      Block_Group == "Healthy" ~ "Healthy", 
      TRUE ~ Status
    ),
    Species = Coral_Host_Species
  ) %>% 
  filter(!is.na(Compartment))

# 4. 数据合并与 TPM 统计
target_genomes <- lyso_clean$Genome
tpm_long <- tpm_df %>% 
  filter(Genome %in% target_genomes) %>% 
  pivot_longer(-Genome, names_to = "Sample", values_to = "TPM")

plot_data_raw <- tpm_long %>% 
  inner_join(lyso_clean, by = "Genome") %>% 
  inner_join(meta_clean, by = "Sample") %>%
  group_by(Sample, Compartment, Ocean_Basin, Block_Group, Simplified_Status, Species, Lifestyle) %>%
  summarise(Total_TPM = sum(TPM, na.rm=TRUE), .groups = "drop")

# 5. 补全缺失值与计算均值
plot_data_filled <- plot_data_raw %>%
  complete(nesting(Sample, Compartment, Ocean_Basin, Block_Group, Simplified_Status, Species), 
           Lifestyle = c("Lysogenic", "Unknown"), 
           fill = list(Total_TPM = 0))

plot_data_final <- plot_data_filled %>% 
  group_by(Sample) %>% 
  mutate(Sample_Total = sum(Total_TPM)) %>%
  filter(Sample_Total > 0) %>% 
  mutate(Relative_Abundance = Total_TPM / Sample_Total * 100) %>% 
  ungroup()

plot_data_mean <- plot_data_final %>%
  group_by(Compartment, Ocean_Basin, Block_Group, Simplified_Status, Species, Lifestyle) %>%
  summarise(Mean_Rel_Abund = mean(Relative_Abundance), .groups = "drop")

# 保存中间结果
write_csv(plot_data_mean, file.path(lyso_results_dir, "phage_lifestyle_stats.csv"))
message(">>> [PART 1] 完成！数据已保存。")


# ==============================================================================
# PART 2: Visualization (Strict Colors)
# ==============================================================================
message(">>> [PART 2] 开始绘图 (应用指定配色)...")

plot_data_mean <- read_csv(file.path(lyso_results_dir, "phage_lifestyle_stats.csv"), show_col_types = FALSE)

# --- 1. 颜色定义 (严格按照你的要求) ---

lyso_colors <- c("Lysogenic" = "#984EA3", "Unknown" = "#E0E0E0") 

# [User Defined] Ocean Colors
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

# [User Defined] Tissue Colors (原始定义)
tissue_colors_raw <- c(
  "Coral tissue"   = "#D95F02", 
  "Coral mucus"    = "#1B9E77", 
  "Coral skeleton" = "#8D40CA", 
  "Coral fragment" = "#E72989", 
  "Other"          = "#999999"
)

# [Mapping] 将 Tissue Colors 映射到数据中的 Compartment 简写
# 确保 Fragment 使用 #E72989
comp_colors <- c(
  "Mucus"    = tissue_colors_raw[["Coral mucus"]],
  "Tissue"   = tissue_colors_raw[["Coral tissue"]],
  "Skeleton" = tissue_colors_raw[["Coral skeleton"]],
  "Fragment" = tissue_colors_raw[["Coral fragment"]]
)

# [User Defined] Block Colors
block_colors <- c(
  "Healthy"  = "#7FC87E", 
  "Bleached" = "#FCC085", 
  "Disease"  = "#E21A1C", 
  "Stressed" = "#3C15D3", 
  "Other"    = "#999999"
)

# [User Defined] Detail Colors
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

# 补全缺失的 Detail 颜色 (防止报错)
all_details <- unique(plot_data_mean$Simplified_Status)
missing_details <- setdiff(all_details, names(detail_colors))
if(length(missing_details) > 0) {
  temp_cols <- rep("grey90", length(missing_details))
  names(temp_cols) <- missing_details
  detail_colors <- c(detail_colors, temp_cols)
}

# --- 2. 排序与因子化 ---
species_levels <- sort(unique(plot_data_mean$Species))
compartment_levels <- c("Mucus", "Tissue", "Fragment", "Skeleton") 
ocean_levels <- sort(unique(plot_data_mean$Ocean_Basin))
block_levels <- c("Healthy", "Bleached", "Disease", "Stressed", "Other")
status_raw <- unique(plot_data_mean$Simplified_Status)
status_levels <- c("Healthy", sort(status_raw[status_raw != "Healthy"]))

plot_data_sorted <- plot_data_mean %>%
  mutate(
    Species = factor(Species, levels = species_levels),
    Compartment = factor(Compartment, levels = compartment_levels),
    Ocean_Basin = factor(Ocean_Basin, levels = ocean_levels),
    Block_Group = factor(Block_Group, levels = block_levels),
    Simplified_Status = factor(Simplified_Status, levels = status_levels)
  ) %>%
  arrange(Species, Compartment, Ocean_Basin, Block_Group, Simplified_Status)

# --- 3. 准备嵌套分面颜色条 ---
unique_combs <- plot_data_sorted %>% 
  select(Species, Compartment, Ocean_Basin, Block_Group, Simplified_Status) %>% 
  distinct()

cols_layer1 <- rep("white", nrow(unique_combs %>% select(Species) %>% distinct()))
cols_layer2 <- comp_colors[as.character(unique_combs %>% group_by(Species, Compartment) %>% slice(1) %>% pull(Compartment))]
cols_layer3 <- ocean_colors[as.character(unique_combs %>% group_by(Species, Compartment, Ocean_Basin) %>% slice(1) %>% pull(Ocean_Basin))]
cols_layer4 <- block_colors[as.character(unique_combs %>% group_by(Species, Compartment, Ocean_Basin, Block_Group) %>% slice(1) %>% pull(Block_Group))]
cols_layer5 <- detail_colors[as.character(unique_combs$Simplified_Status)]

# NA 处理
cols_layer2[is.na(cols_layer2)] <- "grey90"
cols_layer3[is.na(cols_layer3)] <- "grey90"
cols_layer4[is.na(cols_layer4)] <- "grey90"
cols_layer5[is.na(cols_layer5)] <- "grey90"

# --- 4. 绘图 ---
my_big_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 11, color="black"),
    axis.text.y = element_text(size = 12, color="black"),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none", 
    strip.text = element_text(face = "bold", size = 9, color = "black"),
    panel.spacing = unit(0.1, "cm")
  )

p_main <- ggplot(plot_data_sorted, aes(x = Species, y = Mean_Rel_Abund)) +
  geom_bar(aes(fill = Lifestyle), stat = "identity", position = "stack", width = 0.9, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = lyso_colors, name = "Viral Lifestyle") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  
  # 5层嵌套分面
  facet_nested(
    ~ Species + Compartment + Ocean_Basin + Block_Group + Simplified_Status, 
    scales = "free", space = "free",
    nest_line = element_line(linetype = 1, color="grey50"),
    strip = strip_nested(
      background_x = list(
        element_rect(fill = cols_layer1, color = "black"), # Species
        element_rect(fill = cols_layer2, color = "black"), # Compartment (Colors Mapped)
        element_rect(fill = cols_layer3, color = "black"), # Ocean (User Colors)
        element_rect(fill = cols_layer4, color = "black"), # Block (User Colors)
        element_rect(fill = cols_layer5, color = "black")  # Detail (User Colors)
      ),
      by_layer_x = FALSE
    )
  ) +
  labs(title = "Phage Lifestyle (Exact Colors)", y = "Relative Abundance (%)", x = "") +
  my_big_theme

# --- 5. 图例 ---
get_legend_plot <- function(mapping_aes, scale_func, title) {
  p <- ggplot(plot_data_sorted, mapping_aes) +
    geom_tile(aes(x = 1, y = 1), alpha = 1) + 
    scale_func +
    guides(fill = guide_legend(title = title, ncol = 1, override.aes = list(size = 5))) +
    theme_void() +
    theme(legend.position = "right", 
          legend.title = element_text(face="bold", size=12),
          legend.text = element_text(size=10),
          legend.key.size = unit(0.8, "cm"))
  return(get_legend(p))
}

leg_lifestyle <- get_legend_plot(aes(fill = Lifestyle), scale_fill_manual(values = lyso_colors), "Viral Lifestyle")
leg_comp <- get_legend_plot(aes(fill = Compartment), scale_fill_manual(values = comp_colors), "Compartment")
leg_ocean <- get_legend_plot(aes(fill = Ocean_Basin), scale_fill_manual(values = ocean_colors), "Ocean Basin")
leg_block <- get_legend_plot(aes(fill = Block_Group), scale_fill_manual(values = block_colors), "Block Group")
leg_detail <- get_legend_plot(aes(fill = Simplified_Status), scale_fill_manual(values = detail_colors), "Detail Status")

# --- 6. 组合与导出 ---
message(">>> 正在组合图片...")

legends_col <- plot_grid(
  leg_lifestyle, 
  leg_comp, 
  leg_ocean, 
  leg_block, 
  leg_detail,
  ncol = 1, align = "v", 
  rel_heights = c(0.4, 0.5, 0.8, 0.7, 2.0)
)

final_plot <- plot_grid(
  p_main, legends_col,
  ncol = 2, rel_widths = c(1, 0.3)
)

n_cols <- nrow(plot_data_sorted)
calc_width <- ifelse(n_cols > 60, 22, 18) 

ggsave(file.path(lyso_figures_dir, "phage_lifestyle.pdf"), final_plot,
       width = calc_width, height = 12, limitsize = FALSE)
message(">>> 全部完成！颜色已完全同步。")

# Save a compact R-native audit object alongside the reproducible results.
lyso_audit_dir <- file.path(lyso_results_dir, "audit")
dir.create(lyso_audit_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(mean_abundance = plot_data_mean, sorted_abundance = plot_data_sorted,
             sample_abundance = plot_data_final, sample_tpm = plot_data_filled,
             metadata = meta_clean, lifestyle = lyso_clean,
             final_plot = final_plot, main_plot = p_main,
             page_width_inches = calc_width, page_height_inches = 12),
        file.path(lyso_audit_dir, "analysis_objects.rds"))
writeLines(capture.output(sessionInfo()), file.path(lyso_audit_dir, "R_sessionInfo.txt"))
ggsave(file.path(lyso_figures_dir, "phage_lifestyle_preview.png"), final_plot, width = calc_width,
       height = 12, dpi = 140, device = ragg::agg_png, limitsize = FALSE,
       bg = "white")
cat(sprintf("Completed: %d genomes; %d TPM samples; %d metadata samples; %d nonzero-denominator samples; %d plotted groups.\n",
            nrow(tpm_df), ncol(tpm_df) - 1, nrow(meta_clean),
            n_distinct(plot_data_final$Sample),
            nrow(distinct(plot_data_mean, Compartment, Ocean_Basin,
                          Block_Group, Simplified_Status, Species))))
