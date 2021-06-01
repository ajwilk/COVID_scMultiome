suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(magrittr)
  library(Seurat)
})

analysis_dir <- "data/scATAC/" #Set this to the directory with input metadata files
plot_dir <- "." # Set this to the directory where plots should be saved

sample_metadata <- read_csv(glue("{analysis_dir}/atac_sample_metadata.csv"))
cell_metadata <- read_csv(glue("{analysis_dir}/atac_cell_metadata.csv"))
peak_count_metadata <- read_csv(glue("{analysis_dir}/atac_peak_count_metadata.csv"))
motif_zscore_metadata <- read_csv(glue("{analysis_dir}/atac_motif_zscore_metadata.csv"))

cell_type_colors <- c(
  "#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", 
  "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC", "#90D5E4", 
  "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8", 
  "#6E4B9E", "#0C727C", "#7E1416", "#D8A767", "#3D3D3D")

names(cell_type_colors) <-c(
  "CD14 Monocyte", "CD16 Monocyte", "CD4 T", "CD8 T", "NK", 
  "B", "PB", "Neutrophil", "Developing neutrophil", "DC", 
  "pDC", "Eos/Mast Prog", "Prolif Lymph", "Platelet", "Eryth", 
  "HSPC", "Misc T", "Proliferating", "UNUSED", "UNUSED")
# Lookup for converting fine to broad cell types in ATAC UMAP
# celltype.fine_to_broad[fine_types] returns the corresponding broad_types
celltype.fine_to_broad <- c(
  "B intermediate" = "B","B memory" = "B","B naive" = "B", "Plasmablast" = "PB", "B immature" = "B",
  "CD4 CTL" = "CD4 T", "CD4 Naive" = "CD4 T", "CD4 TCM" = "CD4 T", "CD4 TEM" = "CD4 T",
  "CD8 Naive" = "CD8 T","CD8 TCM" ="CD8 T", "CD8 TEM" = "CD8 T",
  "Proliferating" = "Proliferating", "MAIT" = "Misc T","gdT" = "Misc T", "Treg" = "Misc T","dnT" = "Misc T",
  "NK" = "NK", "NK_CD56bright" = "NK", "NK_CD56" = "NK",
  "CD14 Mono" = "CD14 Monocyte","CD16 Mono" = "CD16 Monocyte",
  "pDC" = "pDC","cDC1" = "DC","cDC2" = "DC","ASDC" = "DC",
  "Platelet" = "Platelet",
  "Eryth" = "Eryth",
  "HSPC" = "HSPC"
)

plot_boxplot <- function(df, value_column, value_label, point_size=1, text_size=NULL) {
  shapes <- c(17, 16)
  comparisons <- list(c("0","1-3"),
                      c("0","4-5"),
                      c("0","6-7"))
  df <- mutate(df, Acuity=ifelse(current_severity_bin=="0", "Healthy", "Acute"))
  ggplot(df, aes_string("current_severity_bin", value_column)) +
    geom_boxplot(alpha = 0.25, outlier.color = NA) +
    geom_point(size = point_size, position = position_jitter(width = 0.25, height=0), 
               aes(color=max_severity_bin, shape=Acuity)) +
    ggpubr::stat_compare_means(mapping = aes(current_severity_bin), comparisons = comparisons, label = "p.signif") +
    scale_color_manual(values=c("0" = "steelblue", "1-3"="green4", "4-5"="yellow3", "6-7"="red", "8"="black")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    scale_shape_manual(values = shapes) +
    guides(fill=FALSE) +
    labs(x = "WHO severity score at draw", y = value_label, color="Peak severity\nscore") + 
    theme_bw() + theme(panel.grid.minor = element_blank(), 
                       panel.grid.major = element_blank(), 
                       strip.background = element_rect(fill = NA, color = NA), 
                       strip.text = element_text(face = "bold"),
                       axis.ticks.x = element_blank(), 
                       axis.text = element_text(color = "black"), 
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(panel.grid.major = element_line(color = "grey", size = 0.25), aspect.ratio = 1,
          text = element_text(size = text_size))
}
plot_cell_percentages <- function(samples, cell_types, ncol=4, point_size=1, plot_celltypes = NULL) {
  tibble(
    Sample=samples,
    cell.type=cell_types
  ) %>%
    group_by(Sample, cell.type) %>%
    summarize(count = n(), .groups="drop") %>%
    right_join(expand(., Sample, cell.type), by=c("Sample", "cell.type")) %>%
    mutate(count=replace_na(count, 0)) %>%
    left_join(sample_metadata, by=c("Sample"="arrow_name")) %>%
    filter(!convalescent, !is_fresh, !is.na(current_severity)) %>%
    group_by(Sample) %>%
    mutate(percent = 100*count/sum(count)) %>%
    filter(is.null(plot_celltypes) | cell.type %in% plot_celltypes) %>%
    plot_boxplot("percent", "Proportion of PBMCs (%)", point_size=point_size) + 
      facet_wrap("cell.type", scales="free_y", ncol=ncol)
}

plot_fine_types <- cell_metadata %>%
  filter(!is_doublet) %>%
  {plot_cell_percentages(.$Sample, .$cell_type)}
ggsave(glue("{plot_dir}/boxplot_fine_celltypes.pdf"), plot_fine_types, width =8.5, height = 11)

plot_broad <- cell_metadata %>%
  filter(!is_doublet) %>%
  {plot_cell_percentages(.$Sample, celltype.fine_to_broad[.$cell_type])}
ggsave(glue("{plot_dir}/boxplot_broad_celltypes.pdf"), plot_broad, width =9, height = 6)

plot_nk <- cell_metadata %>%
  filter(!is_doublet) %>%
  {plot_cell_percentages(.$Sample, celltype.fine_to_broad[.$cell_type], plot_celltypes = "NK")}
ggsave(glue("{plot_dir}/boxplot_nk_cells.pdf"), plot_nk)

peaks <- names(peak_count_metadata)[grep("chr*", names(peak_count_metadata))]
for (p in peaks) {
  plot_peak <- cell_metadata %>%
    select(-current_severity_bin) %>%
    inner_join(peak_count_metadata) %>%
    inner_join(sample_metadata, by=c("Sample"="arrow_name")) %>%
    filter(!is_doublet) %>%
    group_by(Sample, current_severity_bin, max_severity_bin) %>%
    summarize(insertions_per_million = sum(.data[[p]])/sum(nFrags*2)*1e6) %>%
    plot_boxplot("insertions_per_million", "Insertions Per Million", point_size = 4, text_size=20) +
    labs(subtitle = p)
  ggsave(glue("{plot_dir}/boxplot_peak_{p}.pdf"), plot_peak, width=6, height=6)
}

tfs <- names(motif_zscore_metadata) %>% setdiff("Sample")

tf <- "NFKB2_714"
tf_name <- str_extract(tf, "^[^_]*")
plot_tf <- inner_join(motif_zscore_metadata, sample_metadata, by=c("Sample"="arrow_name")) %>%
  plot_boxplot(tf, "ChromVAR Motif Z-score", point_size=4, text_size=20) +
  labs(subtitle=tf_name)
ggsave(glue("{plot_dir}/boxplot_motif_{tf_name}.pdf"), plot_tf, width=6, height=6)


umap_cells <-  cell_metadata %>%
  filter(!is_doublet) %>%
  slice_sample(n=nrow(.))
umap_coords <- umap_cells %>%
  {set_colnames(t(cbind(.$UMAP1, .$UMAP2)), .$cell_id)} %>%
  set_rownames(c("UMAP1", "UMAP2"))


seurat <- CreateSeuratObject(umap_coords)
seurat[["umap"]] <- CreateDimReducObject(embeddings = t(umap_coords), key = "UMAP_", assay="RNA")
seurat$max_severity_bin <- umap_cells %>%
  inner_join(sample_metadata, by=c("Sample"="arrow_name")) %>%
  pull(max_severity_bin)
seurat$cell_type <- recode( 
  celltype.fine_to_broad[umap_cells$cell_type],
  "Proliferating" = " ",
  "Misc T" = " ",
  "HSPC" = " "
)

umap_severity <- DimPlot(seurat, group.by = "max_severity_bin") +
  scale_color_manual(values=c("0" = "steelblue", "1-3"="green4", "4-5"="yellow3", "6-7"="red", "8"="black")) +
  labs(x = "UMAP1", y = "UMAP2") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
umap_severity$layers[[1]]$aes_params$alpha <- .2
umap_severity <- umap_severity + guides(colour = guide_legend(override.aes = list(alpha = 1,size=4))) + NoLegend()
ggsave(glue("{plot_dir}/umap_severity.png"), umap_severity)


umap_cell_type <- DimPlot(seurat, group.by = "cell_type", label = TRUE, repel=TRUE) +
  scale_color_manual(values=cell_type_colors) +
  labs(x = "UMAP1", y = "UMAP2") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
umap_cell_type$layers[[1]]$aes_params$alpha <- .2
umap_cell_type <- umap_cell_type + guides(colour = guide_legend(override.aes = list(alpha = 1,size=4))) + NoLegend()
ggsave(glue("{plot_dir}/umap_cell_type.png"), umap_cell_type)
