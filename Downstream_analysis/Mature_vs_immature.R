library(tidyverse)
library(ggpubr)



load_counts <- function(path) {
  counts <- read.delim(path)
  colnames(counts) <- gsub("X", "", colnames(counts)) |> gsub("\\.", "-", x = _)
  counts <- counts[!duplicated(counts$Geneid), ]
  rownames(counts) <- counts$Geneid
  counts[, -c(1, 2)]
}

calculate_ratios <- function(counts,
                             mature   = c("18S", "5.8S", "28S"),
                             immature = c("5'_ETS", "3'_ETS", "ITS1", "ITS2", "Enhancer_Repeats")) {
  data.frame(
    Sample         = colnames(counts),
    Mature_Total   = colSums(counts[rownames(counts) %in% mature,   , drop = FALSE]),
    Immature_Total = colSums(counts[rownames(counts) %in% immature, , drop = FALSE])
  ) |>
    mutate(Ratio = Immature_Total / Mature_Total)
}

assign_conditions <- function(plot_data,
                              conditions = list(DMSO = 1:3, ISD = 4:6)) {
  stopifnot("Not enough samples for condition assignment" =
              nrow(plot_data) >= max(unlist(conditions)))

  plot_data$Condition <- NA_character_
  for (cond in names(conditions)) {
    plot_data$Condition[conditions[[cond]]] <- cond
  }

  plot_data |>
    filter(!is.na(Condition)) |>
    mutate(Condition = factor(Condition, levels = names(conditions)))
}

plot_ratio <- function(plot_data,
                       colors = c(DMSO = "grey70", ISD = "firebrick")) {
  ggplot(plot_data, aes(x = Condition, y = Ratio, fill = Condition)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.1, size = 3, color = "black", alpha = 0.7) +
    stat_compare_means(method = "t.test", label.x.npc = "center",
                       vjust = -1, size = 5) +
    scale_fill_manual(values = colors) +
    theme_classic() +
    labs(
      title = "Ratio of Immature / Mature rRNA Reads (TT-seq)",
      y     = "Ratio (Immature / Mature)",
      x     = "Condition"
    ) +
    theme(
      axis.text    = element_text(size = 12, color = "black"),
      axis.title   = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )
}

# ---------------------------------------------------------
# MAIN
# ---------------------------------------------------------

p <- load_counts("TTseq_counts_processed.txt") |>
  calculate_ratios() |>
  assign_conditions() |>
  plot_ratio()

print(p)
ggsave("ttseq_immature_mature_ratio.pdf", plot = p, width = 5, height = 6)