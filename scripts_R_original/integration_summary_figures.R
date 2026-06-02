\
suppressPackageStartupMessages({
  library(tidyverse)
  source("scripts/00_setup/pipeline_utils.R")
})

cfg <- read_config()
make_dirs(cfg)

same_summary <- file.path(cfg$outputs$tables, "RNA_ATAC_integrated_same_time_summary.csv")
lag_summary <- file.path(cfg$outputs$tables, "RNA_ATAC_integrated_lagged_summary.csv")

if (file.exists(same_summary)) {
  ss <- read_csv(same_summary, show_col_types = FALSE)
  p <- ggplot(ss, aes(x = time, y = n, fill = integrated_class)) +
    geom_col(position = "stack") +
    theme_bw() +
    labs(title = "Same-time RNA-ATAC integration classes", x = "Time", y = "Gene count")
  ggsave(file.path(cfg$outputs$figures, "RNA_ATAC_same_time_class_counts.png"), p, width = 8, height = 5, dpi = 300)
}

if (file.exists(lag_summary)) {
  ls <- read_csv(lag_summary, show_col_types = FALSE) %>%
    mutate(pair = paste(time_atac, "ATAC →", time_rna, "RNA"))
  p <- ggplot(ls, aes(x = pair, y = n, fill = lag_class)) +
    geom_col(position = "stack") +
    theme_bw() +
    labs(title = "Lagged RNA-ATAC integration classes", x = "Lag pair", y = "Gene count") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(cfg$outputs$figures, "RNA_ATAC_lagged_class_counts.png"), p, width = 8, height = 5, dpi = 300)
}

message("Saved integration summary figures where inputs were available.")
