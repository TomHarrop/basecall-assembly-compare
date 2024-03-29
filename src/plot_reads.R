log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(extrafont)
library(ggplot2)
library(scales)

summary_files <- snakemake@input[["summaries"]]
plot_file <- snakemake@output[["plot"]]

# dev
summary_files <- list.files("output/010_basecall",
                            pattern = "sequencing_summary.txt",
                            recursive = TRUE,
                            full.names = TRUE)
plot_file <- "yield_so_far.pdf"

# read the data
names(summary_files) <- basename(dirname(summary_files))
bc_results_list <- lapply(summary_files, fread)
bc_results <- rbindlist(bc_results_list, idcol = "fc", fill = TRUE)

# release memory
rm(bc_results_list)
gc()

# draw plot

ggplot(bc_results, 
       aes(x = sequence_length_template,
           y = mean_qscore_template)) +
  facet_wrap(~ fc) +
  geom_hex()


gp <- ggplot(bc_results, 
             aes(x = sequence_length_template,
                 y = mean_qscore_template)) +
  facet_wrap(~ fc) +
  theme_minimal(base_family = "Lato") +
  ylab("Mean Q score") + xlab("Read length") +
  scale_fill_viridis_c(
    # trans = log_trans(),
    # breaks = log_breaks(),
    guide = guide_colourbar(
      title = "Number of reads"
    )) +
  scale_x_continuous(trans = log_trans(base = 4),
                     breaks = trans_breaks(function(x) log(x, 4),
                                           function(x) 4^x)) +
  geom_hex() +
  geom_hline(yintercept = 9,  linetype = 2)
gp

ggplot(bc_results,
       aes(x = mean_qscore_template)) +
  theme_minimal(base_family = "Lato") +
  facet_grid(fc ~ .) +
  xlab("Mean Q score") +
  geom_density() +
  geom_vline(xintercept = 9, linetype = 2)



# scale_fill_viridis_c(
#   trans = log_trans(),
#   breaks = log_breaks(),
#   guide = guide_colourbar(
#     title = "Number of reads"
#   )) +
scale_x_continuous(trans = log_trans(base = 4),
                   breaks = trans_breaks(function(x) log(x, 4),
                                         function(x) 4^x)) +
  geom_point()

ggsave(plot_file,
       gp,
       width = 10,
       height = 7.5,
       units = "in",
       device=cairo_pdf)

# log
sessionInfo()