log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(ggplot2)

ref_bases <- 659924
window_size <- 1e4
snp_files <- list.files("output/060_dnadiff",
                        pattern = "contigs.snps",
                        recursive = TRUE,
                        full.names = TRUE)


# read the data
names(snp_files) <- basename(dirname(snp_files))
snp_results_list <- lapply(snp_files, fread)
snp_results <- rbindlist(snp_results_list, idcol = "assembly", fill = TRUE)

# release memory
rm(snp_results_list)
gc()

# V1 has to be reference position of the SNP
# TRUE means a SNP at this position
setnames(snp_results, "V1", "ref_base")
snp_positions <- dcast(snp_results,
                       ref_base ~ assembly,
                       value.var = "V12",
                       fun.aggregate = function(x) sum(!is.na(x)) > 0)

# set up references
# have to deal with contig names at some point
ref <- data.table(ref_base = 1:ref_bases)

snp_table <- merge(
  ref,
  snp_positions,
  all.x = TRUE,
  by = "ref_base")

snp_totals <- colSums(snp_table[,!1], na.rm = TRUE)
mean_qscores <- -10 * log(snp_totals / ref_bases, 10)

# NA means no SNP
rolling_count <- snp_table[, lapply(.SD,
                   frollsum,
                   n = window_size, 
                   fill = 0,
                   na.rm = TRUE,
                   hasNA = TRUE),
          .SDcols = names(snp_files)]

rolling_count[, ref_pos := 1:ref_bases]

# plot data
pd <- melt(rolling_count, id.vars = "ref_pos")
pd[, qscore := -10 * log(value / window_size, 10)]
pd[is.infinite(qscore), qscore := NA]
pd[, c("version", "assembly_type") := tstrsplit(variable, ".nano-")]


# release memory
rm(rolling_count, snp_table, snp_positions)
gc()

ggplot(pd,
       aes(x = ref_pos, y = qscore)) +
  facet_grid(version ~ assembly_type) +
  xlab("Reference position") +
  ylab(expression(-10~log[10]("SNP rate"))) +
  geom_line()

biggest_comp <- c("guppy_3.4.4.nano-raw", "guppy_6.1.3.nano-raw")
ggplot(pd[variable %in% biggest_comp],
       aes(x = ref_pos, y = qscore)) +
  facet_grid(version ~ assembly_type) +
  xlab("Reference position") +
  ylab(expression(-10~log[10]("SNP rate"))) +
  geom_line() +
  geom_hline(yintercept = mean_qscores["guppy_6.1.3.nano-raw"],
             linetype = 2)

qscore_table <- data.table(assembly = names(mean_qscores), 
           mean_qscores)
qscore_table[, c("version", "assembly_type") := tstrsplit(assembly, ".nano-")]

ggplot(qscore_table,
       aes(y = mean_qscores, x = version, fill = assembly_type)) +
  scale_fill_viridis_d(guide = guide_legend(title = "Flye mode")) +
  xlab(NULL) +
  ylab(expression(-10~log[10]("SNP rate"))) +
  geom_col(position = "dodge")

