log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(ggplot2)

window_size <- 1e5
vcf_files <- list.files("output/065_minimap-snps",
                        pattern = "out.snps.vcf",
                        recursive = TRUE,
                        full.names = TRUE)
fai_file <- 'output/GCF_003254395.2_Amel_HAv3.1_genomic.fna.fai'

# read the data
fai <- fread(fai_file, col.names = c("#CHROM",
                                     "chr_length",
                                     "offset",
                                     "basesperline",
                                     "bytesperline"))

names(vcf_files) <- basename(dirname(vcf_files))
vcf_results_list <- lapply(vcf_files,
                           fread,
                           skip = "#CHROM")
vcf_results <- rbindlist(vcf_results_list,
                         idcol = "assembly", fill = TRUE)

# release memory
rm(vcf_results_list)
gc()

# V1 has to be reference position of the SNP
# TRUE means a SNP at this position
snp_positions <- dcast(vcf_results,
                       `#CHROM` + POS ~ assembly,
                       value.var = "ALT",
                       fun.aggregate = function(x) sum(!is.na(x)) > 0)

# set up references
# have to deal with contig names at some point
# this will use a lot of RAM
ref <- fai[startsWith(`#CHROM`, "NC_"),
           .(POS = c(1:chr_length)), by = `#CHROM`]

snp_table <- merge(
  ref,
  snp_positions,
  all.x = TRUE,
  by = c("#CHROM", "POS"))

snp_totals <- colSums(snp_table[,!1:2], na.rm = TRUE)
mean_qscores <- -10 * log(snp_totals / ref[, sum(POS)], 10)

# NA means no SNP
rolling_count <- snp_table[, lapply(.SD,
                                    frollsum,
                                    n = window_size, 
                                    fill = 0,
                                    na.rm = TRUE,
                                    hasNA = TRUE),
                           .SDcols = names(vcf_files),
                           by = `#CHROM`]

rolling_count[, POS := ref$POS]

# release memory
rm(snp_table, snp_positions)
gc()


# plot data
pd <- melt(rolling_count,
           id.vars = c("#CHROM", "POS"),
           variable.name = "assembly",
           value.name = "snp_count")

# release memory
rm(rolling_count)
gc()

# tidy data
setnames(pd, "POS", "window_start")
pd[, window_stop := window_start + window_size]
pd[, qscore := -10 * log(snp_count / window_size, 10)]
pd[is.infinite(qscore), qscore := NA]

# write output
