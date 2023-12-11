library(readr)

args = commandArgs(trailingOnly=TRUE)

## debug!
# workdir = "~/vol/pico/depienne/minion_analyses"
# setwd(workdir)
# getwd()
# args = c(
#     "results/004_20221019_FAME3_FGF14/ont-guppy_6.4.6/sup/None/read_counts/input.read_counts.tsv",
#     "results/004_20221019_FAME3_FGF14/ont-guppy_6.4.6/sup/None/read_counts/filtering.read_counts.tsv",
#     "results/004_20221019_FAME3_FGF14/ont-guppy_6.4.6/sup/None/read_counts/repeat_region_spanning.read_counts.tsv",
#     "results/004_20221019_FAME3_FGF14/ont-guppy_6.4.6/sup/None/read_counts/remove_flanking_regions.read_counts.tsv",
#     "results/004_20221019_FAME3_FGF14/ont-guppy_6.4.6/sup/None/read_counts/allele.read_counts.tsv",
#     "results/004_20221019_FAME3_FGF14/ont-guppy_6.4.6/sup/None/read_counts/long.read_counts.tsv"
# )

print(args)
print(length(args))

tsv_files = args[1:length(args)-1]
tsv_files

output_tsv = args[length(args)]
output_tsv

tsv_list <- lapply(tsv_files, read_tsv)
d <- do.call(rbind, tsv_list)
d

write_tsv(d, output_tsv)

