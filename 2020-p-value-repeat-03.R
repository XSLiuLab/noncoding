# Get the mutation frequency
system(
  "cat pred*.tsv.gz | grep -v donor | grep -v mock | cut -f2,3 | sort | uniq -c | sort -n | gzip > final_mutation_freq.tsv.gz"
)

system(
  "cat pred*.tsv.gz | grep -v donor | grep -v mock | gzip > final_mutation.tsv.gz"
)

library(data.table)

freq_dt = fread("final_mutation_freq.tsv.gz", header = FALSE)
head(freq_dt)
freq_dt$chr = sub("[0-9]+ (chr[0-9XYxy]{1,2})", "\\1", freq_dt$V1)
table(freq_dt$chr)

freq_dt$freq = sub("([0-9]+) chr.*", "\\1", freq_dt$V1)
freq_dt$freq = as.integer(freq_dt$freq)

freq_dt = freq_dt[freq > 5]
freq_dt$V1 = NULL
freq_dt$start = freq_dt$V2 - 5
freq_dt$end = freq_dt$V2 + 5
freq_dt$V2 = NULL
freq_dt = freq_dt[, c("chr", "start", "end", "freq"), with = F]
# must be unique
freq_dt$freq = NULL
freq_dt = unique(freq_dt)
saveRDS(freq_dt, file = "recurrent_region.rds")

## Load all mutations
mut_dt = fread("final_mutation.tsv.gz", header = FALSE)
colnames(mut_dt) = c("donor", "chr", "end", "prob")
# must be unique
mut_dt = unique(mut_dt)
mut_dt$start = mut_dt$end
mut_dt = mut_dt[, c("chr", "start", "end", "prob", "donor"), with=F]

str(head(mut_dt))

setkey(freq_dt, chr, start, end)
final_dt = foverlaps(mut_dt, freq_dt, type = "within")
final_dt = final_dt[!is.na(start)]

head(final_dt)
final_dt[, region_midpoint := paste(chr, as.integer(start + 5), sep = ":")]
final_dt$i.start <- NULL
final_dt[, mut_index := paste(chr, as.integer(i.end), sep = ":")]
final_dt$i.end = NULL

head(final_dt)

saveRDS(final_dt, file = "final_mutation_regions.rds")

# gtf_dt = fread("~/zhangjing_20200416/tmp_dat/Homo_sapiens.GRCh37.75.gtf",
#                skip = 5, header = FALSE)

## Add a weight to prob for each donor to get sample-specific prob as Prof.Liu and JingZhang devised.
## TODO


## Get mutation-specific prob
prob_mut <- final_dt[, .(p_val = 1 - poibin::ppoibin(unique(freq) - 1, prob)),
                     by = .(mut_index)]
prob_mut
prob_mut$mut_index[duplicated(prob_mut$mut_index)]
