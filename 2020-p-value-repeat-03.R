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

library(data.table)
final_dt <- readRDS("final_mutation_regions.rds")
## Add a weight to prob for each donor to get sample-specific prob as Prof.Liu and JingZhang devised.
#system("zcat final_mutation.tsv.gz | cut -f 1 | sort | uniq -c | sed 's/^[ \t]*//g' | sed 's/ /\t/g' > donor_noncoding_mut_freq.tsv")
donor_freq <- fread("donor_noncoding_mut_freq.tsv", header = FALSE)
colnames(donor_freq) <- c("freq", "donor")
donor_freq[, weight := freq / sum(freq)]

final_dt <- merge(final_dt, donor_freq, by = "donor", all.x = TRUE)
final_dt[, prob := prob * weight]
any(is.na(final_dt$prob))
final_dt[, c("freq", "weight") := NULL]

## Get mutation-specific prob
## prob x >= K (K is the mutation freq, so here minus 1)
## ppoibin is used to get Pr(x<K)
final_dt2 <- unique(final_dt[, .(donor, prob, mut_index)])
prob_point <- final_dt2[, .(p_val = 1 - poibin::ppoibin(length(donor) - 1, prob),
                           donor_list = paste(donor, collapse = ",")),
                     by = .(mut_index)]
prob_point
prob_point$mut_index[duplicated(prob_point$mut_index)]

## Set 0 to minimal p value
prob_point[, p_val := ifelse(p_val < .Machine$double.xmin, 
                             .Machine$double.xmin,
                             p_val)]

length(unique(final_dt$mut_index))
length(unique(final_dt$region_midpoint))

prob_point = prob_point[order(p_val)]
openxlsx::write.xlsx(prob_point, file = "PointMutationList.xlsx")

## Get region prob
region_df <- merge(unique(final_dt[, .(region_midpoint, mut_index, donor)]),
                   prob_point[, .(mut_index, p_val)], 
                   by = "mut_index", all.x = TRUE)
region_df

cal_region_p = function(p) {
  1 - cumprod(1 - p)[length(p)]
}

prob_region <- region_df[, .(p_val = cal_region_p(p_val),
                             donor_list = paste(donor, collapse = ","),
                             mutation_list = paste(unique(mut_index), collapse = ",")),
                         by = .(region_midpoint)]

## Set 0 to minimal p value
prob_region[, p_val := ifelse(p_val < .Machine$double.xmin, 
                             .Machine$double.xmin,
                             p_val)]
prob_region = prob_region[order(p_val)]
prob_region

openxlsx::write.xlsx(prob_region, file = "RegionMutationList.xlsx")
