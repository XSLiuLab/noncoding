# New Analysis BEGIN ------------------------------------------------------
# 计算1Mbp区间内突变频率与遗传、表观注释值的相关性
# Coding, Noconding, Promoter ...
# Coding为外显子区域

setwd("~/New_Analysis/")
source("functions.R")
library(data.table)
Sys.setenv(PATH=paste0("/home/zhangjing/bedtools2/bin:", Sys.getenv("PATH")))

# 预处理数据 -------------------------------------------------------------------
options(scipen = 30)
original_mut = "../icgc_data/预处理所有突变.tsv"    # 经过预处理过的突变记录文件

if (F) {
  mut <- fread(original_mut, data.table = F)
  pos_mut <- mut[mut$mut == "single base substitution",]
  ###以bed格式输出所有单碱基突变
  mut_bed <- pos_mut[,3:5]
  mut_bed$start <- mut_bed$start - 1L
  mut_bed$start = as.integer(mut_bed$start)
  mut_bed$end = as.integer(mut_bed$end)
  write.table(mut_bed, "single_mut.bed", sep = "\t", row.names = F, quote = F, col.names = F)
  rm(mut, pos_mut, mut_bed); gc()
}

# 编码区 ---------------------------------------------------------------------

get_cds_region("../predict_prob/Homo_sapiens.GRCh37.75.gtf")
gc()
# merge，以防重叠区间
system("sort -k 1,1 -k 2,2n cds_region.bed > sorted_cds_region.bed")
#Sys.getenv("PATH")
system("bedtools merge -i sorted_cds_region.bed > merged_sorted_cds_region.bed")
file.remove("sorted_cds_region.bed")

# 非编码区 --------------------------------------------------------------------
system("sort -k 1,1 -k 2,2n human.hg19.genome > sorted_human.hg19.genome")  #bedtools包中附带文件
system("bedtools complement -i merged_sorted_cds_region.bed -g sorted_human.hg19.genome > region_noncoding.bed")

# Promoter 区 --------------------------------------------------------------
# 选取-2500 到 +500 作为启动子区域
up_site = 2500
down_site = 500
##过滤数据
gtf_ensmble <- fread("../predict_prob/Homo_sapiens.GRCh37.75.gtf", 
                     sep = "\t", colClasses=list(character= 1))
names(gtf_ensmble)[1] <- "V1"
ensmble_pro_cod <- gtf_ensmble[gtf_ensmble$V2 == "protein_coding",] #选择仅仅是蛋白质编码的转录本
transcript_pro_cod_ensmble <- ensmble_pro_cod[ensmble_pro_cod$V3 == "transcript",]          ##########从转录本中选取总的转录本区间
options(scipen=30)
chr_index <- c(c(1:22),"X","Y")
chr_trans_pro <- transcript_pro_cod_ensmble[transcript_pro_cod_ensmble$V1 %in% chr_index, ] #仅仅挑选chr1:22 + X + Y 染色体
##从SST获取up_site到down_site的区间
pos_trans_chr_gtf <- chr_trans_pro[chr_trans_pro$V7 == "+",]
neg_trans_chr_gtf <- chr_trans_pro[chr_trans_pro$V7 == "-",]
pos_chr <- paste("chr",pos_trans_chr_gtf$V1, sep = "")
pos_start <- pos_trans_chr_gtf$V4 - up_site
pos_end <- pos_trans_chr_gtf$V4 + down_site
pos_bed <- data.frame(pos_chr, pos_start, pos_end)
names(pos_bed) <- c("chr","start","end")

neg_chr <- paste("chr",neg_trans_chr_gtf$V1, sep = "")
neg_start <- neg_trans_chr_gtf$V5 + up_site
neg_end <- neg_trans_chr_gtf$V5 - down_site
neg_bed <- data.frame(neg_chr, neg_end, neg_start)
names(neg_bed) <- c("chr","start","end")

promoter_region <- rbind(pos_bed, neg_bed)
promoter_region$index <- 1
write.table(promoter_region,"promoter_region.bed", sep = "\t", row.names = F, col.names = F,quote = F)

# 同样地，排序并merge
system("sort -k 1,1 -k 2,2n promoter_region.bed > sorted_promoter_region.bed")
system("bedtools merge -i sorted_promoter_region.bed > merged_sorted_promoter_region.bed")
file.remove("sorted_promoter_region.bed")

# 非Promoter的noncoding区 ----------------------------------------------------
# noncoding区间除了Promoter之外的区间

system("bedtools subtract -a region_noncoding.bed -b merged_sorted_promoter_region.bed > region_non_promoter.bed")


rm(list = ls()); gc()
source("functions.R")
# 获取不同区间上的突变数目 --------------------------------------------------------
region.coding = "merged_sorted_cds_region.bed"
region.noncoding = "region_noncoding.bed"
region.promoter = "merged_sorted_promoter_region.bed"
region.nonpromoter = "region_non_promoter.bed"

original_bed = "single_mut.bed"

get_reg_mution(original_bed, region.coding, "mut_coding.bed")
get_reg_mution(original_bed, region.noncoding, "mut_noncoding.bed")
get_reg_mution(original_bed, region.promoter, "mut_promoter.bed")
get_reg_mution(original_bed, region.nonpromoter, "mut_nonpromoter.bed")


# 分割基因组区间 -----------------------------------------------------------------
get_region_bed(1000000L)
system("sort -k 1,1 -k 2,2n hg19_mb_1M.bed > sorted_hg19_1M.bed")

# 计算突变数 -------------------------------------------------------------------
get_reg_mution_number("mut_coding.bed", "sorted_hg19_1M.bed", "1M_mut_coding.tsv")
get_reg_mution_number("mut_noncoding.bed", "sorted_hg19_1M.bed", "1M_mut_noncoding.tsv")
get_reg_mution_number("mut_promoter.bed", "sorted_hg19_1M.bed", "1M_mut_promoter.tsv")
get_reg_mution_number("mut_nonpromoter.bed", "sorted_hg19_1M.bed", "1M_mut_nonpromoter.tsv")


# 计算突变的各种注释结果 -------------------------------------------------------------

# 这里取自师兄之前的结果
# 下载原数据
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/mutation_counts_mb.tsv",
              destfile = "mutation_counts_mb.tsv")
# CpG
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/CpG_percent_Mb.bed",
              destfile = "CpG_percent_Mb.bed")
# GC
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/hg19_gcMb.bed",
              destfile = "hg19_gcMb.bed")

# Conservation
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/conservation_Mb.bed",
              destfile = "conservation_Mb.bed")

# Reptime
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/reptime_Mb.bed",
              destfile = "reptime_Mb.bed")

# polII
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/polII_mb.bed",
              destfile = "polII_mb.bed")

# mappability
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/mean_mb_mappability.tsv",
              destfile = "mean_mb_mappability.tsv")

# rec_rate
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/rec_rate.tsv",
              destfile = "rec_rate.tsv")

# tfbs
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/tfbs_mb.bed",
              destfile = "tfbs_mb.bed")



# 收集遗传变异相关性 ---------------------------------------------------------------
rm(list = ls());gc()

mut_coding      = "1M_mut_coding.tsv"
#mut_noncoding  = "1M_mut_noncoding.tsv"
mut_noncoding   = "mutation_counts_mb.tsv" #为保证结果一致，这里使用与师兄一样的文件
mut_promoter    = "1M_mut_promoter.tsv"
mut_nonpromoter = "1M_mut_nonpromoter.tsv"


annot_CpG          = "CpG_percent_Mb.bed"
annot_GC           = "hg19_gcMb.bed"
annot_conservation = "conservation_Mb.bed"
annot_Reptime      = "reptime_Mb.bed"
annot_polII        = "polII_mb.bed"
annot_mappability  = "mean_mb_mappability.tsv"
annot_rec_rate     = "rec_rate.tsv"
annot_tfbs         = "tfbs_mb.bed"

annot = c(annot_CpG, annot_GC, annot_conservation, annot_Reptime,
          annot_polII, annot_mappability, annot_rec_rate, annot_tfbs)

# join by x with V4
# V4, V4, V1, V4, V1, V1, V4, V1

obtain_cor = function(annot_list, join_cols, mut_df) {
  require(data.table)
  out = data.table()
  for (i in seq_along(annot_list)) {
    message("Runing annotation #", i, ": ", annot_list[i])
    value = fread(annot_list[i])
    mut = fread(mut_df)
    result = merge(value, mut, by.x = colnames(value)[join_cols[i]], by.y = "V4")

    # calculate correlation
    if (i == 5) {
      result = result[V3.x != 0]
      tmp = data.table(
        coeff = cor(result$V3.x, result$V5),
        p_val = cor.test(result$V3.x, result$V5)$p.value
      )
    } else if (i == 7) {
      result = result[avg != 0]
      tmp = data.table(
        coeff = cor(result$avg, result$V5.y),
        p_val = cor.test(result$avg, result$V5.y)$p.value
      )
    } else {
      result = result[(V5.x != 0) & (V5.x != ".")]
      result[, V5.x:=as.numeric(V5.x)]
      result[, V5.y:=as.numeric(V5.y)]
      tmp = data.table(
        coeff = cor(result$V5.x, result$V5.y),
        p_val = cor.test(result$V5.x, result$V5.y)$p.value
      )
    }

    out = rbind(out, tmp)
  }
  
  out
}

cor_genetic = list()
cor_genetic[["noncoding"]] = obtain_cor(annot, 
                                        join_cols = c(4, 4, 1, 4, 1, 1, 4, 1), 
                                        mut_noncoding)

cor_genetic[["coding"]] = obtain_cor(annot, 
                                        join_cols = c(4, 4, 1, 4, 1, 1, 4, 1), 
                                        mut_coding)
cor_genetic[["promoter"]] = obtain_cor(annot, 
                                     join_cols = c(4, 4, 1, 4, 1, 1, 4, 1), 
                                     mut_promoter)
cor_genetic[["nonpromoter"]] = obtain_cor(annot, 
                                       join_cols = c(4, 4, 1, 4, 1, 1, 4, 1), 
                                       mut_nonpromoter)


cor_genetic[["noncoding"]]
cor_genetic[["coding"]]
cor_genetic[["promoter"]]
cor_genetic[["nonpromoter"]]

cor_genetic_df = rbindlist(cor_genetic, idcol = TRUE)
setnames(cor_genetic_df, ".id", "region")
cor_genetic_df[, features:=rep(c("CpG island", "GC content",
                                "Conservation", "Replication time",
                                "DNA poly II", "Mappability",
                                "Recombination rate", "TFBS"), 4)]
cor_genetic_df$features = factor(cor_genetic_df$features,
                                 c("Mappability", "Replication time",
                                   "TFBS", "GC content",
                                   "CpG island", "DNA poly II",
                                   "Conservation", 
                                   "Recombination rate"))
library(ggplot2)
library(cowplot)
library(RColorBrewer)

ggplot(cor_genetic_df, aes(x = features, y = coeff, fill=region)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Paired", 
                    labels = c("Coding", "Noncoding", 
                               "Noncoding-promoter", "Promoter")) +
  labs(x = "Genetic features", y = "Correlation coefficient", fill = "Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> p_genetic

save_plot("Genetic_corrplot.pdf", plot = p_genetic, base_aspect_ratio = 1.6)  

cor_genetic_df2 = cor_genetic_df
setnames(cor_genetic_df2,
         colnames(cor_genetic_df2),
         c("Genomic region", "Correlation coefficient",
           "P value", "Genetic feature"))
cor_genetic_df2[, `Genomic region` := 
                  ifelse(`Genomic region` == "nonpromoter",
                         "Noncoding-promoter", `Genomic region`)]

library(gt)
gt_genetic <- gt(data = cor_genetic_df2)
gt_genetic 

# value <- fread("CpG_percent_Mb.bed")
# mut <- fread("mutation_counts_mb.tsv")
# result <- merge(value, mut, by.x = "V4", by.y = "V4")
# result <- result[result$V5.x != 0,]
# cor(result$V5.x, result$V5.y)
# cor.test(result$V5.x, result$V5.y)
# 
# 
# mut_num <- read.table("mutation_counts_mb.tsv", stringsAsFactors = F)
# mut_num <- mut_num[mut_num$V5 != ".",]
# mut_num$V5 <- as.numeric(mut_num$V5)
# 
# reptime <- read.table("reptime_Mb.bed", stringsAsFactors = F)
# reptime <- reptime[reptime$V5 != ".",]
# reptime$V5 <- as.numeric(reptime$V5)
# result <- merge(reptime, mut_num, by.x = "V4", by.y = "V4")
# 
# result <- result[result$V5.x != 0, ]
# cor(result$V5.x, result$V5.y)
# cor.test(result$V5.x, result$V5.y)
