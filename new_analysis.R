# New Analysis BEGIN ------------------------------------------------------
# 计算1Mbp区间内突变频率与遗传、表观注释值的相关性
# Coding, Noconding, Promoter ...
# Coding为外显子区域

setwd("~/New_Analysis/")
source("functions.R")
library(data.table)
Sys.setenv(PATH=paste0("/home/zhangjing/bedtools2/bin:", Sys.getenv("PATH")))

# 预处理数据 -------------------------------------------------------------------

original_mut = "../icgc_data/预处理所有突变.tsv"    # 经过预处理过的突变记录文件
original_bed = "../icgc_data/rmSNP_single_mut.bed"  # SNV并移除超过比例SNP位点的区域
hyperImmuRegion = "../icgc_data/merged_sorted_extended_merged_sort_chr_immu_region_from_gtf.bed"


# 编码区 ---------------------------------------------------------------------
# 外显子移除高变免疫区

# 获取编码区
get_cds_region("../predict_prob/Homo_sapiens.GRCh37.75.gtf")
gc()
# merge，以防重叠区间
system("sort -k 1,1 -k 2,2n cds_region.bed > sorted_cds_region.bed")
#Sys.getenv("PATH")
system("bedtools merge -i sorted_cds_region.bed > merged_sorted_cds_region.bed")
file.remove("sorted_cds_region.bed")

# 减去高免疫突变区
system(paste0("bedtools subtract -a merged_sorted_cds_region.bed -b ",
              hyperImmuRegion, " > clean_cds.bed"))

region.coding = "clean_cds.bed"

# 非编码区 --------------------------------------------------------------------
region.noncoding = "../icgc_data/hg19_map_nonmask_region.bed"

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

# 同样地，去掉重叠以及高免疫突变区
# merge，以防重叠区间
system("sort -k 1,1 -k 2,2n promoter_region.bed > sorted_promoter_region.bed")
system("bedtools merge -i sorted_promoter_region.bed > merged_sorted_promoter_region.bed")
file.remove("sorted_promoter_region.bed")

# 减去高免疫突变区
system(paste0("bedtools subtract -a merged_sorted_promoter_region.bed -b ",
              hyperImmuRegion, " > clean_promoter.bed"))

region.promoter = "clean_promoter.bed"

# 非Promoter的noncoding区 ----------------------------------------------------
# noncoding区间除了Promoter之外的区间

# 减去高免疫突变区
system(paste0("bedtools subtract -a ", region.noncoding, " -b ", 
               region.promoter, " > clean_non_promoter.bed"))
region.nonpromoter = "clean_non_promoter.bed"


# 获取突变 --------------------------------------------------------------------
region.coding = "clean_cds.bed"
region.noncoding = "../icgc_data/hg19_map_nonmask_region.bed"
region.promoter = "clean_promoter.bed"
region.nonpromoter = "clean_non_promoter.bed"

get_reg_mution(original_bed, region.coding, "mut_coding.bed")
get_reg_mution(original_bed, region.noncoding, "mut_noncoding.bed")
get_reg_mution(original_bed, region.promoter, "mut_promoter.bed")
get_reg_mution(original_bed, region.nonpromoter, "mut_nonpromoter.bed")



# 计算突变的各种注释结果 -------------------------------------------------------------

###获取 noncoding_mut_site.bed 文件的三碱基
single_site <- fread("noncoding_mut_site.bed")
single_site$V2 <- single_site$V2 - 1
single_site$V3 <- single_site$V3 + 1
write.table(single_site, "three_base.bed", sep = "\t", row.names = F, col.names = F, quote = F)

pos_fasta("three_base.bed","hg19.fasta","noncoding_three_base.bed")
file.remove("three_base.bed")


three_site <- fread("noncoding_three_base.bed")
three_site$V2 <- three_site$V2 + 1
three_site$V3 <- three_site$V3 - 1
three_site$V6 <- toupper(three_site$V6)

base_type <- fread("base.tsv")
merge_result <- merge(three_site, base_type, by.x = "V6", by.y = "V1")
result <- merge_result[,c(2:6,1,7)]
write.table(result, "noncoding_mut_site_3basetype.bed", sep = "\t", row.names = F, col.names = F, quote = F)


###获取 noncoding_mut_site_3basetype.bed 文件的复制时间
pos_anno("noncoding_mut_site_3basetype.bed","average_reptime_14celllines.bed","noncoding_mut_site_3basetype_reptime.bed")

rep_value <- fread("noncoding_mut_site_3basetype_reptime.bed")
rep_value <- rep_value[rep_value$V8 != ".",]
write.table(rep_value, "noncoding_mut_site_3basetype_reptime.bed", sep = "\t", row.names = F, quote = F, col.names = F)


###获取 noncoding_mut_site_3basetype_reptime.bed 文件的tfbs
tfbs <- fread("wgEncodeRegTfbsClusteredWithCellsV3.bed")
bed <- tfbs[,c(1:3,5)]
write.table(bed,"TFBS_value.bed", sep = "\t", col.names = F, row.names = F, quote = F)

pos_anno("noncoding_mut_site_3basetype_reptime.bed","TFBS_value.bed","noncoding_mut_site_3basetype_reptime_tfbs.bed")
tfbs_value <- fread("noncoding_mut_site_3basetype_reptime_tfbs.bed")
tfbs_value$V9[tfbs_value$V9 == "."] <- 0
write.table(tfbs_value, "noncoding_mut_site_3basetype_reptime_tfbs.bed", sep = "\t", row.names = F, quote = F, col.names = F)



###获取 noncoding_mut_site_3basetype_reptime_tfbs.bed 文件位点的保守性值
pos_anno("noncoding_mut_site_3basetype_reptime_tfbs.bed","hg19.100way.phastCons.bw","noncoding_mut_site_3basetype_reptime_tfbs_cons.bed")


###获取 noncoding_mut_site_3basetype_reptime_tfbs.bed 文件位点的GC含量（1kb）
system("sort -k 1,1 -k 2,2n hg19_gc1kb.bed > sorted_hg19_gc1kb.bed")
gc <- fread("sorted_hg19_gc1kb.bed")
gc <- gc[gc$V5!= ".",]
write.table(gc, "sorted_hg19_gc1kb.bed", sep = "\t", row.names = F, col.names = F, quote = F)
system("bedtools map -a noncoding_mut_site_3basetype_reptime_tfbs_cons.bed -b sorted_hg19_gc1kb.bed -c 5 -o mean > noncoding_mut_site_3basetype_reptime_tfbs_cons_gc.bed")

file.remove("sorted_hg19_gc1kb.bed")

#？？？
system("sort -k 1,1 -k 2,2n promoter_region.bed > sorted_promoter_region.bed")
system("bedtools map -a noncoding_mut_site_3basetype_reptime_tfbs_cons_gc.bed -b sorted_promoter_region.bed -c 4 -o mean > noncoding_mut_site_3basetype_reptime_tfbs_cons_gc_promoter.bed")

pro <- fread("noncoding_mut_site_3basetype_reptime_tfbs_cons_gc_promoter.bed")
pro$V12[pro$V12 == "."] <- "N"
pro$V12[pro$V12 == 1] <- "Y"
write.table(pro, "noncoding_mut_site_3basetype_reptime_tfbs_cons_gc_promoter.bed", sep = "\t", row.names = F, col.names = F, quote = F)


file.remove("promoter_region.bed","sorted_promoter_region.bed")

