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
# 外显子移除高变免疫区

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

#########获取不同区间大小的保守性数值########
~/bigWigAverageOverBed hg19.100way.phastCons.bw sort_hg19_mb.bed conservation_Mb.bed

##########计算每个区间内CpG岛的长度#########
bedtools map -a sort_hg19_mb.bed -b  sort_cpg_length.bed  -c 4 -o sum >  CpG_length_sum_Mb.bed

files <- c("CpG_length_sum_Mb.bed")
read_file <- vector("list", length = 1)
names(read_file) <- files
for(i in files) read_file[[i]] <- fread(i)

for(i in files)
{
  df <- read_file[[i]]
  df <- df[df$V5 != ".",]
  df$V5 <- as.numeric(df$V5)
  read_file[[i]] <- df
}

read_file$CpG_length_sum_Mb.bed$V5 <- read_file$CpG_length_sum_Mb.bed$V5/1000000

write.table(read_file$CpG_length_sum_Mb.bed,"CpG_percent_Mb.bed",
            sep = "\t", row.names = F, col.names = F, quote = F)



###########GC含量###########################
{
  library(data.table)
  library(stringr)
  options(scipen = 30)
  
  gc <- fread("hg19.gc5Base.txt", nThread = 8, sep = ",", header = F)
  chr_index <- grep("variableStep chrom=", gc$V1)
  chr_index_n <- chr_index[1:25]
  
  gc_chr <- vector("list", length = 24)
  chr_names <- gc$V1[chr_index_n][1:24]
  names(gc_chr) <- chr_names
  for(i in 1:24)   gc_chr[[chr_names[i]]] <- gc[chr_index_n[i]:chr_index_n[i+1],]
  
  for(i in chr_names)
  {
    df <- gc_chr[[i]]
    dim_length <- dim(df)[1]
    df_value <- df$V1[2:(dim_length-1)]
    pos <- str_split_fixed(df_value,"\t",2)[,1]
    gc_value <- str_split_fixed(df_value,"\t",2)[,2]
    pos_start <- as.numeric(pos) - 1
    pos_end <- as.numeric(pos) + 4
    chr <- strsplit(i,"variableStep chrom=| span=5")[[1]][2]
    bed_df <- data.frame(chr, pos_start, pos_end, gc_value)
    bed_list[[i]] <- bed_df
    fwrite(bed_df, paste(chr,"_gc5base.bed", sep = ""), sep = "\t", row.names = F, quote = F, col.names = F  )
  }
  
  
  #对 Mb 区间计算GC平均值（开多个终端，同时运算）
  bedtools map -a sort_hg19_mb.bed -b  chr1_gc5base.bed  -c 4 -o mean >  chr1_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr2_gc5base.bed  -c 4 -o mean >  chr2_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr3_gc5base.bed  -c 4 -o mean >  chr3_gcMb.bed
  
  bedtools map -a sort_hg19_mb.bed -b  chr4_gc5base.bed  -c 4 -o mean >  chr4_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr5_gc5base.bed  -c 4 -o mean >  chr5_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr6_gc5base.bed  -c 4 -o mean >  chr6_gcMb.bed
  
  bedtools map -a sort_hg19_mb.bed -b  chr7_gc5base.bed  -c 4 -o mean >  chr7_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr8_gc5base.bed  -c 4 -o mean >  chr8_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr9_gc5base.bed  -c 4 -o mean >  chr9_gcMb.bed
  
  bedtools map -a sort_hg19_mb.bed -b  chr10_gc5base.bed  -c 4 -o mean >  chr10_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr11_gc5base.bed  -c 4 -o mean >  chr11_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr12_gc5base.bed  -c 4 -o mean >  chr12_gcMb.bed
  
  bedtools map -a sort_hg19_mb.bed -b  chr13_gc5base.bed  -c 4 -o mean >  chr13_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr14_gc5base.bed  -c 4 -o mean >  chr14_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr15_gc5base.bed  -c 4 -o mean >  chr15_gcMb.bed
  
  
  bedtools map -a sort_hg19_mb.bed -b  chr16_gc5base.bed  -c 4 -o mean >  chr16_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr17_gc5base.bed  -c 4 -o mean >  chr17_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr18_gc5base.bed  -c 4 -o mean >  chr18_gcMb.bed
  
  bedtools map -a sort_hg19_mb.bed -b  chr19_gc5base.bed  -c 4 -o mean >  chr19_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr20_gc5base.bed  -c 4 -o mean >  chr20_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr21_gc5base.bed  -c 4 -o mean >  chr21_gcMb.bed
  
  bedtools map -a sort_hg19_mb.bed -b  chr22_gc5base.bed  -c 4 -o mean >  chr22_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chrX_gc5base.bed  -c 4 -o mean >  chrX_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chrY_gc5base.bed  -c 4 -o mean >  chrY_gcMb.bed
  
  #合并每条染色体结果（将上述不同长度区间的文件分别放入不同文件夹）
  files <- list.files()
  read_file <- vector("list", length = 24)
  names(read_file) <- files
  for(i in files)  read_file[[i]] <- fread(i)
  
  for(i in files)
  {
    chr <- strsplit(i, "_gc1kb.bed")[[1]][1]
    df <- read_file[[i]]
    read_file[[i]] <- df[df$V1 == chr,]
  }
  
  names(read_file) <- NULL
  result <- do.call("rbind", read_file)
  
  ###分别输出（）
  
}


##########TFBS##############################
{
  ###其中wgEncodeRegTfbsClusteredWithCellsV3.bed文件为
  ###UCSC下载的全基因组转录因子结合位点数据，论文中有链接地址
  tfbs <- fread("wgEncodeRegTfbsClusteredWithCellsV3.bed")
  bed <- tfbs[,c(1:3,5)]
  write.table(bed,"tfbs.bed", sep = "\t", col.names = F, row.names = F, quote = F)
  
  
  ###human_genome_chromosome_length.tsv文件为从human.hg19.genome文件中获取的1:22+X+Y染色体的长度数据
  sort -k1,1 -k2,2n tfbs.bed > sort_tfbs.bed
  sort -k1,1 -k2,2n human_genome_chromosome_length.tsv > hg19_chr_size.bed
  ###不区分转录因子，计算平均值，同时计算重叠部分的平均值
  bedtools merge -i sort_tfbs.bed -c 4 -o mean > mean_tfbs.bed
  ###将bed格式转换为bigwig格式
  ./bedGraphToBigWig mean_tfbs.bed hg19_chr_size.bed tfbs.bigwig
  
  ###分别计算1Mb区间长度上的转录因子结合位点数值
  ./bigWigAverageOverBed tfbs.bigwig sort_hg19_mb.bed tfbs_mb.bed
  
}


##########复制时间##########################
{
  bedtools map -a sort_hg19_mb.bed -b sort_average_reptime_14celllines.bed  -c 4 -o mean >  reptime_Mb.bed
}


##########POLII############################
{
  ./bigWigAverageOverBed  wgEncodeSydhTfbsK562Pol2s2StdSig.bigWig  sort_hg19_mb.bed  k562_pol2_Mb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsMcf10aesPol2Etoh01StdSig.bigWig  sort_hg19_mb.bed Mcf10aes_Etoh01_pol2_Mb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsMcf10aesPol2TamStdSig.bigWig  sort_hg19_mb.bed Mcf10aes_Tam_pol2_Mb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsPbdePol2UcdSig.bigWig  sort_hg19_mb.bed Pbde_pol2_Mb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsRajiPol2UcdSig.bigWig  sort_hg19_mb.bed Raji_pol2_Mb.bed
  
  setwd("/home/zhangjing/paper/icgc_download/polII/10kb")
  files <- list.files()
  read_file <- vector("list", length = 5)
  names(read_file) <- files
  for(i in files)  read_file[[i]] <- fread(i)
  mean_five <- do.call("cbind", read_file)[,c(5,11,17,23,29)] 
  mean_five <- as.matrix(mean_five)
  mean_value <- apply(mean_five, 1, mean)
  mean_result <- data.frame("num" = read_file[[1]][,1], "length" = read_file[[1]][,2], mean_value)
  write.table(mean_result, "polII_10kb.bed", sep = "\t", row.names = F, quote = F, col.names = F)
  
  
  
  setwd("/home/zhangjing/paper/icgc_download/polII/100kb")
  pol <- fread("polII_100kb.bed")
  mut <- fread("mutation_counts_100kb.tsv")
  result <- merge(mut,pol, by.x = "V4", by.y = "V1")
  result <- result[result$V3.y != 0, ]
  cor(result$V5, result$V3.y)
  
}
