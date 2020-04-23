
library(stringr)

setwd("/home/zhangjing/roadmap/blood")
setwd("/home/zhangjing/roadmap/breast")
setwd("/home/zhangjing/roadmap/Esophagus")
setwd("/home/zhangjing/roadmap/Kidney")
setwd("/home/zhangjing/roadmap/liver")
setwd("/home/zhangjing/roadmap/lung")
setwd("/home/zhangjing/roadmap/Ovary")
setwd("/home/zhangjing/roadmap/Pancreas")



broadpeak2bed()
system("mkdir bed")
system("mv *.broadPeak.bed bed")
setwd("bed")


file_names <- list.files()
anno <- str_split_fixed(file_names, "[-|.]",4)[,2]
anno_file <- vector(length = length(file_names) + 1)
anno_file[1] <- "genetic.bed"
for(i in 1:length(file_names))  anno_file[i + 1] <- paste(anno_file[i], anno[i], sep = ".")


****mv genetic.bed to current direct


for(i in 1:length(file_names))    
{
  pos_anno(anno_file[i], file_names[i], anno_file[i+1])
}



result_file <- fread(anno_file[length(file_names) + 1])
names(result_file)[14: (13 + length(file_names))] <- anno
result_file[result_file == "."]  <- 0
write.table(result_file,"final_result.tsv", sep = "\t", row.names = F, quote = F)
system("mv genetic.bed ..")
file.remove(anno_file)







######################################################## 多个数据

###melanoma
setwd("/home/zhangjing/roadmap/E059_Melanoma")


broadpeak2bed()
system("mkdir bed")
system("mv *.broadPeak.bed bed")
setwd("bed")



file_names <- list.files()
anno <- str_split_fixed(file_names, ".broadPeak.bed",4)[,1]
anno_file <- vector(length = length(file_names) + 1)
anno_file[1] <- "genetic.bed"
for(i in 1:length(file_names))  anno_file[i + 1] <- paste(anno_file[i], anno[i], sep = ".")


for(i in 1:length(file_names))    
{
  pos_anno(anno_file[i], file_names[i], anno_file[i+1])
}



result_file <- fread(anno_file[length(file_names) + 1])
names(result_file)[14: (13 + length(file_names))] <- anno
result_file[result_file == "."]  <- 0
write.table(result_file,"final_result.tsv", sep = "\t", row.names = F, quote = F)
system("mv genetic.bed ..")
file.remove(anno_file)




setwd("/home/zhangjing/roadmap/E061_Melanoma")


broadpeak2bed()
system("mkdir bed")
system("mv *.broadPeak.bed bed")
setwd("bed")



file_names <- list.files()
anno <- str_split_fixed(file_names, ".broadPeak.bed",4)[,1]
anno_file <- vector(length = length(file_names) + 1)
anno_file[1] <- "genetic.bed"
for(i in 1:length(file_names))  anno_file[i + 1] <- paste(anno_file[i], anno[i], sep = ".")


for(i in 1:length(file_names))    
{
  pos_anno(anno_file[i], file_names[i], anno_file[i+1])
}



result_file <- fread(anno_file[length(file_names) + 1])
names(result_file)[14: (13 + length(file_names))] <- anno
result_file[result_file == "."]  <- 0
write.table(result_file,"final_result.tsv", sep = "\t", row.names = F, quote = F)
system("mv genetic.bed ..")
file.remove(anno_file)








install.packages("pROC")
install.packages("caret")


#十折交叉验证和ROC曲线
###肺癌
{
  source("/home/zhangjing/paper/lasso/logit_function20181112.R")
  setwd("/home/zhangjing/roadmap/lung/bed")
  system("mv final_result.tsv cal_model_lung_final_result.tsv")
  
  
  
  project <- c("LUSC-CN","LUSC-KR","mock")
  lung_mut <- logit_form(project, "cal_model_lung_final_result.tsv")
  
  
  ###多重共线性检验
  mt <- lung_mut[,c(8:11,14:19)] 
  x <- cor(mt)         
  kappa(x, exact = TRUE)         # 8.947881                   
  mt <- lung_mut[,c(8:11)]
  x <- cor(mt)
  kappa(x, exact = TRUE)         # 3.16644        
  mt <- lung_mut[,c(14:19)]
  x <- cor(mt)
  kappa(x, exact = TRUE)		   # 5.184799
  
  
  #检验数值类型
  is.character(lung_mut[,7])
  is.numeric(as.matrix(lung_mut[,c(8:11,14:19)]))
  
  
  #整理格式及结果计算
  lung_form <- lung_mut[,c(5,7:19)]
  cross_val(lung_form,"lung_logit.RData")                       # 0.6850658           
}                                  



##食道癌
{
  source("/home/zhangjing/paper/lasso/logit_function20181112.R")
  setwd("/home/zhangjing/roadmap/Esophagus/bed")
  system("mv final_result.tsv cal_model_esophagus_final_result.tsv")
  
  project <- c("ESAD-UK","ESCA-CN","mock")
  esophagus_mut <- logit_form(project, "cal_model_esophagus_final_result.tsv")
  
  
  ###多重共线性检验
  mt <- esophagus_mut[,c(8:11,14:19)] 
  x <- cor(mt)
  kappa(x, exact = TRUE)            # 8.529881           
  mt <- esophagus_mut[,c(8:11)]
  x <- cor(mt)
  kappa(x, exact = TRUE)            # 3.178326      
  mt <- esophagus_mut[,c(14:19)]
  x <- cor(mt)
  kappa(x, exact = TRUE)			  # 4.565268
  
  
  
  #检验数值类型
  is.character(esophagus_mut[,7])
  is.numeric(as.matrix(esophagus_mut[,c(8:11,14:19)]))
  
  #整理格式
  esophagus_form <- esophagus_mut[,c(5,7:19)]
  cross_val(esophagus_form,"esophagus_logit.RData")
}



####肝癌
{
  source("/home/zhangjing/paper/lasso/logit_function20181112.R")
  setwd("/home/zhangjing/roadmap/liver/bed")
  system("mv final_result.tsv cal_model_liver_final_result.tsv")
  
  
  project <- c("LIAD-FR", "LICA-CN", "LICA-FR", "LINC-JP", "LIRI-JP","mock")
  liver_mut <- logit_form(project, "cal_model_liver_final_result.tsv")
  
  
  ###多重共线性检验
  mult_validate(liver_mut)	 # 8.991961		
  
  #检验数据类型
  type_validate(liver_mut)   
  
  #整理格式
  liver_form <-format_trans(liver_mut)
  cross_val(liver_form,"liver_logit.RData")   
}



{#######乳腺癌
  source("/home/zhangjing/paper/lasso/logit_function20181112.R")
  setwd("/home/zhangjing/roadmap/breast/bed")
  system("mv final_result.tsv cal_model_breast_final_result.tsv")
  
  project <- c("BRCA-EU","BRCA-FR","BRCA-US","mock")
  breast_mut <- logit_form(project, "cal_model_breast_final_result.tsv")
  
  #添加dnase数据
  a <- fread("E028-DNase.hotspot.broad.bed")[,-4]
  write.table(a,"E028-DNase.hotspot.broad.bed", sep = "\t", row.names = F, quote = F, col.names = F)
  
  
  add_res <- add_dnase(breast_mut, "E028-DNase.hotspot.broad.bed")
  breast_mut <- add_res
  
  ###多重共线性检验
  mult_validate(breast_mut)
  
  #检验数据类型
  type_validate(breast_mut)
  
  
  #整理格式
  breast_form <-format_trans(breast_mut)
  cross_val(breast_form,"breast_logit.RData")
}



## 胰腺癌
{
  source("/home/zhangjing/paper/lasso/logit_function20181112.R")
  setwd("/home/zhangjing/roadmap/Pancreas/bed")
  system("mv final_result.tsv cal_model_pancreas_final_result.tsv")
  
  
  project <- c("PACA-AU","PACA-CA","PAEN-AU","PAEN-IT","mock")
  pancreas_mut <- logit_form(project, "cal_model_pancreas_final_result.tsv")
  
  #添加dnase数据
  a <- fread("E098-DNase.hotspot.broad.bed")[,-4]
  write.table(a,"E098-DNase.hotspot.broad.bed", sep = "\t", row.names = F, quote = F, col.names = F)
  
  add_res <- add_dnase(pancreas_mut, "E098-DNase.hotspot.broad.bed")
  pancreas_mut <- add_res
  
  ###多重共线性检验
  mult_validate(pancreas_mut)
  
  #检验数据类型
  type_validate(pancreas_mut)
  
  #整理格式
  pancreas_form <-format_trans(pancreas_mut)
  cross_val(pancreas_form,"pancreas_logit.RData")
}



###肾癌
{
  source("/home/zhangjing/paper/lasso/logit_function20181112.R")
  setwd("/home/zhangjing/roadmap/Kidney/bed")
  system("mv final_result.tsv cal_model_kidney_final_result.tsv")
  
  
  project <- c("RECA-EU","mock")
  kidney_mut <- logit_form(project, "cal_model_kidney_final_result.tsv")
  
  #添加dnase数据
  a <- fread("E086-DNase.hotspot.broad.bed")[,-4]
  write.table(a,"E086-DNase.hotspot.broad.bed", sep = "\t", row.names = F, quote = F, col.names = F)
  
  add_res <- add_dnase(kidney_mut, "E086-DNase.hotspot.broad.bed")
  kidney_mut <- add_res
  
  ###多重共线性检验
  mult_validate(kidney_mut)			
  
  #检验数据类型
  type_validate(kidney_mut)
  
  #整理格式
  kidney_form <- format_trans(kidney_mut)
  # cross_val(kidney_form,"kidney_logit.RData")    #AUC 0.626 不合格,优化
  
  ##调整    res ~ nbase + reptime + conser + cpg + H3K27me3 + H3K36me3 + H3K9me3  
  kidney_form1 <- kidney_form[,c(1,2,3,5,8,9,10,14)]     ###使用该结果
  cross_val(kidney_form1,"kidney1_logit.RData") 
}



{####血癌
  source("/home/zhangjing/paper/lasso/logit_function20181112.R")
  setwd("/home/zhangjing/roadmap/blood/bed")
  system("mv final_result.tsv cal_model_blood_final_result.tsv")
  
  
  project <- c("ALL-US","CLLE-ES","MALY-DE","NKTL-SG","mock")
  blood_mut <- logit_form(project, "cal_model_blood_final_result.tsv")
  
  
  ###多重共线性检验
  mult_validate(blood_mut)			
  
  #检验数据类型
  type_validate(blood_mut)
  
  #整理格式
  blood_form <-format_trans(blood_mut)
  cross_val(blood_form,"blood_logit.RData")  
}



###肝癌
{
  source("/home/zhangjing/paper/lasso/logit_function20181112.R")
  setwd("/home/zhangjing/roadmap/liver/bed")
  system("mv final_result.tsv cal_model_liver_final_result.tsv")
  
  
  project <- c("LIAD-FR", "LICA-CN", "LICA-FR", "LINC-JP", "LIRI-JP","mock")
  liver_mut <- logit_form(project, "cal_model_liver_final_result.tsv")
  
  
  ###多重共线性检验
  mult_validate(liver_mut)			
  
  #检验数据类型
  type_validate(liver_mut)
  
  #整理格式
  liver_form <-format_trans(liver_mut)
  cross_val(liver_form,"liver_logit.RData")  
  
}



###卵巢癌
{
  source("/home/zhangjing/paper/lasso/logit_function20181112.R")
  setwd("/home/zhangjing/roadmap/Ovary/bed")
  system("mv final_result.tsv cal_model_ovary_final_result.tsv")
  
  project <- c("OV-AU","mock")
  ovary_mut <- logit_form(project, "cal_model_ovary_final_result.tsv")
  
  
  a <- fread("E097-DNase.hotspot.broad.bed")[,-4]
  write.table(a,"E097-DNase.hotspot.broad.bed", sep = "\t", row.names = F, quote = F, col.names = F)
  
  add_res <- add_dnase(ovary_mut, "E097-DNase.hotspot.broad.bed")
  ovary_mut <- add_res
  
  
  ###多重共线性检验
  mult_validate(ovary_mut)			
  
  #检验数据类型
  type_validate(ovary_mut)
  
  #整理格式
  ovary_form <-format_trans(ovary_mut)
  #cross_val(ovary_form,"ovary_logit.RData")   
  ###更改
  ovary_form1 <- ovary_form[,-c(6:8,15)]
  cross_val(ovary_form1,"ovary1_logit.RData")  
  
  
  
  
  
  
  
}



####黑色素瘤
{
  ###获取E059和E061的平均值
  setwd("/home/zhangjing/roadmap/E059_Melanoma/bed")
  system("mv final_result.tsv E059_final_result.tsv")
  setwd("/home/zhangjing/roadmap/E061_Melanoma/bed")
  system("mv final_result.tsv E061_final_result.tsv")
  
  df59 <- fread("/home/zhangjing/roadmap/E059_Melanoma/bed/E059_final_result.tsv", data.table = F)
  df61 <- fread("/home/zhangjing/roadmap/E061_Melanoma/bed/E061_final_result.tsv", data.table = F)
  e059 <- as.matrix(df59[,14:19])
  e061 <- as.matrix(df61[,14:19])
  
  mean_value <- (e059 + e061)/2
  mean_value <- as.data.frame(mean_value)
  
  names(mean_value) <- c("H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3")
  res <- data.frame(df59[,1:13], mean_value)
  write.table(res, "final_result.tsv", sep = "\t", row.names = F, quote = F)
  
  
  
  source("/home/zhangjing/paper/lasso/logit_function20181112.R")
  setwd("/home/zhangjing/roadmap/E061_Melanoma/bed")
  system("mv final_result.tsv cal_model_mela_final_result.tsv")
  
  project <- c("MELA-AU","SKCA-BR","SKCM-US","mock")
  mela_mut <- logit_form(project, "cal_model_mela_final_result.tsv")
  
  #添加dnase数据
  a <- fread("E059-DNase.hotspot.broad.bed")[,-4]
  write.table(a,"E059-DNase.hotspot.broad.bed", sep = "\t", row.names = F, quote = F, col.names = F)
  
  add_res <- add_dnase(mela_mut, "E059-DNase.hotspot.broad.bed")
  mela_mut <- add_res
  
  ###多重共线性检验
  mult_validate(mela_mut)			
  
  #检验数据类型
  type_validate(mela_mut)
  
  #整理格式
  mela_form <-format_trans(mela_mut)
  cross_val(mela_form,"mela_logit.RData")   	
}






###取所有的突变并删除不符合条件的位点
setwd("/home/zhangjing/predict_prob")
{###选取 "预处理所有突变.tsv" 文件中的单位点突变
  mut <- fread("预处理所有突变.tsv", data.table = F)
  pos_mut <- mut[mut$mut == "single base substitution",]
  
  
  ###以bed格式输出所有单碱基突变
  mut_bed <- pos_mut[,c(3:5,2)]
  mut_bed$start <- mut_bed$start - 1
  write.table(mut_bed, "single_mut.bed", sep = "\t", row.names = F, quote = F, col.names = F)
}


###删除maf0.05
rm_maf(0.05, "single_mut.bed")



###获取免疫高变区，并扩增指定bp
{
  ###从Ensmbles gtf 文件 Homo_sapiens.GRCh37.75.gtf 中获取免疫高变区区间
  immu_id <- c('IG_C_gene',
               'IG_D_gene',
               'IG_J_gene',
               'IG_LV_gene',
               'IG_V_gene',
               'TR_C_gene',
               'TR_J_gene',
               'TR_V_gene',
               'TR_D_gene',
               'IG_pseudogene',
               'IG_C_pseudogene',
               'IG_J_pseudogene',
               'IG_V_pseudogene',
               'TR_V_pseudogene',
               'TR_J_pseudogene')
  
  gtf <- fread("Homo_sapiens.GRCh37.75.gtf", sep = "\t", colClasses=list(character= 1))
  names(gtf)[1] <- "V1"
  
  chr_index <- c(c(1:22),"X","Y")
  gtf <- gtf[gtf$V1 %in% chr_index,]
  
  df_immu_region <- gtf[gtf$V2 %in% immu_id,]
  bed_immu_region <- data.frame(df_immu_region$V1,df_immu_region$V4, df_immu_region$V5)
  names(bed_immu_region) <- c("chr","start","end")
  bed_immu_region$chr <- paste("chr",as.character(bed_immu_region$chr),sep = "")
  write.table(bed_immu_region,"immu_region_from_gtf.bed", sep = "\t", row.names = F, col.names = F,quote = F)
  
  ###将免疫高变区merge
  system("sort -k 1,1 -k 2,2n immu_region_from_gtf.bed > sort_chr_immu_region_from_gtf.bed")
  system("bedtools merge -i sort_chr_immu_region_from_gtf.bed > merged_sort_chr_immu_region_from_gtf.bed")
  
  
  ###将免疫高变区前后扩增ext KB
  ext = 50000
  merged <- fread("merged_sort_chr_immu_region_from_gtf.bed")
  merged$V2 <- merged$V2 - ext
  merged$V3 <- merged$V3 + ext
  write.table(merged,"extended_merged_sort_chr_immu_region_from_gtf.bed", sep = "\t", row.names = F ,col.names = F, quote = F)
  
  ###再次merge
  system("sort -k 1,1 -k 2,2n extended_merged_sort_chr_immu_region_from_gtf.bed > sorted_extended_merged_sort_chr_immu_region_from_gtf.bed")
  system("bedtools merge -i sorted_extended_merged_sort_chr_immu_region_from_gtf.bed > merged_sorted_extended_merged_sort_chr_immu_region_from_gtf.bed")
  
  file.remove("sort_chr_immu_region_from_gtf.bed",
              "merged_sort_chr_immu_region_from_gtf.bed",
              "immu_region_from_gtf.bed",
              "sorted_extended_merged_sort_chr_immu_region_from_gtf.bed",
              "extended_merged_sort_chr_immu_region_from_gtf.bed")
  
}				


###获取免疫和外显子两个区域的并集区间
{
  cds <- fread("cds_region.bed")
  immu <- fread("merged_sorted_extended_merged_sort_chr_immu_region_from_gtf.bed")
  result <- rbind(cds,immu)
  write.table(result, "rbind_three.bed", sep = "\t", row.names = F, col.names = F, quote = F)	
  
  system("sort -k 1,1 -k 2,2n rbind_three.bed > sorted_rbind_three.bed")
  system("bedtools merge -i sorted_rbind_three.bed > merged_sorted_rbind_three.bed")
  
  file.remove("rbind_three.bed","sorted_rbind_three.bed")
}


###获取基因组上除去免疫和外显子的剩下的非编码区间
{
  system("sort -k 1,1 -k 2,2n merged_sorted_rbind_three.bed > sorted_merged_sorted_rbind_three.bed")
  system("sort -k 1,1 -k 2,2n human.hg19.genome > sorted_human.hg19.genome")
  ######获取全基因组mask区间的补集######
  system("bedtools complement -i sorted_merged_sorted_rbind_three.bed -g sorted_human.hg19.genome > hg19_nonmask_region.bed")
  
  file.remove("sorted_merged_sorted_rbind_three.bed","sorted_human.hg19.genome")
}



###获取非编码区的所有体细胞突变
{
  options(scipen = 30)
  df <- fread("rmSNP_single_mut.bed")
  write.table(df,"rmSNP_single_mut.bed", sep = "\t", row.names = F, quote = F, col.names = F)
  
  region_mut("rmSNP_single_mut.bed","hg19_nonmask_region.bed","noncoding_mut")
}

















