# 获取区间突变
get_reg_mution = function(mutation_file, region_file, result_file)
{
  temp1 <- paste("sorted_", mutation_file, sep = "")
  temp2 <- paste("sorted_", region_file, sep = "")
  
  cmd1 <- paste("sort -k 1,1 -k 2,2n ",mutation_file," > ",temp1, sep = "")
  cmd2 <- paste("sort -k 1,1 -k 2,2n ",region_file," > ",temp2, sep = "")
  system(cmd1)
  system(cmd2)
  cmd3 <- paste("bedtools intersect -a ",temp1," -b ",temp2," -sorted > ",result_file, sep = "")
  system(cmd3)
  file.remove(temp1, temp2)
}


#获取CDS的外显子区域
get_cds_region <- function(gtf)
{
  gtf <- fread(gtf, skip = 5)
  
  chr_index <- c(c(1:22),"X","Y")
  gtf <- gtf[gtf$V1 %in% chr_index,]
  gtf_exon <- gtf[gtf$V3 == "exon",]
  
  cds <- gtf_exon[,c(1,4,5)]
  cds$V1 <- paste("chr",as.character(cds$V1),sep = "")
  write.table(cds,"cds_region.bed", sep = "\t", row.names = F, col.names = F, quote = F)
}
