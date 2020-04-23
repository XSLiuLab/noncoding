###����ʮ�۽�����֤���������߼��ع�
library(caret)
cross_val <- function(df, result_name)
{
	#�趨���������
	set.seed(100)
	
	names(df)[1] <- "Y"
	folds <- createFolds(y=df$Y,k=10)

	max=0
	num=0
	
	for(i in 1:10)
	{
		fold_test <- df[folds[[i]],]   
		fold_train <- df[-folds[[i]],]
		
		print("*******Number*******")
		#ģ�����
		fold_pre <- glm(Y ~., family=binomial(link='logit'), data=fold_train)
		fold_predict <- predict(fold_pre, type='response', newdata=fold_test)
		
		#������Լ���ȷ��
		fold_predict <- ifelse(fold_predict>0.5,1,0)
		fold_test$predict <- fold_predict
		fold_error = fold_test$predict - fold_test$Y
		fold_accuracy = (nrow(fold_test)-sum(abs(fold_error)))/nrow(fold_test)
		print(i)
		print("***test_df_accuracy***")
		print(fold_accuracy)
		
		#����ѵ������ȷ��
		fold_predict2 <- predict(fold_pre,type='response',newdata=fold_train)
		fold_predict2 =ifelse(fold_predict2>0.5,1,0)
		fold_train$predict = fold_predict2
		fold_error2 = fold_train$predict - fold_train$Y
		fold_accuracy2 = (nrow(fold_train)-sum(abs(fold_error2)))/nrow(fold_train) 
		print("***train_df_accuracy***")
		print(fold_accuracy2)
		
 		
		if(fold_accuracy>max)
		{
			max=fold_accuracy  
			num=i
		}	
	}
	
	###���Ч����õ�һ������
	print("***Best_result: ***")
	print(max)
	print(num)
	
	
	###ʹ����õĽ��
	testi <- df[folds[[num]],]
	traini <- df[-folds[[num]],]
	
	prei <- glm(Y ~., family=binomial(link='logit'), data=traini)
	predicti. <- predict.glm(prei,type='response',newdata=testi)
	
	predicti =ifelse(predicti. >0.5,1,0)
    testi$predict = predicti
    errori = testi$predict-testi$Y
    accuracyi = (nrow(testi)-sum(abs(errori)))/nrow(testi) 
	print("****Best_test_df_accuracy:**")
	print(accuracyi)
	
	
	###��ROC����
	library(pROC)
	true_value <- testi$Y
	plot.roc(true_value, predicti., legacy.axes=TRUE, grid=c(0.1, 0.2),
	print.auc=TRUE, max.auc.polygon=TRUE, print.thres= F, ci = TRUE) 

	###����ѵ��ģ��
	save(prei, file = result_name)
}


###�߼��ع�������������
var_scale <- function(df)
{
	df_names <- c("donor","chr","start","end","res","base","nbase","reptime","tfbs","conser","gc","pro","cpg")
	if(!all(names(df)[1:13] == df_names))  break("Error!")
	df$reptime <- scale(df$reptime, center = F)
	df$tfbs <- scale(df$tfbs, center = F)
	df$gc <- scale(df$gc, center = F)
	df$nbase <- as.character(df$nbase)
	return(df)
}

samp_neg_df <- function(df)
{
	df_pos <- df[df$res == 1, ]
	df_neg <- df[df$res == 0, ]
	pos_num <- dim(df_pos)[1]
	neg_num <- dim(df_neg)[1]
	neg_index <- sample(1:neg_num, pos_num)
	samp_neg <- df_neg[neg_index,]
	rbind(df_pos, samp_neg)
}

logit_form <- function(project,file_name, donor = "~/paper/lasso/id_pro_index.tsv")
{	
	###ƥ����Ŀ��Ϣ
	id_index <- fread(donor)
	mock <- data.frame("icgc_donor_id" = "mock", "project_code" ="mock")
	id_index <- rbind(id_index, mock)
	result <- fread(file_name)
	merge_result <- merge(result, id_index, by.x = "V4", by.y = "icgc_donor_id", all.x = T)
	
	###ѡȡ�ض�����
	mut <- merge_result[merge_result$project_code %in% project,]
	
	###���ݹ�һ��
	names(mut)[1:13] <- c("donor","chr","start","end","res","base","nbase","reptime","tfbs","conser","gc","pro","cpg")
	mut <- var_scale(mut)
	
	###ѡȡ�����Ķ�������
	mut <- samp_neg_df(mut)
	
	return(as.data.frame(mut))
}


pos_anno <- function(pos_file, anno_file, result_file)
{
    temp1 <- paste("sorted_", pos_file, sep = "")
	temp2 <- paste("sorted_", anno_file, sep = "")
	cmd1 <- paste("sort -k 1,1 -k 2,2n ",pos_file," > ",temp1, sep = "")
	cmd2 <- paste("sort -k 1,1 -k 2,2n ",anno_file," > ",temp2, sep = "")
	system(cmd1)
	system(cmd2)
	cmd3 <- paste("bedtools map -a ",temp1," -b ",temp2," -c 4 -o mean > ",result_file, sep = "")
	system(cmd3)
	file.remove(temp1, temp2)
}




###����DNase����
add_dnase <- function(df,file_names)
{
	df <- df[,-c(1,dim(df)[2])]
	write.table(df,"pos.bed", sep = "\t", row.names = F, quote = F, col.names = F)
	pos_anno("pos.bed", file_names, "raw_add.bed")
	
	raw_add <- fread("raw_add.bed", stringsAsFactors = F, data.table = F)
	col_n <- dim(raw_add)[2]
	new_value <- raw_add[,col_n]
	raw_add[,col_n][new_value == "."] <- 0
	names(raw_add)[1:dim(df)[2]] <- names(df)
	names(raw_add)[col_n] <- "dnase"
	raw_add$dnase <- as.numeric(raw_add$dnase)
	raw_add$nbase <- as.character(raw_add$nbase)
	
	return(raw_add)
}


###���ع����Լ���
mult_validate <- function(df)
{
	anno <- c("reptime","tfbs","conser","gc","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3","H3K9ac","dnase")  
	df <- df[,names(df) %in% anno]
	x <- cor(df)
	
	print("***multiple validata is***")
	print(kappa(x, exact = TRUE))
}

###���ͼ���
type_validate <- function(df)
{
	print("is the base type character:")
	if(is.character(df$nbase))
	{	
		print("TRUE")
	}else{   
		print("FALSE")
	}
	
	anno <- c("reptime","tfbs","conser","gc","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3","H3K9ac","dnase")  
	
	print("is other type is numeric:")
	if(is.numeric(as.matrix(df[,names(df) %in% anno])))
	{
		print("TRUE")
	}else{
		print("FALSE")
	}
}


###������ʽ����
format_trans <- function(df)
{
	anno <- c("res","nbase","pro","cpg","reptime","tfbs","conser","gc","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3","H3K9ac","dnase")  
	return(df[ ,names(df) %in% anno])
}




