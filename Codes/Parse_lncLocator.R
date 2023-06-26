#Import required library
library(caret)
library(stringr)
library(seqinr)
library(data.table)

get_cleaned_output <- function(lncLocator_file){
  print(head(lncLocator_file))
  lncLocator_preds = lncLocator_file[seq(from=1, to=nrow(lncLocator_file), by=7),]
  lncLocator_preds$label = ifelse(lncLocator_preds$location=='Exosome', 'EV', 'IC')
  head(lncLocator_preds)
  dim(lncLocator_preds)
  return(lncLocator_preds)
}


draw_confusion_matrix <- function(cm, Class1, Class2, title_val) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title(title_val, cex.main=2)
  
  # create the matrix 
  rect(150, 430, 240, 370, col='#3F97D0')
  text(195, 435, Class1, cex=1.2)
  rect(250, 430, 340, 370, col='#F7AD50')
  text(295, 435, Class2, cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='#F7AD50')
  rect(250, 305, 340, 365, col='#3F97D0')
  text(140, 400, Class1, cex=1.2, srt=90)
  text(140, 335, Class2, cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  



#####################################################################
################### Parsing out EV sequences  ###################
#####################################################################
lncLocator_files_list = list.files('Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/EV/lncLocator/', full.names = T)
lncLocator_files_names = gsub('.csv', '',list.files('Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/EV/lncLocator/'))

lncLocator_files = lapply(lncLocator_files_list[1:3], read.table, header = T)
lncLocator_files[[4]] = read.csv(lncLocator_files_list[4])

names(lncLocator_files) = lncLocator_files_names

lapply(lncLocator_files, head)
lapply(lncLocator_files, dim)

lncLocator_preds_list <- lapply(lncLocator_files, get_cleaned_output)
names(lncLocator_preds_list) = lncLocator_files_names

lapply(lncLocator_preds_list, head)
lapply(lncLocator_preds_list, dim)
lncLocator_preds_all <- do.call(rbind, lncLocator_preds_list)

write.csv(lncLocator_preds_all, 
          file = 'Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/EV/lncLocator/test_EV_subsample_lncLocator_all.csv',
          row.names = T, quote = F)



#####################################################################
################### Parsing out IC sequences  ###################
#####################################################################

lncLocator_files_list = list.files('Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/IC/lncLocator/', full.names = T)
lncLocator_files_names = gsub('.csv', '',list.files('Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/IC/lncLocator/'))

lncLocator_files = lapply(lncLocator_files_list, read.csv, header = T)
lapply(lncLocator_files, head)

names(lncLocator_files) = lncLocator_files_names

lapply(lncLocator_files, head)
lapply(lncLocator_files, dim)

lncLocator_preds_list <- lapply(lncLocator_files, get_cleaned_output)
names(lncLocator_preds_list) = lncLocator_files_names

lapply(lncLocator_preds_list, head)
lapply(lncLocator_preds_list, dim)
lncLocator_preds_all <- do.call(rbind, lncLocator_preds_list)
dim(lncLocator_preds_all)
head(lncLocator_preds_all)

write.csv(lncLocator_preds_all, 
          file = 'Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/IC/lncLocator/test_IC_subsample_lncLocator_all.csv',
          row.names = T, quote = F)

#####################################################################


lncLocator_preds_IC <- read.csv('Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/IC/lncLocator/test_IC_subsample_lncLocator_all.csv')
lncLocator_preds_EV <- read.csv('Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/EV/lncLocator/test_EV_subsample_lncLocator_all.csv')

table(lncLocator_preds_IC$label)
table(lncLocator_preds_EV$label)

lncLocator_preds_IC$true_label = 'IC'
lncLocator_preds_EV$true_label = 'EV'
lncLocator_preds_total <- rbind(lncLocator_preds_EV, lncLocator_preds_IC)
colnames(lncLocator_preds_total) <- c('file_seq_id', 'seq_id', 'pred_location', 'pred_label', 'true_label')
head(lncLocator_preds_total)
dim(lncLocator_preds_total)


write.csv(lncLocator_preds_total, 
          file = 'Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/lncLocator_predictions_total.csv',
          row.names = T, quote = F)




################################################
######## Creating confusion matrix
lncLocator_preds_total <- read.csv('Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/lncLocator_predictions_total.csv')

cm <- caret::confusionMatrix(data=factor(lncLocator_preds_total$pred_label,levels=c('IC', 'EV')), 
                             reference = factor(lncLocator_preds_total$true_label, levels=c('IC', 'EV')))
cm
draw_confusion_matrix(cm, Class1='IC', Class2='EV', title_val = 'lncLocator')




