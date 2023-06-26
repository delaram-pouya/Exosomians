library(stringr)

draw_confusion_matrix <- function(cm, Class1, Class2) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)
  
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



get_cleaned_output <- function(predictions_text){
  string_pattern = 'location is'
  pred_lines = predictions_text[grep(string_pattern, predictions_text)]
  
  pred_lines = gsub("The sequence ","",pred_lines)
  pred_lines = gsub(" in the query sequence "," ",pred_lines)
  pred_lines = gsub('. Its subcellular location is',"",pred_lines)
  pred_lines = gsub(". The probability score is ","",pred_lines)
  pred_lines = gsub("Cytoplasm, Cytosol","Cytoplasm,Cytosol",pred_lines)
  
  
  line_split_list <- str_split(pred_lines, ' ')
  pred_df = data.frame(seq_ids = unlist(lapply(line_split_list, '[[', 2)),
                       seq_num=unlist(lapply(line_split_list, '[[', 4)),
                       loc_pred=unlist(lapply(line_split_list, '[[', 5)),
                       loc_prob=unlist(lapply(line_split_list, '[[', 6)))
  
  pred_df$loc_prob = as.numeric(substr(pred_df$loc_prob,1,nchar(pred_df$loc_prob)-1))
  pred_df$seq_ids = gsub('>','',pred_df$seq_ids)
  
  table(pred_df$loc_pred)
  
  return(pred_df)
}


#####################################################################
################### Parsing out EV sequences  ###################
#####################################################################

iLocLnc_files_list = list.files('Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/EV/iLoc-LncRNA/', full.names = T)
iLocLnc_files_names = gsub('.txt', '',list.files('Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/EV/iLoc-LncRNA/'))

iLocLnc_files = lapply(iLocLnc_files_list, readLines)
predictions_text = iLocLnc_files[[4]]


iLocLnc_preds_list <- lapply(iLocLnc_files, get_cleaned_output)
names(iLocLnc_preds_list) <- iLocLnc_files_names
lapply(iLocLnc_preds_list, head)
lapply(iLocLnc_preds_list, function(x) table(x$loc_pred))


iLocLnc_preds_all <- do.call(rbind, iLocLnc_preds_list)
dim(iLocLnc_preds_all)
head(iLocLnc_preds_all)
table(iLocLnc_preds_all$loc_pred)

write.csv(iLocLnc_preds_all, 
          file = 'Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/EV/iLoc-LncRNA/test_EV_subsample_iLocLncRNA_all.csv',
          row.names = T, quote = T)


#####################################################################
################### Parsing out IC sequences  ###################
#####################################################################

### toDO: fix the bug

iLocLnc_files_list = list.files('Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/IC/iLoc-LncRNA/', full.names = T)
iLocLnc_files_names = gsub('.txt', '',list.files('Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/IC/iLoc-LncRNA/'))

iLocLnc_files = lapply(iLocLnc_files_list, readLines)
iLocLnc_preds_list <- lapply(iLocLnc_files, get_cleaned_output)
names(iLocLnc_preds_list) <- iLocLnc_files_names
lapply(iLocLnc_preds_list, head)
lapply(iLocLnc_preds_list, function(x) table(x$loc_pred))


iLocLnc_preds_all <- do.call(rbind, iLocLnc_preds_list)
dim(iLocLnc_preds_all)
head(iLocLnc_preds_all)
table(iLocLnc_preds_all$loc_pred)

write.csv(iLocLnc_preds_all, 
          file = 'Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/IC/iLoc-LncRNA/test_IC_subsample_iLocLncRNA_all.csv',
          row.names = T, quote = T)


#####################################################################


iLocLnc_preds_IC <- read.csv('Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/IC/iLoc-LncRNA/test_IC_subsample_iLocLncRNA_all.csv')
iLocLnc_preds_EV <- read.csv('Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/EV/iLoc-LncRNA/test_EV_subsample_iLocLncRNA_all.csv')

table(iLocLnc_preds_IC$loc_pred)
table(iLocLnc_preds_EV$loc_pred)


iLocLnc_preds_IC$pred_label = ifelse(iLocLnc_preds_IC$loc_pred=='Exosome', 'EV', 'IC')
iLocLnc_preds_EV$pred_label = ifelse(iLocLnc_preds_EV$loc_pred=='Exosome', 'EV', 'IC')


iLocLnc_preds_IC$true_label = 'IC'
iLocLnc_preds_EV$true_label = 'EV'

iLocLnc_preds_total <- rbind(iLocLnc_preds_EV, iLocLnc_preds_IC)
head(iLocLnc_preds_total)
colnames(iLocLnc_preds_total)
colnames(iLocLnc_preds_total) <- c('file_seq_id', 'seq_id','seq_num' ,'pred_location','pred_prob' ,'pred_label', 'true_label')
head(iLocLnc_preds_total)
dim(iLocLnc_preds_total)


write.csv(iLocLnc_preds_total, 
          file = 'Data/exoGRU_seqs/split_test_seqs/local_pred_outputs/iLocLnc_predictions_total.csv',
          row.names = T, quote = F)




#Import required library
library(caret)
#Creating confusion matrix
cm <- caret::confusionMatrix(data=factor(iLocLnc_preds_total$pred_label,levels=c('IC', 'EV')), 
                      reference = factor(iLocLnc_preds_total$true_label, levels=c('IC', 'EV')))
cm
draw_confusion_matrix(cm, Class1='IC', Class2='EV')

