# In this file, we will read the data,
#   Try to filter the noisy (false) regions (samples/sequences) based on their IC expression value,
#   Then label the data, i.e. which sample/sequence is exported via exosome

#   Input: Raw data files (smRNA_counts) 
#   Output: Labels Matrix


source('Codes/Functions.R')
Initialize()


#### Reading raw data from file ####

COUNTS_DATA_DIR = 'Data/smRNA_counts'
filesPath = list.files(COUNTS_DATA_DIR, full.names = T)
countsFiles = lapply(filesPath, read.delim, header=T)
names(countsFiles) <- substr(list.files(COUNTS_DATA_DIR),1,2)

countsFiles <- lapply(countsFiles, function(aFile) {
  rownames(aFile) = aFile$X
  aFile$X = NULL
  aFile
})



## SKBR==SKBR3  & CNLM1A==LM1A
### HCC1395 needs to be removed 

RefineCellLineNames <-  function(Files) {
  
  EVcellLines <- toupper(gsub('\\..*','',colnames(Files[['EV']])))
  ICcellLines <- toupper(gsub('\\..*','',gsub('\\_.*','',colnames(Files[['IC']]))))
  
  ICcellLines[ICcellLines=='SKBR'] <- 'SKBR3'
  EVcellLines[EVcellLines=='CNLM1A'] <- 'LM1A' 
  
  colnames(Files[['EV']]) <- EVcellLines
  colnames(Files[['CM']]) <- EVcellLines
  colnames(Files[['IC']]) <- ICcellLines
  
  Files <- lapply(Files, function(x) x[,!colnames(x) %in% c('HCC1395','HCC1395.1')] )
  Files
}

Counts <- RefineCellLineNames(countsFiles)

### filtering Counts based on Valid IDs defined by Ali
ValidIDs <- designMat2$id
Counts <- lapply(Counts, function(x) x[rownames(x) %in% ValidIDs,])



#### Filtering noisy (false) regions (samples) based on IC expression value ####

IC_NOISE_COUNTS_THRESHOLD <- quantile(rowSums(Counts[['IC']]), 0.75)
IC_FILTER_PASSED_INDICES <- rowSums(Counts[['IC']]) > IC_NOISE_COUNTS_THRESHOLD

print(paste0('Percentage of Passed Regions: %' ,
       round(sum(IC_FILTER_PASSED_INDICES)/length(IC_FILTER_PASSED_INDICES),3)*100))

FilterCounts <- lapply(Counts, function(x) x[IC_FILTER_PASSED_INDICES,])
names(FilterCounts) <- names(Counts)

EV_EXPORTED_THRESHOLD <- 20
EV_FILTER_PASSED_INDICES <- rowSums(FilterCounts[['EV']]) > EV_EXPORTED_THRESHOLD


LabelMat <- data.frame(ic=rowSums(FilterCounts[['IC']]), ev=rowSums(FilterCounts[['EV']]), 
           label=ifelse(EV_FILTER_PASSED_INDICES,'YES','NO'))

write.csv(LabelMat, file = 'Data/TertiaryLabelsMat.csv', quote = F, na = 'NA', row.names = T, col.names = T)






##### making the final design matrix based on new IC_threshold ####

designMatNgram <- read.csv('Data/NgramLen4DesignMatrix.csv',stringsAsFactors = F)
FilterDesignMat <- designMatNgram[IC_FILTER_PASSED_INDICES,]
write.csv(FilterDesignMat, 'Data/filterNgramDesignMat.csv',quote = F, na = 'NA', row.names = T, col.names = T)

## checking if filtered regions are compatible between DesignMatrix and label-Matrix
sum(rownames(FilterCounts[[1]])!=rownames(FilterDesignMat))
sum(rownames(LabelMat)!=rownames(FilterDesignMat))





