#  In this file, we will try to find the ngrams(K-mers) 
#     of RNA sequences to use them as features. 
#     We'll make all the possible 4-grams and 2-grams
#     then, we'll check the count of the ngrams in each RNA sequence


source('Codes/Functions.R')
Initialize()

N_GRAM_SIZE = 4

permut <- gtools::permutations(n=4, r=N_GRAM_SIZE, v=c('A','T','C','G'),repeats=T)
All_Possible_ngrams <- apply( permut, 1, function(x)paste0(x, collapse=''))



MakeSingleNgrams <- function(String){
  sapply( 1:(nchar(String) - N_GRAM_SIZE+1 ), 
         function(i)
           substr(String, i, i+N_GRAM_SIZE-1 )
  )
}

designMat = read.csv('Data/PrimaryDesignMat.csv')
Seqs <- designMat$seq

Feature_Ngrams <- sapply(Seqs,MakeSingleNgrams )
summary(sapply(Feature_Ngrams, length))

#ZeroOneMat <- sapply(1:length(Feature_Ngrams), 
#              function(i)ifelse(All_Possible_ngrams %in% Feature_Ngrams[[i]],1,0),simplify = F)

ngramMat <- sapply(1:length(Feature_Ngrams), 
             function(i) match(All_Possible_ngrams, Feature_Ngrams[[i]], nomatch = 0),simplify = F)

ngramMat <- as.data.frame(do.call(rbind, ngramMat))
colnames(ngramMat) <- All_Possible_ngrams
rownames(ngramMat) <- designMat$id

summary(rowSums(ngramMat))
write.csv(ngramMat,'Data/ngramMatrix.csv',quote = F)

designMatNgram <- cbind(designMat, ngramMat)
write.csv(designMatNgram,'Data/NgramLen4DesignMatrix.csv',quote = F)

###################


