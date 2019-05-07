#  In this file, we will try to find the ngrams of RNA sequences
# to use them as features. We'll make all the possible 4-grams
# and align them to the sequence. The alignment scores will be used as features

#install.packages('gtools')
library(gtools)
library(Biostrings)

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
head(ngramMat)
dim(ngramMat)
summary(rowSums(ngramMat))
write.csv(ngramMat,'/media/pgdrive/users/delaram/data/ExosomeProj/Data/ngramMatrix.csv',quote = F)

designMat2 <- cbind(designMat, ngramMat)
write.csv(designMat2,'/media/pgdrive/users/delaram/data/ExosomeProj/Data/NgramLen4DesignMatrix.csv',quote = F)
###################



tmp <- pairwiseAlignment(pattern=ngrams[2:4], subject=data[4],
                  patternQuality=PhredQuality(22L),
                  subjectQuality=PhredQuality(22L),
                  type="global-local",
                  substitutionMatrix=NULL, fuzzyMatrix=NULL,
                  gapOpening=10, gapExtension=4,
                  scoreOnly=FALSE)
