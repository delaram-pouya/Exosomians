## In this file we write the sequences of our regions  
##  of interest in a fasta file, which is needed by 
##  ViennaRNA Package for RNA Secondary structure prediction

#   input: merged label and design matrix
#   output: fasta file


TotalMatrix <- read.csv('Data/MergedDesignMatLabel_LenFilter.csv',stringsAsFactors=F )
Sequence <- data.frame(id=TotalMatrix$id,seq=TotalMatrix$seq)


writeFasta<-function(data, filename){
  fastaLines = c()
  
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"id"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

head(Sequence)
writeFasta(Sequence, "Data/RegionsLenFilter.fasta")





##################################
# Secondary structure prediction using ViennaRNA tool

### bash 
# RNAfold  -i RegionsLenFilter.fasta -o RNAseconStructPredict.txt
# grep '(\|>' RNAseconStructPredict.txt > RNAseconStructPredict_withoutSeq.txt

RNAseqFasta <- readDNAStringSet('Data/RegionsLenFilter.fasta')
SeqID = names(RNAseqFasta)
RNAseqLength <- nchar(paste(RNAseqFasta))


SecondaryStructureWithEnergy = readBStringSet("Data/RNAseconStructPredict_withoutSeq.txt")
sum(names(SecondaryStructureWithEnergy)!=SeqID)

SecondaryStructureWithEnergy <- paste(SecondaryStructureWithEnergy)  
DotBracketSecondStruct  <-  substr(SecondaryStructureWithEnergy,1, RNAseqLength)

FreeEnergy <- gsub('[() ]','', 
                   substr(SecondaryStructureWithEnergy, RNAseqLength,nchar(SecondaryStructureWithEnergy)))

extraDotIndex <- substr(FreeEnergy,1,1)=='.'
FreeEnergy[extraDotIndex] = as.numeric(substr(FreeEnergy[extraDotIndex],2,
                                              nchar(FreeEnergy[extraDotIndex])))
head(FreeEnergy)





##################### 
# RNA motif annotation

##  convert dot_bracket notation > Ct
# bash
# RNAfold < ~/delaram/data/ExosomeProj/Data/RNAseconStructPredict.txt | b2ct > ct_file.ct


library(RNAstructureModuleMiner)
library(RRNA)

CtSecondStruct <- sapply(1:nrow(Sequence), 
                        function(i) makeCt(DotBracketSecondStruct[i], Sequence[i,'seq']),
                        simplify = F)


names(CtSecondStruct) <- Sequence$id




## test 
ct=makeCt("(((...(((...)))...(((...)))...)))","AAAUUUCCCAAAGGGUUUAAAGGGUUUCCCUUU")
coord=ct2coord(ct)
RNAPlot(coord,hl=c("GGGUUU","AAAUUU"),seqcols=c(2,4),labTF=TRUE)

col1 <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
col2 <- c("G","C","C","A","C","C","C","U","G","C","A","G","G","G","U","C","G","G","C")
col3 <- c(0,1,2,3,4,5,6,7,8,0,10,11,12,13,14,15,16,17,18)
col4 <- c(2,3,4,5,6,7,8,9,0,11,12,13,14,15,16,17,18,19,0)
col5 <- c(19,18,17,15,14,13,12,11,10,9,8,7,6,5,4,0,3,2,1)
col6 <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
data <- matrix(c(col1,col2,col3,col4,col5,col6),byrow=FALSE,ncol =6)

bulge_loop(data)
hairpin_loop(data)
internal_loop(data)
multi_branch_loop(data)
stem(data)





