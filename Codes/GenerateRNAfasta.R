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
