#17000 test sequences - the tools only take 500 seq per time
#Need to copy-paste sequences
#Take 2000 sample from IC and 2000 from EV
#Feed 4 times for each label
#Parse out the sequence

library(seqinr)
library(data.table)
set.seed(10)
SAMPLING_NUM = 2000
NUM_SEQ = 500

write_split_seqs_to_fasta <- function(seq, names, num_seq, file_name){
  start = 1
  for(i in 1:(length(seq)/num_seq)){
    
    sequences_to_write = seq[start:(start+num_seq-1) ]
    seq_names_to_write = names[start:(start+num_seq-1) ]
    
    file_name_to_write = paste0('Data/exoGRU_seqs/split_test_seqs/test_',file_name,'_subsample_', start,'_' ,(start+num_seq-1), '.fasta')
    start = start+num_seq
    
    write.fasta(sequences = as.list(sequences_to_write), 
                names=seq_names_to_write, 
                file=file_name_to_write)
  }
}

data_df = read.csv('Data/exoGRU_seqs/ExoGRU_predictions_splitted.csv')
table(data_df$split)
test_EV <- data_df[data_df$split == 'test' & data_df$label == 'EV',]
test_IC <- data_df[data_df$split == 'test' & data_df$label == 'IC',]


#######################################################
################ splitting the IC sequences and writing into fasta files
test_IC.dt <- data.table(test_IC)
test_IC_subsample = data.frame(test_IC.dt[sample(.N, SAMPLING_NUM)])
dim(test_IC_subsample)

write_split_seqs_to_fasta(seq = test_IC_subsample$seq,
                          names = test_IC_subsample$id,
                          num_seq = NUM_SEQ,
                          file_name = 'IC')


#######################################################
################ splitting the EV sequences and writing into fasta files
test_EV.dt <- data.table(test_EV)
test_EV_subsample = data.frame(test_EV.dt[sample(.N, SAMPLING_NUM)])
dim(test_EV_subsample)

start = 1
write_split_seqs_to_fasta(seq = test_EV_subsample$seq,
                          names = test_EV_subsample$id,
                          num_seq = NUM_SEQ,
                          file_name = 'EV')


head(test_IC_subsample)
head(test_EV_subsample)
