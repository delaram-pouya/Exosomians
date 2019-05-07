
designMat <- read.csv('Data/NgramLen4DesignMatrix.csv',stringsAsFactors = F)
LabelMat <- read.csv('Data/TertiaryLabelsMat.csv', stringsAsFactors = F)


set.seed(123) 
split = sample.split(DElabels$label , SplitRatio = 0.75) 
training_set = subset(DElabels, split == TRUE ) 
test_set = subset(DElabels, split == FALSE) 
training_set = training_set[,-1]
test_set = test_set[,-1]



classifier = svm(formula = label ~ rowsum, 
                 data = training_set, 
                 type = 'C-classification', 
                 kernel = 'linear') 
classifier