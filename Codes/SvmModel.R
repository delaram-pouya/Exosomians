source('Codes/Functions.R')
Initialize()

### loading data
# designMat <- read.csv('Data/filterNgramDesignMat.csv',stringsAsFactors = F)
# LabelMat <- read.csv('Data/TertiaryLabelsMat.csv', stringsAsFactors = F)

TotalMatrixWithStruct <- read.csv('Data/MergedDesignMatLabel_SecondStruct_LenFilter.csv', stringsAsFactors = F)

LabelMat <- subset(TotalMatrixWithStruct, select=c('id','ic','ev','label'))
ColumnsToDrop <- c('id','X','X.1','ic','ev','seq','DB','annotation','rnaType')
designMat <- TotalMatrixWithStruct[,!colnames(TotalMatrixWithStruct) %in% ColumnsToDrop ]
colnames(designMat)



##### Data Preperation 
set.seed(123) 
split = sample.split(designMat$label , SplitRatio = 0.9) 
training_set = subset(designMat, split == TRUE ) 
test_set = subset(designMat, split == FALSE)


colsToBeFactorized = c('chr', 'label', 'strand')
training_set[colsToBeFactorized] = lapply(training_set[colsToBeFactorized], factor)
test_set[colsToBeFactorized] = lapply(test_set[colsToBeFactorized], factor)



################### SVM 
## cross = 10 > cross-validation in the future

SVMlin_Bal = svm(formula= label~., data=training_set , type ='C-classification', kernel = 'linear',
           class.weights= c("NO" = 1, "YES" = 10), probability= T, cachesize=10000)  ## balancing labels

SVMrad_Bal = svm(formula= label~., data=training_set , type ='C-classification', kernel = 'radial',
           class.weights= c("NO" = 1, "YES" = 10), probability= T, cachesize=10000) ## balancing labels + non-linear kernel

SVMlin_NoBal = svm(formula= label~., data=training_set , type ='C-classification', kernel = 'linear',
                   class.weights= c("NO" = 1, "YES" = 1), probability= T, cachesize=10000)

SVMlin_NoStruct = svm(formula= label~., data=training_set[,1:264], type ='C-classification', kernel = 'linear',
                      class.weights= c("NO" = 1, "YES" = 10), probability= T, cachesize=10000)  ## balancing labels

SVMrad_NoStruct = svm(formula= label~., data=training_set[,1:264], type ='C-classification', kernel = 'radial',
                      class.weights= c("NO" = 1, "YES" = 10), probability= T, cachesize=10000)  ## balancing labels




#### loading the trained models 
SVMlin_Bal <- readRDS('models/SVMlin_Bal.rds')
SVMrad_Bal <- readRDS('models/SVMrad_Bal.rds')
SVMlin_NoBal <- readRDS('models/SVMlin_NoBalrds')
SVMlin_NoStruct <- readRDS('models/SVMlin_NoStruct.rds')
SVMrad_NoStruct <- readRDS('models/SVMrad_NoStruct.rds')



### predicion on test data
SVM_linBal_pred <- predict(SVMlin_Bal, test_set)
SVM_radBal_pred <- predict(SVMrad_Bal, test_set)
SVM_linNoBal_pred <- predict(SVMlin_NoBal, test_set)
SVM_linNoStruct_pred <- predict(SVMlin_NoStruct, test_set)
SVM_radNoStruct_pred <- predict(SVMrad_NoStruct, test_set)



#### SVM Model Evaluation
pdf('plots/SVMresults.pdf')
draw_confusion_matrix(confusionMatrix(data = SVM_linBal_pred , reference = test_set$label), 'SVM-weighted(1:10)-linearKernel', 'No','Yes')
draw_confusion_matrix(confusionMatrix(data = SVM_radBal_pred , reference = test_set$label), 'SVM-weighted(1:10)-RadialKernel', 'No','Yes')
draw_confusion_matrix(confusionMatrix(data = SVM_linNoBal_pred , reference = test_set$label), 'SVM-unWeighted-linearKernel', 'No','Yes')
draw_confusion_matrix(confusionMatrix(data = SVM_linNoStruct_pred , reference = test_set$label), 'SVM-weighted(1:10)-noStructFeature-linearKernel', 'No','Yes')
draw_confusion_matrix(confusionMatrix(data = SVM_radNoStruct_pred , reference = test_set$label), 'SVM-weighted(1:10)-noStructFeature-radialKernel', 'No','Yes')
dev.off()

## note: linear kernels did not converge > not able to increase max number of iterations





######### Tuning SVM 

### why doesn't this work??? can't handle the factors :/
SVMrad_Bal_tune <- tune(method = svm, train.x=subset(training_set,select=-label) , train.y=training_set$label, 
                 kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(0.001,0.005,0.01,0.5,1,2))) #epsilon = seq(0,1,0.01)

print(SVMrad_Bal_tune)
svm_model_after_tune <- svm(Species ~ ., data=iris, kernel="radial", cost=1, gamma=0.5)
summary(svm_model_after_tune)








##### Parallel tuning for svm


pkgs <- c('foreach', 'doParallel')
lapply(pkgs, require, character.only = T)
registerDoParallel(cores = 4)
### PREPARE FOR THE DATA ###
df1 <- read.csv("credit_count.txt")
df2 <- df1[df1$CARDHLDR == 1, ]
x <- paste("AGE + ACADMOS + ADEPCNT + MAJORDRG + MINORDRG + OWNRENT + INCOME + SELFEMPL + INCPER + EXP_INC")
fml <- as.formula(paste("as.factor(label) ~ ", x))
### SPLIT DATA INTO K FOLDS ###
set.seed(2016)
df2$fold <- caret::createFolds(1:nrow(df2), k = 4, list = FALSE)
### PARAMETER LIST ###
cost <- c(10, 100)
gamma <- c(1, 2)
parms <- expand.grid(cost = cost, gamma = gamma)
### LOOP THROUGH PARAMETER VALUES ###
result <- foreach(i = 1:nrow(parms), .combine = rbind) %do% {
  c <- parms[i, ]$cost
  g <- parms[i, ]$gamma
  ### K-FOLD VALIDATION ###
  out <- foreach(j = 1:max(df2$fold), .combine = rbind, .inorder = FALSE) %dopar% {
    deve <- df2[df2$fold != j, ]
    test <- df2[df2$fold == j, ]
    mdl <- e1071::svm(fml, data = deve, type = "C-classification", kernel = "radial", cost = c, gamma = g, probability = TRUE)
    pred <- predict(mdl, test, decision.values = TRUE, probability = TRUE)
    data.frame(y = test$label, prob = attributes(pred)$probabilities[, 2])
  }
  ### CALCULATE SVM PERFORMANCE ###
  roc <- pROC::roc(as.factor(out$y), out$prob) 
  data.frame(parms[i, ], roc = roc$auc[1])
}


