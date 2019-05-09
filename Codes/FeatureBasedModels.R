source('Codes/Functions.R')
Initialize()
h2o.init()

designMat <- read.csv('Data/filterNgramDesignMat.csv',stringsAsFactors = F)
LabelMat <- read.csv('Data/TertiaryLabelsMat.csv', stringsAsFactors = F)

### bed-file > range of strings 
#### length of strings

set.seed(123) 
split = sample.split(LabelMat$label , SplitRatio = 0.9) 
training_set_label = subset(LabelMat, split == TRUE ) 
training_set_features = subset(designMat, split== TRUE)

test_set_label = subset(LabelMat, split == FALSE) 
test_set_features = subset(designMat, split == FALSE)

training_set = cbind(training_set_features[,-c(1,2,4,11,12)], training_set_label[,'label'])
colnames(training_set)[ncol(training_set)] <- 'label'

test_set = cbind(test_set_features[,-c(1,2,4,11,12)], test_set_label[,'label'])
colnames(test_set)[ncol(test_set)] <- 'label'


colsToBeFactorized = c('chr', 'label', 'strand')
training_set[colsToBeFactorized] = lapply(training_set[colsToBeFactorized], factor)
test_set[colsToBeFactorized] = lapply(test_set[colsToBeFactorized], factor)


pdf('initFeatureModels.pdf')

##### SVM 

#### cross = 10 > cross-validation in the future

svm1 = svm(formula= label~., data=training_set , type ='C-classification', kernel = 'linear') 
svm2 = svm(formula= label~., data=training_set , type ='C-classification', kernel = 'linear',
           class.weights= c("NO" = 1, "YES" = 10), probability= T, cachesize=10000)  ## balancing labels

svm3 = svm(formula= label~., data=training_set , type ='C-classification', kernel = 'radial',
           class.weights= c("NO" = 1, "YES" = 10), probability= T, cachesize=10000) ## balancing labels + non-linear kernel

### prediction on test data
SVM1_pred <- predict(svm1, test_set)
SVM2_pred <- predict(svm2, test_set)
SVM3_pred <- predict(svm3, test_set)

#### SVM Model Evaluation
draw_confusion_matrix(confusionMatrix(data = SVM1_pred , reference = test_set$label))
draw_confusion_matrix(confusionMatrix(data = SVM2_pred , reference = test_set$label))
draw_confusion_matrix(confusionMatrix(data = SVM3_pred , reference = test_set$label))

#fourfoldplot(table(SVM1_pred, test_set$label),color = c("#CC6666", "#99CC99"),conf.level = 0, margin = 1, main = "SVM1-model Confusion Matrix")




##### Naive-Bayes

nb1 <- naiveBayes(formula= label~., data=training_set)
NB1_pred <- predict(nb1, test_set)
confusionMatrix(data = NB1_pred , reference = test_set$label)
draw_confusion_matrix(confusionMatrix(data = NB1_pred , reference = test_set$label))

### h2o implementation:
features = colnames(training_set)[-ncol(training_set)]
label = 'label'
training_set.h2o = training_set
test_set.h2o = test_set
training_set.h2o$weight  <- ifelse(training_set$label=='YES',10,1)
test_set.h2o$weight <- ifelse(test_set$label=='YES', 10, 1)
training_set.h2o <- as.h2o(training_set)
test_set.h2o <- as.h2o(test_set)
nb1 = h2o.naiveBayes(x=features, y = label, training_frame = training_set.h2o, validation_frame = test_set.h2o, balance_classes = T, seed = 1398)
######





####### Random Forest

rf1 = h2o.randomForest(x =features, y = label, training_frame = training_set.h2o,
                       validation_frame = test_set.h2o, balance_classes = F, ntrees = 50, max_depth = 20, seed = 1398)

rf2 = h2o.randomForest(x =features, y = label, training_frame = training_set.h2o,
                      validation_frame = test_set.h2o, balance_classes = T, ntrees = 50, max_depth = 20, seed = 1398)

rf3 = h2o.randomForest(x =features, y = label, training_frame = training_set.h2o,
                       validation_frame = test_set.h2o, balance_classes = T, ntrees = 200, max_depth = 20, seed = 1398)


cm <- as.matrix(h2o.confusionMatrix(rf3))[1:2,1:2]
sensitivity = as.numeric(cm[2, 2])/(as.numeric(cm[2, 2])+ as.numeric(cm[2, 1]))
specificity = as.numeric(cm[1, 1])/(as.numeric(cm[1, 1])+ as.numeric(cm[1, 2]))
Accuracy  = (as.numeric(cm[1, 1])+as.numeric(cm[2, 2]))/sum(as.numeric(cm))

sensitivity
specificity
Accuracy

#fourfoldplot(as.matrix(as.numeric(cm),ncol=2,nrow=2),color = c("#CC6666", "#99CC99"),conf.level = 0, margin = 1, main = "SVM1-model Confusion Matrix")








######## parameter tuning  > SVM 
# perform a grid search
svm_tune <- tune(svm, y ~ x, data = train,
                 ranges = list(epsilon = seq(0,1,0.01), cost = 2^(2:9))
)
print(svm_tune)

#Parameter tuning of ‘svm’:

# - sampling method: 10-fold cross validation 

#- best parameters:
# epsilon cost
#0 8

#- best performance: 2.872047 

#The best model
best_mod <- svm_tune$best.model
best_mod_pred <- predict(best_mod, train) 

error_best_mod <- train$y - best_mod_pred 

# this value can be different on your computer
# because the tune method randomly shuffles the data
best_mod_RMSE <- sqrt(mean(error_best_mod^2)) # 1.290738 

plot(svm_tune)
plot(train,pch=16)
points(train$x, best_mod_pred, col = "blue", pch=4)




