source('Codes/Functions.R')
Initialize()

designMat <- read.csv('Data/filterNgramDesignMat.csv',stringsAsFactors = F)
LabelMat <- read.csv('Data/TertiaryLabelsMat.csv', stringsAsFactors = F)



##### SVM 
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


#### cross = 10 > cross-validation in the future

svm1 = svm(formula= label~., data=training_set , type ='C-classification', kernel = 'linear') 
svm2 = svm(formula= label~., data=training_set , type ='C-classification', kernel = 'linear',
           class.weights= c("NO" = 1, "YES" = 10), probability= T, cachesize=8000)  ## balancing labels

svm3 = svm(formula= label~., data=training_set , type ='C-classification', kernel = 'radial',
           class.weights= c("NO" = 1, "YES" = 10), probability= T, cachesize=10000) ## balancing labels + non-linear kernel


#### Model Evaluation
SVM1_pred <- predict(svm1, test_set)
confusionMatrix(data = SVM1_pred , reference = test_set$label)


#### write the code for cross validation 
###### which parameters to change for svm? kernel, gamma, c , epsilon, ...
#### can we train them in parallel?
### write the code for naive bayes > H2o? 

## naive bayes> h2o?
## RF



























######## parameter tuning 
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


