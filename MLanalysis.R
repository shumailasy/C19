library(MASS)
library(tidyverse)
library(caret)
library(ggstatsplot)
library(ggplot2)
library(pROC)
library(randomForest)
library(pROC)
library(ROCR)

install.packages("ROCR")
#install.packages("randomForest")
#install.packages("caret")
#install.packages("pROC")
#install.packages("ggstatsplot")

trainData <- read.table("trainData_Coh1_gr.txt", sep="\t", header = T, dec=".", row.names=1)
head(trainData)
testData <- read.table("testData_coh2.txt", sep="\t", header = T, dec=".", row.names=1)
head(testData)

testDataGr <- read.table("testDataGr.txt", sep="\t", header = T, dec=".", row.names=1)
head(testDataGr)

# Linear Discriminant Analysis (LDA) to predict classes
fit <- lda(Group~., data=trainData)
#observed.classes <- testDataGr$Groups
predictions <- predict(fit, testData)
prediction.probabilities <- predictions$posterior[,2]
predicted.classes <- predictions$class
result.lda <- cbind(predictions$class, predictions$posterior)

# Calculate the Confusion Matrix and AUC
observed.classes <- factor(testDataGr$Group)
# Find the accuracy and error of the model
accuracy <- mean(observed.classes == predicted.classes)
accuracy
error <- mean(observed.classes != predicted.classes)
error

# Calculate the ConfusionMatrix form number of cases and proportion of cases
table(observed.classes, predicted.classes)
table(observed.classes, predicted.classes) %>% 
  prop.table() %>% round(digits = 3)

# Print the confusion matrix
confusionMatrix(predicted.classes, 
                observed.classes, 
                positive = "Covid")


# Calculate the AUC 
res.roc <- roc(controls=observed.classes["Healthy"], cases = observed.classes["Covid"])
res.roc <- roc(observed.classes,prediction.probabilities, quiet = T)
res.roc
# Read the ROC data
roc.data <- data_frame(
  thresholds = res.roc$thresholds,
  sensitivity = res.roc$sensitivities,
  specificity = res.roc$specificities
)
roc.data

# Plot the ROC-AUC 
#plot.roc(res.roc, print.auc = TRUE, print.thres = "best")
#legend("topleft", "CD3 vs Esterhuyse TB-contacts data", inset = c(0.05, 0.05))

#plot.roc(res.roc, print.auc = TRUE, print.thres = "best")
#legend("topleft", "CD3 vs  data", inset = c(0.05, 0.05))

# CD3 vs CysFib
plot.roc(res.roc, 
         print.auc = TRUE, 
         print.auc.x = 0.3, 
         print.auc.y = 0.2, 
         print.auc.cex = 2, 
         #print.thres = "best", 
         #print.thres.cex = 2, 
         lwd = 3, cex.lab = 2, cex.axis = 2, family = "Times New Roman")
ClassTable <- cbind(testDataGr$Group, predict(fit, testData)$posterior, data.frame(predict(fit, testData)$class))
ClassTable
write.table(ClassTable, "ML_result_LDA.txt", sep = "\t", quote = F, row.names = T)
# RandomForest

set.seed(123)
model <- train(
  Group ~., data = trainData, method = "rf",
  trControl = trainControl("cv", number = 10),
  importance = TRUE
)
# Best tuning parameter
model$bestTune
model$finalModel
predicted.classes_rf <- model %>% predict(testData)
head(predicted.classes_rf)

predicted.classes_rf

importance(model$finalModel)
ImportantVariables <- varImp(model)
ImportantVariables
ImpVar <- as.matrix(ImportantVariables$importance)
write.table(ImpVar, "Trained_Model_Important_Variables_using_RandomForest.txt", quote = F, row.names = T)

par(mfrow=c(1,2))
# Plot MeanDecreaseAccuracy
varImpPlot(model$finalModel, type = 1)
# Plot MeanDecreaseGini
varImpPlot(model$finalModel, type = 2)


# Find the accuracy and error of the model
accuracy.rf <- mean(observed.classes == predicted.classes_rf)
accuracy.rf
error.rf <- mean(observed.classes != predicted.classes_rf)
error.rf

# Calculate the ConfusionMatrix form number of cases and proportion of cases
table(observed.classes, predicted.classes_rf)
table(observed.classes, predicted.classes_rf) %>% 
  prop.table() %>% round(digits = 3)

confusionMatrix(predicted.classes_rf, 
                observed.classes, 
                positive = "Covid")

result.predicted.prob = predict(model, testData, type = "prob")
result.predicted.prob
covid <- result.predicted.prob[,1]
healthy <- result.predicted.prob[,2]
pROC::roc(testDataGr$Group, covid, quiet = FALSE)
res.roc.rf <- pROC::roc(testDataGr$Group, healthy, quiet = TRUE)
res.roc.rf
plot.roc(res.roc.rf, 
         print.auc = TRUE, 
         print.auc.x = 0.3, 
         print.auc.y = 0.2, 
         print.auc.cex = 2, 
         #print.thres = "best", 
         #print.thres.cex = 2, 
         lwd = 3, cex.lab = 2, cex.axis = 2, family = "Times New Roman")


ClassTable.rf <- cbind(testDataGr$Group, 
                       data.frame(result.predicted.prob), data.frame(predicted.classes_rf))
ClassTable.rf

write.table(ClassTable.rf, "ML_result_randomForest.txt", sep = "\t", quote = F, row.names = T)
