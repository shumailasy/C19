### TB Biosignature analysis with missMethyl

library(missMethyl)
library(limma)
library(minfi)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("minfiData")
library(minfiData)

baseDir <- "/proj/methylation/SignatureAnalysis/HLA-DR/"
targets <- read.metharray.sheet(baseDir)

rgSet <- read.metharray.exp(targets = targets)
mSet <- preprocessRaw(rgSet)
mSetSw <- SWAN(mSet,verbose=TRUE)

detP <- detectionP(rgSet)
keep <- rowSums(detP < 0.01) == ncol(rgSet)
mSetSw <- mSetSw[keep,]

mset_reduced <- mSetSw[sample(1:nrow(mSetSw)),]
meth <- getMeth(mset_reduced)
unmeth <- getUnmeth(mset_reduced)
Mval <- log2((meth + 100)/(unmeth + 100))
beta <- getBeta(mset_reduced)
dim(Mval)
targets$Sample_Group
group <- factor(targets$Sample_Group,levels=c("Act_LatTB","NE_negIGRA"))
#id <- factor(targets$person)
design <- model.matrix(~group)
design

fit.reduced <- lmFit(Mval,design)
fit.reduced <- eBayes(fit.reduced)
summary(decideTests(fit.reduced))
top<-topTable(fit.reduced,coef=2)
top
head(beta)
colnames(beta) <- targets$Sample_Name
write.table(beta, "BetaValues_LinköpingLima_7DMCs_CtrlIGRA-VsPatIGRA+.txt", sep = "\t", quote = F)


cpgs <- rownames(top)
par(mfrow=c(2,4))
for(i in 1:8){
  stripchart(beta[rownames(beta)==cpgs[i],]~design[,2],method="jitter",
             group.names=c("LA_TB","CtrlNeg"),pch=c(2, 19),cex=1.5,col=c(4,2),ylab=expression(beta~ "values"),
             vertical=TRUE,cex.axis=1.5,cex.lab=1.5)
  title(cpgs[i],cex.main=1.5)
}



# ChAMP analysis
#library("ChAMP")
myLoad <- champ.load(arraytype="450K")

myNorm <- champ.norm(beta=myLoad$beta, rgSet=myLoad$rgSet,
                     mset=myLoad$mset,
                     resultsDir="./CHAMP_Normalization/", method="BMIQ",
                     plotBMIQ=FALSE,
                     arraytype="450K",
                     cores=24)
# Save the normalized file
write.table(myNorm, "Norm_LinköpingLima_CtrlIGRA-VsPatIGRA+.txt",
            sep = "\t",
            quote = F,
            row.names = T)
# Calculate the Differential Methylation CpG sites
myDMC <- champ.DMP(beta = myNorm,
                   pheno = myLoad$pd$Sample_Group,
                   compare.group = NULL,
                   adjPVal = 0.05,
                   adjust.method = "BH",
                   arraytype = "450K")

write.table(myDMC[[1]], "DMCs_LinköpingLima_CtrlIGRA-VsPatIGRA+.txt",
            sep = "\t",
            quote = F,
            row.names = T)




#test
group <- factor(targets$Sample_Group,levels=c("NE_negIGRA", "Act_LatTB"))
#id <- factor(targets$person)
design <- model.matrix(~group)
design

fit.reduced <- lmFit(Mval,design)
fit.reduced <- eBayes(fit.reduced)
summary(decideTests(fit.reduced))
top<-topTable(fit.reduced,coef=2)
top

cpgs <- rownames(top)
par(mfrow=c(2,4))
for(i in 1:8){
  stripchart(beta[rownames(beta)==cpgs[i],]~design[,2],method="jitter",
             group.names=c("CtrlNeg", "LA_TB"),pch=c(2, 19),cex=1.5,col=c(4,2),ylab=expression(beta~ "values"),
             vertical=TRUE,cex.axis=1.5,cex.lab=1.5)
  title(cpgs[i],cex.main=1.5)
}


# AUC from the ROC curve
library(MASS)
library(tidyverse)
library(caret)
library(ggstatsplot)
library(ggplot2)
library(pROC)

trainData <- read.table("trainData20_all3520.txt", sep="\t", header = T, dec=".", row.names=1)
#head(trainData)
testData <- read.table("testData20_all3520.txt", sep="\t", header = T, dec=".", row.names=1)
#head(testData)
testDataGr <- read.table("testDataGr20_all3520.txt", sep="\t", header = T, dec=".", row.names=1)
#head(testDataGr)
fit <- lda(Groups~., data=trainData)
observed.classes <- testDataGr$Groups
predictions <- predict(fit, testData)
prediction.probabilities <- predictions$posterior[,2]
predicted.classes <- predictions$class
table(observed.classes, predicted.classes)
#                predicted.classes
#observed.classes neg pos
#neg               14   6
#pos               0   7

ClassTable <- cbind(testDataGr$Groups, predict(fit, testData)$posterior, data.frame(predict(fit, testData)$class))
ClassTable

rsAVG <- multiclass.roc(testDataGr$Groups, prediction.probabilities)
rsAVG

confusionMatrix(predictions$class, observed.classes, positive = "pos")
#confusionMatrix(predictions$class, testDataGr$Groups)
accuracy <- mean(observed.classes == predicted.classes)
accuracy
error <- mean(observed.classes != predicted.classes)
error

res.roc <- roc(observed.classes, prediction.probabilities)
res.roc
#Call:
#  roc.default(response = observed.classes, predictor = prediction.probabilities)

#Data: prediction.probabilities in 20 controls (observed.classes neg) < 7 cases (observed.classes pos).
#Area under the curve: 0.9

plot.roc(res.roc, 
         print.auc = TRUE, 
         print.auc.x = 0.3, 
         print.auc.y = 0.2, 
         print.auc.cex = 2, 
         #print.thres = "best", 
         #print.thres.cex = 2, 
         lwd = 3, cex.lab = 2, cex.axis = 2, family = "Times New Roman")




set.seed(123)
model <- train(
  Groups ~., data = trainData, method = "rf",
  trControl = trainControl("cv", number = 10),
  importance = TRUE
)
# Best tuning parameter
model$bestTune
model$finalModel

# Make predictions on the test data
predicted.classes <- model %>% predict(testData)
head(predicted.classes)
mean(predicted.classes == testDataGr$Groups)
RMSE(predicted.classes, testDataGr$Groups)

res.roc <- roc(observed.classes, predicted.classes)



# Variable Importance
library(tidyverse)
library(caret)
library(randomForest)
importance(model$finalModel)

par(mfrow=c(1,2))
# Plot MeanDecreaseAccuracy
varImpPlot(model$finalModel, n.var = 30, type = 1)
# Plot MeanDecreaseGini
varImpPlot(model$finalModel, type = 2)

varImp(model$finalModel)
nrow(varImp(model)$importance)
varImpList <- varImp(model)$importance

varImp(model)
, n = 30)
, scale = F)


importance(model$finalModel, type = 1)
head(model$finalModel$importance)





# AUC calculations

library(MASS)
library(tidyverse)
library(caret)
#library(ggstatsplot)
library(ggplot2)
library(pROC)

trainData <- read.table("trainData20_all3520.txt", sep="\t", header = T, dec=".", row.names=1)
testData <- read.table("testData20_Data2.txt", sep="\t", header = T, dec=".", row.names=1)
testDataGr <- read.table("testDataGr20_Data2.txt", sep="\t", header = T, dec=".", row.names=1)
fit <- lda(Groups~., data=trainData)
observed.classes <- testDataGr$Groups
predictions <- predict(fit, testData)
prediction.probabilities <- predictions$posterior[,2]
predicted.classes <- predictions$class
table(observed.classes, predicted.classes)
ClassTable <- cbind(testDataGr$Groups, predict(fit, testData)$posterior, data.frame(predict(fit, testData)$class))
ClassTable
confusionMatrix(predictions$class, observed.classes, positive = "pos")
accuracy <- mean(observed.classes == predicted.classes)
accuracy
error <- mean(observed.classes != predicted.classes)
error
res.roc <- roc(observed.classes, prediction.probabilities)
res.roc
plot.roc(res.roc, 
         print.auc = TRUE, 
         print.auc.x = 0.3, 
         print.auc.y = 0.2, 
         print.auc.cex = 2, 
         #print.thres = "best", 
         #print.thres.cex = 2, 
         lwd = 3, cex.lab = 2, cex.axis = 2, family = "Times New Roman")






# Interactive Heatmap 
library(heatmaply)
heatmaply(as.matrix(t(hm)), 
          col_side_colors = heatmap$Groups,
          showticklabels = c(TRUE, FALSE),
          scale = "none")












