#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("missMethyl")
#BiocManager::install("minfiData")
library(missMethyl)
library(limma)
library(minfi)
library(minfiData)

baseDir <- "/Users/jyoda68/Documents/ML/COVID-19/NewData_201025/Analysis_201026/"
baseDir <- "/Users/jyoda68/Documents/ML/For_Patent/ServerDoc/SignatureAnalysis/HLA-DR"
baseDir <- "/Users/jyoda68/Documents/ML/COVID-19/NewData20201113/Analysis"
baseDir <- "/Users/jyoda68/Documents/ML/COVID-19/NewData_201025/Analysis_201026/"
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

group <- factor(targets$Sample_Group,levels=c("CON","COVID"))
#id <- factor(targets$person)
design <- model.matrix(~group)
design

fit.reduced <- lmFit(Mval,design)
fit.reduced <- eBayes(fit.reduced)

summary(decideTests(fit.reduced))

top<-topTable(fit.reduced,coef=2)
top
tab1 <- topTable(fit.reduced, n=Inf, coef=2)
head(tab1)
tab <- as.data.frame(fit.reduced)
write.table(tab1, "DMCs_COVID19-1.txt", sep = "\t",quote = F )
# Plot beta values for two CpG sites
cpgs <- rownames(top)
par(mfrow=c(1,2))
for(i in 1:2){
  stripchart(beta[rownames(beta)==cpgs[i],]~design[,2],method="jitter",
             group.names=c("Control","COVID-19"),pch=19,cex=1.2,col=c(3.5,2),
             ylab= c(expression(paste(beta)~"values"),family="Times New Roman"),
             vertical=TRUE,cex.axis=1.2,cex.lab=1.2, family="Times New Roman")
  title(cpgs[i],cex.main=1.2, family="Times New Roman")
}


# Save the Beta and M value files
colnames(beta) <- targets$Sample_Name
colnames(Mval) <- targets$Sample_Name
write.table(beta, "BetaValueData_COVID-19.txt", sep = "\t", quote = F)
write.table(Mval, "MethylationValueData_COVID-19.txt", sep = "\t", quote = F)


# Violin plot of significant CpGs
dat880 <- read.table("cg880.txt", stringsAsFactors = F, header = T)
dat070 <- read.table("cg070.txt", stringsAsFactors = F, header = T)
con <- grepl("Control", dat880$group)
p <- ggplot(dat880, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)+
  stat_compare_means(method = "t.test")
p1 <- 
  ggplot(dat070, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=NULL, x=NULL)+
  stat_compare_means(method = "t.test")

#install.packages("ggpubr")
#library(ggpubr)

ggarrange(p, p1, 
          labels = c("cg...880", "cg...070"),
          ncol = 2, nrow = 1, align = "v", hjust = -2.5, vjust = 2, common.legend = TRUE)


p <- ggplot(dat880, aes(x=group, y=beta, fill=group)) + 
  geom_boxplotplot(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)+
  stat_compare_means(method = "t.test", label = "p.signif", label.y = 0.9)+
  stat_compare_means(label.y = 0.7)
p1 <- 
  ggplot(dat070, aes(x=group, y=beta, fill=group)) + 
  geom_boxplotplot(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=NULL, x=NULL)+
  stat_compare_means(method = "t.test")
ggarrange(p, p1, 
          labels = c("cg...880", "cg...070"),
          ncol = 2, nrow = 1, align = "v", hjust = -2.5, vjust = 2, common.legend = TRUE)

### Boxplot For TB-signature analysis
#### HLA-DR cell types:: Control+IGRA negative Vs IGRA+

boxData <- read.table("ForBox.txt", sep = ",", header = T, row.names = 1)
head(boxData)
boxplot(boxData)

#library(ggplot2)
#library(ggpubr)
# Reading CpGs beta values
boxData1 <- read.table("ForBox1_008.txt", stringsAsFactors = F, header = T)
boxData2 <- read.table("ForBox1_012.txt", stringsAsFactors = F, header = T)
boxData3 <- read.table("ForBox1_300.txt", stringsAsFactors = F, header = T)
boxData4 <- read.table("ForBox1_375.txt", stringsAsFactors = F, header = T)
boxData5 <- read.table("ForBox1_559.txt", stringsAsFactors = F, header = T)
boxData6 <- read.table("ForBox1_627.txt", stringsAsFactors = F, header = T)
boxData7 <- read.table("ForBox1_642.txt", stringsAsFactors = F, header = T)
boxData8 <- read.table("ForBox1_810.txt", stringsAsFactors = F, header = T)
boxData9 <- read.table("ForBox1_822.txt", stringsAsFactors = F, header = T)
boxData10 <- read.table("ForBox1_894.txt", stringsAsFactors = F, header = T)
boxData11 <- read.table("ForBox1_956.txt", stringsAsFactors = F, header = T)
con <- grepl("Control_IGRAneg", boxData1$group)
p1 <- ggplot(boxData1, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p2 <- ggplot(boxData2, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p3 <- ggplot(boxData3, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p4 <- ggplot(boxData4, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p5 <- ggplot(boxData5, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p6 <- ggplot(boxData6, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p7 <- ggplot(boxData7, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p8 <- ggplot(boxData8, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p9 <- ggplot(boxData9, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p10 <- ggplot(boxData10, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p11 <- ggplot(boxData11, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)

ggarrange(p1, p2,p3, p4, p5, p6, p7, p8, p9, p10,p11, 
          labels = c("cg...008", "cg...012", "cg...300", "cg...375", "cg...559",
                     "cg...627", "cg...642", "cg...810", "cg...822", "cg...894", "cg...956"),
          ncol = 6, nrow = 2, align = "v", hjust = -2.5, vjust = 2)

#### HLA-DR cell types:: Latent TB,IGRA+ Vs Active TB, patients
box_618 <- read.table("ForBox_618.txt", stringsAsFactors = F, header = T)
box_445 <- read.table("ForBox_445.txt", stringsAsFactors = F, header = T)
expp <- grepl("Latent_TB", box_445$group)
p618 <- ggplot(box_618, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(expp, "red", "blue"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("red","orange"))+
  labs(y=beta~"value", x=NULL)+
  stat_compare_means(method = "t.test")
p445 <- ggplot(box_445, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(expp,  "red", "blue"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("red","orange"))+
  labs(y=beta~"value", x=NULL)+
  stat_compare_means(method = "t.test")
ggarrange(p618, p445,
          labels = c("cg...618", "cg...445"),
          ncol = 2, nrow = 1, align = "v", hjust = -2.5, vjust = 2, common.legend = TRUE)

#### CD3 cell types:: Latent TB, IGRA+ vs Active TB, patients

box_685 <- read.table("box_685.txt", stringsAsFactors = F, header = T)
box_766 <- read.table("box_766.txt", stringsAsFactors = F, header = T)
box_956 <- read.table("box_956.txt", stringsAsFactors = F, header = T)
expp <- grepl("Latent_TB", box_685$group)
p685 <- ggplot(box_685, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(expp, "red", "blue"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("red","orange"))+
  labs(y=beta~"value", x=NULL)+
  stat_compare_means(method = "t.test", label.y = 0.75)
p766 <- ggplot(box_766, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(expp, "red", "blue"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("red", "orange"))+
  labs(y=beta~"value", x=NULL)+
  stat_compare_means(method = "t.test", label.y = 0.75)
p956 <- ggplot(box_956, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(expp, "red", "blue"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("red", "orange"))+
  labs(y=beta~"value", x=NULL)+
  stat_compare_means(method = "t.test", label.y = 0.78)
ggarrange(p685, p766,p956,
          labels = c("cg...685", "cg...766", "cg...956"),
          ncol = 3, nrow = 1, align = "v", hjust = -2.5, vjust = 2, common.legend = TRUE)


## New analysis on 2020-11-23 with the pre-destined genes
dmcs <- read.table("CpGs_GeneSymbol_sorted1.txt", stringsAsFactors = F, header = T, fill = T)
head(dmcs)
gs <- read.table("GeneSymbol_curated.txt", stringsAsFactors = F, header = T)
head(gs)
gs_dmcs <- merge(gs, dmcs, by = "GeneSymbol")
head(gs_dmcs)
write.table(gs_dmcs, "curated_dmcs.txt", sep = "\t", row.names = F, quote = F)

beta <- read.table("BetaValueData_COVID-19.txt", stringsAsFactors = F, header = T)
head(beta)
beta1 <- cbind(rownames(beta), beta) 
head(beta1)
rownames(beta1) <- NULL
colnames(beta1)[1] <- "CpGID"
head(beta1)
gs_dmcs_beta <- merge(gs_dmcs, beta1, by = "CpGID")
head(gs_dmcs_beta)
write.table(gs_dmcs_beta, "CuratedDataBetaVal.txt", sep="\t", row.names = F, quote = F)


baseDir <- "/Users/jyoda68/Documents/ML/COVID-19/NewData_201025/Analysis_201026/"
targets <- read.metharray.sheet(baseDir)
rgSet <- read.metharray.exp(targets = targets)
mSet <- preprocessRaw(rgSet)
detP <- detectionP(rgSet)
keep <- rowSums(detP < 0.01) == ncol(rgSet)
table(keep)
mSetSw <- mSet[keep,]
mset_reduced <- mSetSw[sample(1:nrow(mSetSw)),]
meth <- getMeth(mset_reduced)
unmeth <- getUnmeth(mset_reduced)
Mval <- log2((meth + 100)/(unmeth + 100))
beta <- getBeta(mset_reduced)
dim(Mval)
targets$Sample_Group
group <- factor(targets$Sample_Group,levels=c("CON","COVID"))
#id <- factor(targets$person)
design <- model.matrix(~group)
design
fit.reduced <- lmFit(Mval,design)
fit.reduced <- eBayes(fit.reduced)
summary(decideTests(fit.reduced))
top<-topTable(fit.reduced,coef=2)
top
tab <- as.data.frame(fit.reduced)
head(tab)
tab1 <- topTable(fit, n=Inf, coef=2)
tab1 <- topTable(fit.reduced, n=Inf, coef=2)
head(tab1)
anno <- getAnnotation(tab1)
anno <- getAnnotation(fit.reduced)
anno <- getManifest(fit.reduced)
rSet <- ratioConvert(mSet, what="both", keepCN =T)
rSet
beta <- getBeta(rSet)
GRset <- mapToGenome(rSet)
GRset
beta <- getBeta(GRset)
M <- getM(GRset)
CN <- getCN(GRset)
annot <- getAnnotation(GRset)
head(annot)
dim(annot)
dim(tab1)
head(tab1, n =10)
tab1_annot <- cbind(tab1, annot)
tab1_annot <- annotation(tab1)
tab1_annot <- annotation[match(rownames(tab1), rownames(annot)),]
tab1_annot <- annotation[match(rownames(tab1), rownames(annot))]
tab1_anno <- annot[match(rownames(M), annot$Name), c(1:4,12:19,24:ncol(annot))]
head(tab1_anno)
topAnno <- topTable(fit.reduced, num =Inf, coef = 2, genelist = tab1_anno)
head(topAnno)
dim(topAnno)
write.table(topAnno, "DMCs_COVID19_Data1_ALL.txt", sep = "\t", quote = F)


####
betaAll <- read.table("../../COVIDdata-12/BetaVal_ALL_COVID_12.txt", stringsAsFactors = F, header = T)
head(betaAll)
betaAll1 <- cbind(rownames(betaAll), betaAll) 
head(betaAll1)
rownames(betaAll1) <- NULL
colnames(betaAll1)[1] <- "CpGID"
head(betaAll1)
head(gs_dmcs)
gs_betaAll <- merge(gs_dmcs, betaAll1, by = "CpGID")
head(gs_betaAll)
write.table(gs_dmcs_beta, "CuratedDataBetaVal.txt", sep="\t", row.names = F, quote = F)


# Machine Learning analysis
## AUC from the ROC Curve
library(MASS)
library(tidyverse)
library(caret)
library(ggstatsplot)
library(ggplot2)
library(pROC)

trainData <- read.table("trainDataAll.txt", sep="\t", header = T, dec=".", row.names=1)
#head(trainData)
testData <- read.table("testDataAll.txt", sep="\t", header = T, dec=".", row.names=1)
#head(testData)
testDataGr <- read.table("testDataGr.txt", sep="\t", header = T, dec=".", row.names=1)
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

confusionMatrix(predictions$class, testDataGr$Groups)
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


# Assign CpG ID to the Gene Symbol of Marias List
## Take all DMCs (with P = 1)
### Read Illumina EPIC array annotation file
ilEPICdata <- read.csv("../../infinium-methylationepic-v-1-0-b5-manifest-file.csv")
head(ilEPICdata)
ilmData <- ilEPICdata[,c(1, 16)]
head(ilmData)

# Lollypop plot
library(ggplot2)
library(dplyr)
data <- read.table("varImportance_GeneSymbol.txt", stringsAsFactors = F, header = T)
data %>%
  arrange(Importance) %>%
  mutate(name=factor(CpGID_GeneSymbol, levels = CpG_GeneSymbol)) %>%
  ggplot(name=factor(x=CpGID_GeneSymbol, y = Importance)) %>%
  geom_segment(aes(xend=CpGID_GeneSymbol, yend=0)) +
  geom_point(size =4, color= "grey") +
  coord_flip()+
  theme_bw()+
  xlab("")
  
ggplot(data, aes(CpGID_eneSymbol, Importance)) +
  geom_segment(aes(x = 0, y = Importance, xend = CpGID_eneSymbol, yend = Importance), color = "grey50") +
  geom_point()

data2 <- read.table("Top20_dat12LPIpbmc.txt", stringsAsFactors = F, header = T)
head(data2)
ggplot(data2, aes(x=CpGsID_GeneSymbol, y=Importance)) +
  geom_segment( aes(x=CpGsID_GeneSymbol, xend=CpGsID_GeneSymbol, y=0, yend=Importance), color="grey") +
  geom_point( color="black", size=4, alpha=1) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )


# Fine tuned Lollypop Plot
library(tidyverse)
library(hrbrthemes)
library(kableExtra)
options(knitr.table.format = "html")
library(patchwork)

# Load dataset from github
data <- read.table("varImportance_GeneSymbol.txt", stringsAsFactors = F, header = T)

# Plot
data2 %>%
  filter(!is.na(Importance)) %>%
  arrange(Importance) %>%
  tail(20) %>%
  mutate(CpGsID_GeneSymbol=factor(CpGsID_GeneSymbol, CpGsID_GeneSymbol)) %>%
  ggplot( aes(x=CpGsID_GeneSymbol, y=Importance) ) +
  geom_segment( aes(x=CpGsID_GeneSymbol ,xend=CpGsID_GeneSymbol, y=0, yend=Importance), color="black") +
  geom_point(size=3, color="black") +
  coord_flip() +
  theme_ipsum() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("CpG ID and Gene Symbol") +
  ylab("Variable Importance") +
  theme(text=element_text(size = 12, family = "Times New Roman"))+
  theme(axis.text.x = element_text(size = 12, family = "Times New Roman", colour = "black"),
        axis.text.y = element_text(size = 12, family = "Times New Roman", colour = "black"))


# Mathced top 20 CpGs from 3520 data (mixture of 850K and 450K)
cpg20 <- read.table("CpGsID20_from3520.txt", stringsAsFactors = F, header = T)
allData <- read.table("AllData_dat12_LPIPBMC3520.txt", stringsAsFactors = F, header = T)
cg20AllData <- merge(cpg20, allData, by = "CpGID")
head(cpg20)
head(allData)
colnames(cpg20) <- "CpGID"
dim(cg20AllData)
write.table(cg20AllData, "CpG20_AllData.txt", sep = "\t", row.names = F, quote = F)




##### Violin plot 20 CpGs (derived from 3520 CpGs after matched with 450K dataset) of COVID-19 dataset 1.
library(ggplot2)
library(ggpubr)
# Reading CpGs beta values
boxData1 <- read.table("cg00141548.txt", stringsAsFactors = F, header = T)
boxData2 <- read.table("cg00954771.txt", stringsAsFactors = F, header = T)
boxData3 <- read.table("cg01593751.txt", stringsAsFactors = F, header = T)
boxData4 <- read.table("cg02940070.txt", stringsAsFactors = F, header = T)
boxData5 <- read.table("cg03466841.txt", stringsAsFactors = F, header = T)
boxData6 <- read.table("cg04920428.txt", stringsAsFactors = F, header = T)
boxData7 <- read.table("cg09015880.txt", stringsAsFactors = F, header = T)
boxData8 <- read.table("cg10171267.txt", stringsAsFactors = F, header = T)
boxData9 <- read.table("cg12149606.txt", stringsAsFactors = F, header = T)
boxData10 <- read.table("cg12478978.txt", stringsAsFactors = F, header = T)
boxData11 <- read.table("cg13847858.txt", stringsAsFactors = F, header = T)
boxData12 <- read.table("cg14070162.txt", stringsAsFactors = F, header = T)
boxData13 <- read.table("cg14415283.txt", stringsAsFactors = F, header = T)
boxData14 <- read.table("cg23307218.txt", stringsAsFactors = F, header = T)
boxData15 <- read.table("cg23757825.txt", stringsAsFactors = F, header = T)
boxData16 <- read.table("cg24321706.txt", stringsAsFactors = F, header = T)
boxData17 <- read.table("cg24874557.txt", stringsAsFactors = F, header = T)
boxData18 <- read.table("cg26181929.txt", stringsAsFactors = F, header = T)
boxData19 <- read.table("cg26458816.txt", stringsAsFactors = F, header = T)
boxData20 <- read.table("cg27123975.txt", stringsAsFactors = F, header = T)
con <- grepl("neg", boxData1$group)
p1 <- ggplot(boxData1, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p2 <- ggplot(boxData2, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p3 <- ggplot(boxData3, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p4 <- ggplot(boxData4, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p5 <- ggplot(boxData5, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p6 <- ggplot(boxData6, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p7 <- ggplot(boxData7, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p8 <- ggplot(boxData8, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p9 <- ggplot(boxData9, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p10 <- ggplot(boxData10, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p11 <- ggplot(boxData11, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p12 <- ggplot(boxData12, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p13 <- ggplot(boxData13, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p14 <- ggplot(boxData14, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p15 <- ggplot(boxData15, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p16 <- ggplot(boxData16, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p17 <- ggplot(boxData17, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p18 <- ggplot(boxData18, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p19 <- ggplot(boxData19, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p20 <- ggplot(boxData20, aes(x=group, y=beta, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)

ggarrange(p1, p2,p3, p4, p5, p6, p7, p8, p9, p10,p11, p12,p13, p14, p15, p16, p17, p18, p19, p20,
          labels = c("cg00141548", "cg00954771", "cg01593751", "cg02940070", "cg03466841", "cg04920428",
                     "cg09015880", "cg10171267", "cg12149606", "cg12478978", "cg13847858", "cg14070162",
                     "cg14415283", "cg23307218", "cg23757825", "cg24321706", "cg24874557", "cg26181929",
                     "cg26458816", "cg27123975"),
          ncol = 5, nrow = 4, align = "v", hjust = -2.5, vjust = 2)



#### Boxplot
p1 <- ggplot(boxData1, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p2 <- ggplot(boxData2, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p3 <- ggplot(boxData3, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p4 <- ggplot(boxData4, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p5 <- ggplot(boxData5, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p6 <- ggplot(boxData6, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p7 <- ggplot(boxData7, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p8 <- ggplot(boxData8, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p9 <- ggplot(boxData9, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p10 <- ggplot(boxData10, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p11 <- ggplot(boxData11, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p12 <- ggplot(boxData12, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p13 <- ggplot(boxData13, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p14 <- ggplot(boxData14, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p15 <- ggplot(boxData15, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p16 <- ggplot(boxData16, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p17 <- ggplot(boxData17, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p18 <- ggplot(boxData18, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p19 <- ggplot(boxData19, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)
p20 <- ggplot(boxData20, aes(x=group, y=beta, fill=group)) + 
  geom_boxplot()+
  geom_jitter(position = position_jitter(0.1),color = ifelse(con, "darkgreen", "red"), size =3)+
  theme_classic()+
  scale_fill_manual(values = c("lightgreen", "orange"))+
  labs(y=beta~"value", x=NULL)

ggarrange(p1, p2,p3, p4, p5, p6, p7, p8, p9, p10,p11, p12,p13, p14, p15, p16, p17, p18, p19, p20,
          labels = c("cg00141548", "cg00954771", "cg01593751", "cg02940070", "cg03466841", "cg04920428",
                     "cg09015880", "cg10171267", "cg12149606", "cg12478978", "cg13847858", "cg14070162",
                     "cg14415283", "cg23307218", "cg23757825", "cg24321706", "cg24874557", "cg26181929",
                     "cg26458816", "cg27123975"),
          ncol = 5, nrow = 4, align = "v", hjust = -2.5, vjust = 2)
