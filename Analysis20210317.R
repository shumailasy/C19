

dataDir <- setwd(getwd())
list.files(dataDir)
suppressMessages(c("library(knitr)","library(limma)",
                   "library(minfi)","library(IlluminaHumanMethylation450kanno.ilmn12.hg19)",
                   "library(IlluminaHumanMethylationEPICmanifest)",
                   "library(IlluminaHumanMethylation450kmanifest)","library(RColorBrewer)",
                   "library(missMethyl)","library(minfiData)","library(Gviz)","library(DMRcate)",
                   "library(stringr)"))
library(RColorBrewer)
# annotation file
library(devtools)
#install_github("achilleasNP/IlluminaHumanMethylationEPICmanifest")
#install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")
library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
head(annEPIC)
pal = brewer.pal(8,"Dark2")

# read the sample sheet
samples <- read.metharray.sheet(dataDir, pattern = "Sample_sheet.csv")
rgSet <- read.metharray.exp(targets = samples)
rgSet
# give the samples descriptive names
samples$ID <- paste(samples$Sample_Group,samples$Sample_Name,sep=".")
sampleNames(rgSet) <- samples$ID
rgSet

# Quality Control of the dataset
detP <- detectionP(rgSet)
qcReport(rgSet, sampNames=samples$ID, sampGroups=samples$Sample_Group,
         pdf="qcReport.pdf")
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
rgSet

## remove poor quality samples from targets data
samples <- samples[keep,]
samples[,1:5]

## remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)

# Normalization
mSetRaw <- preprocessRaw(rgSet)
mSetSq <- preprocessSWAN(rgSet)

# visualise what the data looks like before and after normalisation
NormPlot <- function(normPlot){
  par(mfrow=c(1,2))
  densityPlot(rgSet, sampGroups=samples$Sample_Group,main="Raw", legend=FALSE)
  legend("top", legend = levels(factor(samples$Sample_Group)), 
         text.col=brewer.pal(8,"Dark2"))
  densityPlot(getBeta(mSetSq), sampGroups=samples$Sample_Group,
              main="SWAN.Normalized", legend=FALSE)
  legend("top", legend = levels(factor(samples$Sample_Group)), 
         text.col=brewer.pal(8,"Dark2"))
}

NormPlot()

MDSplot <-function(MDSplot){
  par(mfrow=c(1,2))
  plotMDS(getM(mSetSq), top=10000, gene.selection="common", 
          col=as.integer(factor(samples$Sample_Group)))
  legend("topright", legend=levels(factor(samples$Sample_Group)), 
         bg="white", cex=0.7)
  
  plotMDS(getM(mSetSq), top=10000, gene.selection="common",  
          col=as.integer(factor(samples$Sample_Source)))
  legend("topright", legend=levels(factor(samples$Sample_Source)),
         bg="white", cex=0.7)
}

MDSplot()


# Filtering of the data
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

## remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

## if your data includes males and females, remove probes on the sex chromosomes

keep <- !(featureNames(mSetSqFlt) %in% annEPIC$Name[annEPIC$chr %in%
                                                      c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]

## remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

## exclude cross reactive probes
library(RCurl)
nsp <- getURL("https://raw.githubusercontent.com/sirselim/illumina450k_filtering/master/48639-non-specific-probes-Illumina450k.csv")
xReactiveProbes <- read.csv(text = nsp)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]
mSetSqFlt

MDSplot.filt <-function(MDSplotFilt){
  par(mfrow=c(1,2))
  plotMDS(getM(mSetSq), top=10000, gene.selection="common", 
          col=as.integer(factor(samples$Sample_Group)))
  legend("topright", legend=levels(factor(samples$Sample_Group)), 
         bg="white", cex=0.7)
  
  plotMDS(getM(mSetSq), top=10000, gene.selection="common",  
          col=as.integer(factor(samples$Sample_Source)))
  legend("topright", legend=levels(factor(samples$Sample_Source)),
         bg="white", cex=0.7)
}
MDSplot.filt()

# Calculation of M values and beta values
## calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)

Mval.bValplot <-function(mval.bVal){
  par(mfrow=c(1,2))
  densityPlot(bVals, sampGroups=samples$Sample_Group, main="Beta values", 
              legend=FALSE, xlab="Beta values")
  legend("top", legend = levels(factor(samples$Sample_Group)), 
         text.col=brewer.pal(8,"Dark2"))
  densityPlot(mVals, sampGroups=samples$Sample_Group, main="M-values", 
              legend=FALSE, xlab="M values")
  legend("topleft", legend = levels(factor(samples$Sample_Group)), 
         text.col=brewer.pal(8,"Dark2"))
}
Mval.bValplot()


# Cell deconvolution
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("FlowSorted.Blood.EPIC")
suppressMessages(suppressWarnings(library(FlowSorted.Blood.EPIC)))

countsEPIC<-estimateCellCounts2(rgSet, compositeCellType = "Blood",
                                processMethod = "preprocessNoob",
                                probeSelect = "IDOL",
                                cellTypes = c("CD8T", "CD4T", "NK", "Bcell",
                                              "Mono", "Neu"),
                                referencePlatform =
                                  "IlluminaHumanMethylationEPIC",
                                referenceset = NULL,
                                IDOLOptimizedCpGs =IDOLOptimizedCpGs,
                                returnAll = FALSE)
par(mfrow=c(1,1))
a = countsEPIC$counts[samples$Sample_Group == "Infection",]
b = countsEPIC$counts[samples$Sample_Group == "Mock",]
boxplot(a, at=0:5*3 + 1, xlim=c(0, 18), ylim=range(a, b), xaxt="n",
        col=pal[1], main="", ylab="Cell type proportion")
boxplot(b, at=0:5*3 + 2, xaxt="n", add=TRUE, col=pal[2])
axis(1, at=0:5*3 + 1.5, labels=colnames(a), tick=TRUE)
legend("topleft", legend=c("Infection","Mock"), fill=pal)


# Testing Differential Methylation using Limma
groups <- factor(samples$Sample_Group)

## Using the Group only
design.groups <- model.matrix(~0+groups, data = samples)
colnames(design.groups) <- c(levels(groups))

fit.groups <- lmFit(mVals, design.groups)
contMatrix.groups <- makeContrasts(Infection-Mock, levels = design.groups)
fit2.groups <- contrasts.fit(fit.groups, contMatrix.groups)
fit2.groups.ebayes <- eBayes(fit2.groups)
summary(decideTests(fit2.groups.ebayes))

topTable(fit2.groups.ebayes)
# get the table of results for the group

annoEPIC <- as.data.frame(annEPIC)

annoEPICSub <- annoEPIC[match(rownames(mVals),annoEPIC$Name),
                        c(1:4,12:19,24:ncol(annoEPIC))]

DMCs.groups.ann <- topTable(fit2.groups.ebayes, 
                            num = Inf, 
                            coef = 1, 
                            genelist = annoEPICSub)
head(DMCs.groups.ann)

write.table(DMCs.groups.ann, "DMC_result_groupsInfectionvsMOCK.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = T)
write.table(mVals, "mVals_Group_BEA21P030_C19.txt", sep = "\t", quote = F, row.names = T)
write.table(bVals, "Beta_Vals_Group_BEA21P030_C19.txt", sep = "\t", quote = F, row.names = T)
save.image("EPIC_C19invitro.RData")

## Using Individuals
individual <- factor(samples$Sample_Name)
design.ind <- model.matrix(~0+individual, data = samples)
colnames(design.ind) <- c(levels(individual))
design.ind
fit.ind <- lmFit(mVals, design.ind)
contMatrix.ind.D1 <- makeContrasts(D1_Covid, D1_Mock, Diff= (D1_Covid - D1_Mock), 
                                   levels = design.ind)

fit2.ind.D1 <- contrasts.fit(fit.ind, contMatrix.ind.D1)
fit2.D1.ebayes <- eBayes(fit2.ind.D1) # will not work for single sample

# How to do the analysis
head(bVals)
anno.bVals <- cbind(bVals, annoEPICSub)
head(anno.bVals)

df.InfD1 <- anno.bVals[,c(9,10,1)]
head(df.InfD1)
df.InfD1$Unmeth = 1 - abs(df.InfD1$Infection.D1_Covid)
df.InfD1$N <- df.InfD1$Infection.D1_Covid + df.InfD1$Unmeth
df.InfD1$X <- df.InfD1$Infection.D1_Covid
df.IndD1.data <- df.InfD1[c(1,2,5,6)]
head(df.IndD1.data)

df.mockD1 <- anno.bVals[,c(9,10,5)]
head(df.mockD1)
df.mockD1$Unmeth = 1 - abs(df.mockD1$Mock.D1_Mock)
df.mockD1$N <- df.mockD1$Mock.D1_Mock + df.InfD1$Unmeth
df.mockD1$X <- df.mockD1$Mock.D1_Mock
df.mockD1.data <- df.mockD1[c(1,2,5,6)]

diff.obj <- makeBSseqData(list(df.IndD1.data, df.mockD1), c("Inf", "Mock"))
dmlTest.sm = DMLtest(diff.obj, group1 = "Inf", group2 = "Mock", smoothing = TRUE,
                     smoothing.span = 500)
dmls = callDML(dmlTest.sm, p.threshold = 1)
head(dmls)
dmrs = callDMR(dmlTest.sm, p.threshold=1)


## An approach to do the single individual analysis from mVals
head(mVals)
### select only D1 samples
d1.InfMock <- as.data.frame(mVals[,c(1,5)])
### Calculate log2 values for absolute numeric M-values
d1.InfMock$log2Inf <- -log(as.numeric(abs(d1.InfMock$Infection.D1_Covid)), 2)
d1.InfMock$log2Mock <- -log(as.numeric(abs(d1.InfMock$Mock.D1_Mock)), 2)
### Get the log2FC from the Infection - Mock data
d1.InfMock$log2FC <- d1.InfMock$log2Inf -d1.InfMock$log2Mock
### check to remove ±Inf from the data, if any
d1.InfMock <- d1.InfMock[is.finite(d1.InfMock$log2FC), ] 
head(d1.InfMock)
### Calculate the t.test
d1.InfMock$t <- (mean(d1.InfMock$log2FC))/(sd(d1.InfMock$log2FC)/sqrt(length(d1.InfMock$log2FC)))
### Calculate the p.value
d1.InfMock$p.Value <- 2*pt(-abs(d1.InfMock$t), df=length(d1.InfMock$log2FC)-1)


## Generate a line density plot to set the log2FC cut-off
### From M val
library(ggplot2)
p.dens.mVals <- ggplot(d1.InfMock, aes(x=log2FC)) + 
  geom_density()+
  coord_cartesian(xlim=c(-2, 2))
p.dens.mVals + 
  geom_vline(aes(xintercept=-0.5),color="blue", linetype="dashed", size=0.5) +
  geom_vline(aes(xintercept=0.5),color="blue", linetype="dashed", size=0.5)

### From beta value
head(bVals)
Bvals <- as.data.frame(bVals)
d1.bVals <- Bvals[c(1,5)]
d1.bVals$log2Inf <- -log(as.numeric(abs(d1.bVals$Infection.D1_Covid)), 2)
d1.bVals$log2mock <- -log(as.numeric(abs(d1.bVals$Mock.D1_Mock)), 2)
d1.bVals$log2FC <- d1.bVals$log2Inf - d1.bVals$log2mock
d1.bVals <- d1.bVals[is.finite(d1.bVals$log2FC), ] 
head(d1.bVals)

p.dens.bVals <- ggplot(d1.bVals, aes(x=log2FC)) + 
  geom_density()+
  coord_cartesian(xlim=c(-1, 1))
p.dens.bVals + 
  geom_vline(aes(xintercept=-0.25),color="blue", linetype="dashed", size=0.5) +
  geom_vline(aes(xintercept=0.25),color="blue", linetype="dashed", size=0.5)

## Select CpGs basd on the density plot cut-off, |log2FCmVals| = 2
d1.InfMock.CpGs.hypo <- d1.InfMock[d1.InfMock$log2FC < -2,]
d1.InfMock.CpGs.hyper <- d1.InfMock[d1.InfMock$log2FC > 2,]


### select only D2 samples
d2.InfMock <- as.data.frame(mVals[,c(2,6)])
### Calculate log2 values for absolute numeric M-values
d2.InfMock$log2Inf <- -log(as.numeric(abs(d2.InfMock$Infection.D2_Covid)), 2)
d2.InfMock$log2Mock <- -log(as.numeric(abs(d2.InfMock$Mock.D2_Mock)), 2)
### Get the log2FC from the Infection - Mock data
d2.InfMock$log2FC <- d2.InfMock$log2Inf -d2.InfMock$log2Mock
### check to remove ±Inf from the data, if any
d2.InfMock <- d2.InfMock[is.finite(d2.InfMock$log2FC), ] 
head(d2.InfMock)
### Calculate the t.test
d2.InfMock$t <- (mean(d2.InfMock$log2FC))/(sd(d2.InfMock$log2FC)/sqrt(length(d2.InfMock$log2FC)))
### Calculate the p.value
d2.InfMock$p.Value <- 2*pt(-abs(d2.InfMock$t), df=length(d2.InfMock$log2FC)-1)

d2.InfMock.CpGs.hypo <- d2.InfMock[d2.InfMock$log2FC < -2,]
d2.InfMock.CpGs.hyper <- d2.InfMock[d2.InfMock$log2FC > 2,]



### select only D3 samples
d3.InfMock <- as.data.frame(mVals[,c(3,7)])
### Calculate log2 values for absolute numeric M-values
d3.InfMock$log2Inf <- -log(as.numeric(abs(d3.InfMock$Infection.D3_Covid)), 2)
d3.InfMock$log2Mock <- -log(as.numeric(abs(d3.InfMock$Mock.D3_Mock)), 2)
### Get the log2FC from the Infection - Mock data
d3.InfMock$log2FC <- d3.InfMock$log2Inf -d3.InfMock$log2Mock
### check to remove ±Inf from the data, if any
d3.InfMock <- d3.InfMock[is.finite(d3.InfMock$log2FC), ] 
head(d3.InfMock)
### Calculate the t.test
d3.InfMock$t <- (mean(d3.InfMock$log2FC))/(sd(d3.InfMock$log2FC)/sqrt(length(d3.InfMock$log2FC)))
### Calculate the p.value
d3.InfMock$p.Value <- 2*pt(-abs(d3.InfMock$t), df=length(d3.InfMock$log2FC)-1)

d3.InfMock.CpGs.hypo <- d3.InfMock[d3.InfMock$log2FC < -2,]
d3.InfMock.CpGs.hyper <- d3.InfMock[d3.InfMock$log2FC > 2,]


### select only D4 samples
d4.InfMock <- as.data.frame(mVals[,c(4,8)])
### Calculate log2 values for absolute numeric M-values
d4.InfMock$log2Inf <- -log(as.numeric(abs(d4.InfMock$Infection.D4_Covid)), 2)
d4.InfMock$log2Mock <- -log(as.numeric(abs(d4.InfMock$Mock.D4_Mock)), 2)
### Get the log2FC from the Infection - Mock data
d4.InfMock$log2FC <- d4.InfMock$log2Inf -d4.InfMock$log2Mock
### check to remove ±Inf from the data, if any
d4.InfMock <- d4.InfMock[is.finite(d4.InfMock$log2FC), ] 
head(d4.InfMock)
### Calculate the t.test
d4.InfMock$t <- (mean(d4.InfMock$log2FC))/(sd(d4.InfMock$log2FC)/sqrt(length(d4.InfMock$log2FC)))
### Calculate the p.value
d4.InfMock$p.Value <- 2*pt(-abs(d4.InfMock$t), df=length(d4.InfMock$log2FC)-1)

d4.InfMock.CpGs.hypo <- d4.InfMock[d4.InfMock$log2FC < -2,]
d4.InfMock.CpGs.hyper <- d4.InfMock[d4.InfMock$log2FC > 2,]

# Venn from 4 samples, hyper & hypo separately
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
install.packages("sf" , dependencies = TRUE)
install.packages("ggVennDiagram")

library("ggVennDiagram")

Venn

d1.hyper <- read.table("d1_InfMock_CpGs_Hyper.txt", stringsAsFactors = F, 
                       header = T)
d2.hyper <- read.table("d2_InfMock_CpGs_Hyper.txt", stringsAsFactors = F, 
                       header = T)
d3.hyper <- read.table("d3_InfMock_CpGs_Hyper.txt", stringsAsFactors = F, 
                       header = T)
d4.hyper <- read.table("d4_InfMock_CpGs_Hyper.txt", stringsAsFactors = F, 
                       header = T)

vData <- list(
  D1 = sample(rownames(d1.hyper),8349),
  D2 = sample(rownames(d2.hyper), 5038),
  D3 = sample(rownames(d3.hyper), 6010),
  D4 = sample(rownames(d4.hyper), 12129)
) 
ggVennDiagram(vData, label_alpha = 0) +
  ggplot2::scale_fill_gradient(low="white",high = "pink")


if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library("ggvenn")
ggvenn(
  vData, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 6
)



d1.hypo <- read.table("d1_InfMock_CpGs_Hypo.txt", stringsAsFactors = F, 
                       header = T)
d2.hypo <- read.table("d2_InfMock_CpGs_Hypo.txt", stringsAsFactors = F, 
                       header = T)
d3.hypo <- read.table("d3_InfMock_CpGs_Hypo.txt", stringsAsFactors = F, 
                       header = T)
d4.hypo <- read.table("d4_InfMock_CpGs_Hypo.txt", stringsAsFactors = F, 
                       header = T)

vData.hypo <- list(
  D1 = sample(rownames(d1.hypo),10110),
  D2 = sample(rownames(d2.hypo), 5802),
  D3 = sample(rownames(d3.hypo), 6329),
  D4 = sample(rownames(d4.hypo), 32461)
)

ggVennDiagram(vData.hypo, label_alpha = 0)

ggvenn(
  vData.hypo, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 6
)



venn.hyper <- ggvenn(
  vData, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 6
)
venn.hypo <- ggvenn(
  vData.hypo, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 6
)

par(mfrow=c(1,2))
library(ggpubr)

ggarrange(venn.hyper, venn.hypo, labels = c("A. Hypermethylated", "B. Hypomethylated"))


d1.InfMock.anno.hyper <- read.table("d1_InfMock_anno_Hyper.txt",
                                    stringsAsFactors = F, header = T, fill = T)

d2.InfMock.anno.hyper <- read.table("d2_InfMock_anno_Hyper.txt",
                                    stringsAsFactors = F, header = T, fill = T)
d3.InfMock.anno.hyper <- read.table("d3_InfMock_anno_Hyper.txt",
                                    stringsAsFactors = F, header = T, fill = T)
d4.InfMock.anno.hyper <- read.table("d4_InfMock_anno_Hyper.txt",
                                    stringsAsFactors = F, header = T, fill = T)

d1.InfMock.anno.hypo <- read.table("d1_InfMock_anno_Hypo.txt",
                                    stringsAsFactors = F, header = T, fill = T)

d2.InfMock.anno.hypo <- read.table("d2_InfMock_anno_Hypo.txt",
                                    stringsAsFactors = F, header = T, fill = T)
d3.InfMock.anno.hypo <- read.table("d3_InfMock_anno_Hypo.txt",
                                    stringsAsFactors = F, header = T, fill = T)
d4.InfMock.anno.hypo <- read.table("d4_InfMock_anno_Hypo.txt",
                                    stringsAsFactors = F, header = T, fill = T)

vData.anno.hyper <- list(
  D1 = sample(d1.InfMock.anno.hyper$GencodeBasicV12_NAME, 8346),
  D2 = sample(d2.InfMock.anno.hyper$GencodeBasicV12_NAME, 5036),
  D3 = sample(d3.InfMock.anno.hyper$GencodeBasicV12_NAME, 6007),
  D4 = sample(d4.InfMock.anno.hyper$GencodeBasicV12_NAME, 12124)
)
vData.anno.hypo <- list(
  D1 = sample(d1.InfMock.anno.hypo$GencodeBasicV12_NAME, 10102),
  D2 = sample(d2.InfMock.anno.hypo$GencodeBasicV12_NAME, 5793),
  D3 = sample(d3.InfMock.anno.hypo$GencodeBasicV12_NAME, 6327),
  D4 = sample(d4.InfMock.anno.hypo$GencodeBasicV12_NAME, 32450)
)


venn.hyper <- ggvenn(
  vData.anno.hyper, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 6
)
venn.hypo <- ggvenn(
  vData.anno.hypo, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 6
)


ggarrange(venn.hyper, venn.hypo, labels = c("A. Hypermethylated", "B. Hypomethylated"))



### arrange genes
head(d1.InfMock.anno.hyper)

# reading single lines GeneSymbols and put them into venn analysis
library("ggVennDiagram")
library("ggpubr")
d1.anno.hyper.sort <- read.table("d1_anno_hyper_sort.txt", stringsAsFactors = F)
d2.anno.hyper.sort <- read.table("d2_anno_hyper_sort.txt", stringsAsFactors = F)
d3.anno.hyper.sort <- read.table("d3_anno_hyper_sort.txt", stringsAsFactors = F)
d4.anno.hyper.sort <- read.table("d4_anno_hyper_sort.txt", stringsAsFactors = F)

d1.anno.hypo.sort <- read.table("d1_anno_hypo_sort.txt", stringsAsFactors = F)
d2.anno.hypo.sort <- read.table("d2_anno_hypo_sort.txt", stringsAsFactors = F)
d3.anno.hypo.sort <- read.table("d3_anno_hypo_sort.txt", stringsAsFactors = F)
d4.anno.hypo.sort <- read.table("d4_anno_hypo_sort.txt", stringsAsFactors = F)

vData.gene.anno.hyper <- list(
  D1 = sample(d1.anno.hyper.sort$V1, 9058),
  D2 = sample(d2.anno.hyper.sort$V1, 5983),
  D3 = sample(d3.anno.hyper.sort$V1, 7052),
  D4 = sample(d4.anno.hyper.sort$V1, 14574)
)

vData.gene.anno.hypo <- list(
  D1 = sample(d1.anno.hypo.sort$V1, 11231),
  D2 = sample(d2.anno.hypo.sort$V1, 6785),
  D3 = sample(d3.anno.hypo.sort$V1, 7226),
  D4 = sample(d4.anno.hypo.sort$V1, 37082)
)
hyper.gene.anno <- ggVennDiagram(vData.gene.anno.hyper, label_alpha = 0)
hypo.gene.anno <- ggVennDiagram(vData.gene.anno.hypo, label_alpha = 0)

library(ggpubr)
ggarrange(hyper.gene.anno, hypo.gene.anno, labels = c("A. Hypermethylated", "B. Hypomethylated"),
          common.legend = TRUE, legend = "right",
          font.label = list(size = 14, color = "black", family = "Times New Roman"))


venn.gene.hyper <- ggvenn(
  vData.gene.anno.hyper, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 6
)
venn.gene.hypo <- ggvenn(
  vData.gene.anno.hypo, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 6
)


ggarrange(venn.gene.hyper, venn.gene.hypo, labels = c("A. Hypermethylated", "B. Hypomethylated"))
