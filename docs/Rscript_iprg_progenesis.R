##### iprg-progenesis

library(MSstats)
?MSstats

# First, read output of Progenesis
raw.progenesis <- read.csv("data/data_Progenesis/ABRF2015_Progenesis_raw.csv", stringsAsFactors=F) # the data file
head(raw.progenesis)

## Read in annotation including condition and biological replicates: ABRF2015_Progenesis_annotation.csv
annot.progenesis <- read.csv("data/data_Progenesis/ABRF2015_Progenesis_annotation.csv", header = TRUE)
annot.progenesis

?ProgenesistoMSstatsFormat

# reformating and pre-processing for Progenesis output.
input.progenesis <- ProgenesistoMSstatsFormat(raw.progenesis, 
                                              annotation=annot.progenesis,
                                              removeProtein_with1Peptide=TRUE)

## now 'input.progenesis' is ready for MSstats
head(input.progenesis)

#### Preliminary check for preprocessed data

length(unique(input.progenesis$ProteinName)) 
sum(is.na(input.progenesis$Intensity)) 
sum(!is.na(input.progenesis$Intensity) & input.progenesis$Intensity==0)
table(input.progenesis$Run)

# save the work
save(input.progenesis, file='data/data_Progenesis/input.progenesis.rda')


load(file='data/data_Progenesis/input.progenesis.rda')

head(input.progenesis)
sum(is.na(input.progenesis$Intensity)) 
sum(!is.na(input.progenesis$Intensity) & input.progenesis$Intensity==0)

# **Note!** Progenesis output has only `0`. **censoredInt='0'** should be used for Progenesis output.

quant.progenesis <- dataProcess(raw = input.progenesis, 
                                logTrans=2, 
                                #normalization = 'quantile',
                                summaryMethod = 'TMP', 
                                MBimpute=TRUE,
                                censoredInt='0',
                                cutoffCensored='minFeature',
                                maxQuantileforCensored = 0.999)


dataProcessPlots(data = quant.progenesis, 
                 type="Profileplot", 
                 featureName="NA",
                 width=7, height=7,
                 which.Protein = 'sp|P55249|ZRT4_YEAST',
                 address="data/data_Progenesis/ABRF_progenesis_equalMed_P55249_")

# save your work
save(quant.progenesis, file='data/data_Progenesis/quant.progenesis.rda')


## inference
load(file='data/data_Progenesis/quant.progenesis.rda')

unique(quant.progenesis$ProcessedData$GROUP_ORIGINAL)

comparison1<-matrix(c(-1,1,0,0),nrow=1)
comparison2<-matrix(c(-1,0,1,0),nrow=1)
comparison3<-matrix(c(-1,0,0,1),nrow=1)
comparison4<-matrix(c(0,-1,1,0),nrow=1)
comparison5<-matrix(c(0,-1,0,1),nrow=1)
comparison6<-matrix(c(0,0,-1,1),nrow=1)
comparison<-rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6)
row.names(comparison)<-c("C2-C1","C3-C1","C4-C1","C3-C2","C4-C2","C4-C3")

comparison


test.progenesis <- groupComparison(contrast.matrix=comparison, data=quant.progenesis)
Progenesis.result <- test.progenesis$ComparisonResult

# save your work
save(Progenesis.result, file='data/data_Progenesis/Progenesis.result.rda')
write.csv(Progenesis.result, file='data/data_Progenesis/testResult_ABRF_progenesis.csv')

groupComparisonPlots(data = Progenesis.result, 
                     type = 'VolcanoPlot',
                     address = 'data/data_Progenesis/testResult_ABRF_progenesis_')

groupComparisonPlots(data = Progenesis.result, 
                     type = 'ComparisonPlot',
                     which.Protein = 'P55249',
                     address = 'data/data_Progenesis/testResult_ABRF_progenesis_')

groupComparisonPlots(data = Progenesis.result, 
                     type = 'ComparisonPlot',
                     which.Protein = 'sp|P55249|ZRT4_YEAST',
                     address = 'data/data_Progenesis/testResult_ABRF_progenesis_')
