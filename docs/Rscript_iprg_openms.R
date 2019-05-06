##### iprg-OpenMS

library(MSstats)
?MSstats

# First, read output of OpenMS
raw.openMS <- read.csv("data/data_OpenMS/ABRF2015_OpenMS_raw.csv", stringsAsFactors=F) # the data file
head(raw.openMS)

# Set annotation file
# The output from `OpenMS` already includes `Run`, `BioReplicate`, `Condition` information. Let's check it.

unique(raw.openMS[, c('Run', 'BioReplicate', 'Condition')])

?OpenMStoMSstatsFormat

# reformating and pre-processing for OpenMS output.
input.openms <- OpenMStoMSstatsFormat(raw.openMS,
                                      removeProtein_with1Feature=TRUE)

## now 'input.openms' is ready for MSstats
head(input.openms)

length(unique(input.openms$ProteinName)) 
sum(is.na(input.openms$Intensity)) 
sum(!is.na(input.openms$Intensity) & input.openms$Intensity==0)
table(input.openms$Run)

## save the work
save(input.openms, file='data/data_OpenMS/input.openms.rda')

load(file='data/data_OpenMS/input.openms.rda')

quant.openms <- dataProcess(raw = input.openms, 
                            logTrans=2, 
                            #normalization = 'quantile',
                            summaryMethod = 'TMP', 
                            MBimpute=TRUE,
                            censoredInt='NA',
                            cutoffCensored='minFeature',
                            maxQuantileforCensored = 0.999)

dataProcessPlots(data = quant.openms, 
                 type="QCplot", 
                 width=7, height=7,
                 which.Protein = 'allonly',
                 address='data/data_OpenMS/ABRF_openms_equalMed_')

dataProcessPlots(data = quant.openms, 
                 type="Profileplot", 
                 featureName="NA",
                 width=7, height=7,
                 which.Protein = 'sp|P55249|ZRT4_YEAST',
                 address="data/data_OpenMS/ABRF_openms_equalMed_P55249_")

# save your work
save(quant.openms, file='data/data_OpenMS/quant.openms.rda')


## inference
load(file='data/data_OpenMS/quant.openms.rda')

unique(quant.skyline$ProcessedData$GROUP_ORIGINAL)

comparison1<-matrix(c(-1,1,0,0),nrow=1)
comparison2<-matrix(c(-1,0,1,0),nrow=1)
comparison3<-matrix(c(-1,0,0,1),nrow=1)
comparison4<-matrix(c(0,-1,1,0),nrow=1)
comparison5<-matrix(c(0,-1,0,1),nrow=1)
comparison6<-matrix(c(0,0,-1,1),nrow=1)
comparison<-rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6)
row.names(comparison)<-c("C2-C1","C3-C1","C4-C1","C3-C2","C4-C2","C4-C3")

comparison


test.openms <- groupComparison(contrast.matrix=comparison, data=quant.openms)
OpenMS.result <- test.openms$ComparisonResult

# save your work
save(OpenMS.result, file='data/data_OpenMS/OpenMS.result.rda')
write.csv(OpenMS.result, file='data/data_OpenMS/testResult_ABRF_openms.csv')

groupComparisonPlots(data = OpenMS.result, 
                     type = 'VolcanoPlot',
                     address = 'data/data_OpenMS/testResult_ABRF_openms_')

groupComparisonPlots(data = OpenMS.result, 
                     type = 'ComparisonPlot',
                     which.Protein = 'sp|P55249|ZRT4_YEAST',
                     address = 'data/data_OpenMS/testResult_ABRF_openms_')


