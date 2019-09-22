setwd("/users/jleibowitz/Dropbox (Partners Healthcare)/CEL")
library(oligo)
library(annotate)

celFiles<-list.celfiles()
rawData<-read.celfiles(celFiles)

probes<-probeNames(rawData)
uniqueprobeID<-unique(probes)

bgData1<-backgroundCorrect(rawData)
normData <- normalize(bgData1)
int<-intensity(rawData)
colnames(norm)<-colnames(int)

test<-rma(rawData)
head(exprs(test))

annot2<-read.table("MoEx-1_0-st-v1.na33.1.mm9.transcript.csv",sep=",",header=FALSE)





ID<-featureNames(rawData)

######## start here
celFiles <- list.celfiles()
affyRaw <- oligo::read.celfiles(celFiles)
ph<-affyRaw@phenoData
feat <- affyRaw@featureData
test<-feat@data
ID<-featureNames(affyRaw)

probeID<-probeNames(affyRaw)
uniqueprobeID<-probesetNames(affyRaw)
uniqueprobeID<-unique(probeID)


pmSeq<-pmSequence(affyRaw)
data(pmSequence)

platDesign<-getPlatformDesign(affyRaw)
probeInfo<-getProbeInfo(affyRaw,target="probeset")

#expr = exprs(affyRaw)
#int = intensity(affyRaw)
int<-pm(affyRaw)

bgData1<-oligo::backgroundCorrect(int)

norm<-oligo::normalize(bgData1)
colnames(norm)<-colnames(int)

#sum<-summarize(norm,probes=probeInfo$man_fsetid)

sum<-matrix(NA,nrow=length(uniqueprobeID),ncol=ncol(int))
rownames(sum)<-uniqueprobeID
for (i in 1:length(uniqueprobeID)){
  temp<-norm[grep(uniqueprobeID[i],probeID),]
  temp<-t(temp)
  medtemp<-medpolish(temp)
  subResiduals<-temp-medtemp$residuals
  mean<-rowMeans(subResiduals)
  sum[i,]<-t(mean)
}
colnames(sum)<-colnames(norm)



annot<-read.table("MoEx-1_0-st-v1.na33.1.mm9.transcript.csv",sep=",",header=FALSE)
AffyIds<-rownames(sum)
genenames<-annot$V8[match(AffyIds,annot$V2)]

Final<-cbind(as.character(genenames),sum)

write.csv(Final,"final.csv")

eset <- rma(affyRaw)

head(exprs(eset))
