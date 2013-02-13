reproduceMappability <-
function(){
path <- system.file("extdata", package="targetSeqViewFigureReproduction")
setwd(path)

sureselect=read.delim('sureselectformat.bed',header=FALSE)

rr=read.delim(file='repeatmaskerbaits.txt',header=F)
segd=read.delim(file='genomicSuperDups.txt',header=F)
repm=GRanges(seqnames=rr[,5],ranges=IRanges(start=as.numeric(as.character(rr[,6])),end=as.numeric(as.character(gsub('[^0-9]','',rr[,7])))),elementMetadata=as.character(rr[,11]))
segdGR=GRanges(seqnames=as.character(segd[,2]),ranges=IRanges(start=as.numeric(as.character(segd[,3])),end=as.numeric(as.character(segd[,4]))))
baits2=cbind(as.character(sureselect[,1]),sureselect[,2]+1e3,sureselect[,3]-1e3)
baits=cbind(as.character(sureselect[,1]),sureselect[,2]+1e2,sureselect[,3]-1e2)
baitR=IRanges(start=as.numeric(baits[,2]),end=as.numeric(baits[,3]),names=baits[,1])
baitR2=IRanges(start=as.numeric(baits2[,2]),end=as.numeric(baits2[,3]),names=baits2[,1])
baitGR=GRanges(seqnames=names(baitR),IRanges(baitR))
baitGR=resize(baitGR,fix="center",width=2000)
baitGR2=GRanges(seqnames=names(baitR2),IRanges(baitR2))


## first do calculation for all baits
allbaits=processbaits(baitR,baitGR,repm,segdGR)


## read in caught baits ##
path <- system.file("extdata", package="targetSeqViewFigureReproduction")
setwd(path)

ir=read.delim(file='irregVDJ.txt',header=F);
nr=read.delim(file='regVDJ.txt',header=F)
irGR=GRanges(seqnames=paste('chr',gsub('hs','',ir[,2]),sep=''),IRanges(start=ir[,3],end=ir[,4]))
nrGR=GRanges(seqnames=paste('chr',gsub('hs','',nr[,2]),sep=''),IRanges(start=nr[,3],end=nr[,4]))
###
mare=IRanges::match(nrGR,baitGR2)
mare[is.na(mare)]=nearest(nrGR,baitGR2)[is.na(mare)]
mareGR=baitGR[mare]

## this needs to be slightly tailored to accomodate chr18 and chr8
mair=IRanges::match(irGR,baitGR2)
mair[which(is.na(mair))][1:2]=subjectHits(findOverlaps(irGR[is.na(mair)][1:2],baitGR2,maxgap=100))
mairGR=c(baitGR[na.omit(mair)],irGR[which(is.na(mair))])
mairGR=resize(mairGR,fix="center",width=2000)


## do calculation for only baits involved in events
baitRcaptured=IRanges::sort(c(mareGR,mairGR))
names(baitRcaptured)=seqnames(baitRcaptured)
baitGRcaptured=baitRcaptured
capturedevents=processbaits(baitRcaptured,baitGRcaptured,repm,segdGR)

newrept=rbind(allbaits,capturedevents)
newrept2=cbind(newrept,c(rep('All Baits',nrow(allbaits)),rep( 'Reported Events',nrow(capturedevents))))
colnames(newrept2)[4]="Type"

return(newrept2)
}
