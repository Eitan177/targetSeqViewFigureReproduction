processbaits <-
function(baitR,baitGR,repm,segdGR){
baitw=list()
for(ii in 1:length(baitR)){
prebaitw=(system(paste('./bigWigToWig wgEncodeCrgMapabilityAlign100mer.bigWig stdout -chrom=',names(baitR)[ii]," -start=",start(baitR)[ii]," -end=",end(baitR)[ii],sep=''),intern=T))
prebaitw=prebaitw[-grep('^#',prebaitw)]
prebaits=strsplit(prebaitw,'\t')
rb2=unlist(prebaits)
rb3=cbind(rb2[seq(1,length(rb2),4)],rb2[seq(2,length(rb2),4)],rb2[seq(3,length(rb2),4)],rb2[seq(4,length(rb2),4)])
ba=(unlist(lapply(prebaits,function(x){z=GRanges(seqnames=x[1],IRanges(start=as.numeric(x[2]),end=as.numeric(x[3])-1));elementMetadata(z)=x[4];z})))
ba=do.call(c,ba)
baitw[[ii]]=ba
cat(ii)
}


for(ii in 1:length(baitw)){
  baitw[[ii]]=restrict(baitw[[ii]],start=start(baitR)[ii],end=end(baitR)[ii])
  cat(ii)
}

mapabil=vector()
for(ii in 1:length(baitw)){
  mapabil[ii]=sum(width(baitw[[ii]])*floor(as.numeric(elementMetadata(baitw[[ii]])[,1])))/sum(width(baitw[[ii]]))
}

mapabil=round(mapabil,2)

allrepl=list();allreplong=list()
repclasses= unique(elementMetadata(repm)[,1])
  for(ii in 1:(length(repclasses)+1)){
    if(ii <=length(repclasses)){
    repclass=repclasses[ii]
    reps=repm[elementMetadata(repm)[,1] %in% repclass]
  }else{
    reps=segdGR
  }

ff=findOverlaps(baitGR,reps)
ss=split(subjectHits(ff),queryHits(ff))
s2=mapply(function(x,y){sum(width(reduce(restrict(reps[x],start=start(baitGR[as.numeric(y)]),end=end(baitGR[as.numeric(y)])))))},ss,as.list(names(ss)))
s3=mapply(function(x,y){sh=shift(reduce(restrict(reps[x],start=start(baitGR[as.numeric(y)]),end=end(baitGR[as.numeric(y)]))),-(start(baitGR[as.numeric(y)])-1));
   cb=coverage(sh);nn=as.numeric(IRanges::as.list(cb)[[which(as.character(names(cb))==as.character(seqnames(sh)))]]);if(length(nn)<width(baitGR[as.numeric(y)])){nn[length(nn):width(baitGR[as.numeric(y)])]=0};nn},
   ss,as.list(names(ss)))
filv=rep(0,length(baitGR))
filvtotal=matrix(0,ncol=length(baitGR),nrow=width(baitGR)[1])
if(length(ff)>0){
filv[as.numeric(names(s2))]=(s2/(width(baitGR)[as.numeric(names(s2))]))
filvtotal[,as.numeric(colnames(s3))]=s3
}
allrepl[[ii]]=filv
allreplong[[ii]]=filvtotal

    cat(ii)
  }
ra=sapply(allrepl,function(x){f=sum(x[mapabil<.8])/(length(which(mapabil<.8)));se=sum(x)/(length(mapabil));c(f/se,f,se)})

cc=cbind(c(repclasses,'segdup'),round(t(ra),3))
cc[,2]=as.numeric(cc[,2]);cc[,4]=as.numeric(cc[,4])*100;cc[,3]=as.numeric(cc[,3])*100
fp1=cc[which(cc[,2]>1 & cc[,3]>2),]

maprep=rep(mapabil,each=2000)

rept=vector()
leg=c('SegDup','/L1','/Erv1')
hh=1
for(ii in c(31,1,6)){
tomelt=allreplong[[ii]]
colnames(tomelt)=1:ncol(tomelt)
repelement=reshape::melt.data.frame(as.data.frame(tomelt))
repelement[repelement[,2]==1,2]=leg[hh]
repelement[repelement[,2]=='0',2]=''
repe=cbind(maprep,repelement)
#repe=repe[,c(1,3,4)]
colnames(repe)=c('map','sample','val')
if(hh==1) rept=repe else rept[,3]=paste(rept[,3],repe[,3],sep='')
hh=hh+1
}
rept[,3]=gsub('^/','',rept[,3])
rept[rept[,3]=='',3]='0/3'
return(rept)
}
