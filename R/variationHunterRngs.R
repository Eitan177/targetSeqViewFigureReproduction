variationHunterRngs <-
function(getTop=FALSE){
vevents=list()
cutoff=2
path <- system.file("extdata", package="targetSeqViewFigureReproduction")
setwd(path)

clustf=list.files(pattern="New1320KB00NoNs.+rFast_DIVET.vh.SV$")
for(ii in seq_along(clustf)){
    totestclust=read.delim(file=clustf[ii],header=F,skip=4,sep=" ")

    totestclust=totestclust[-nrow(totestclust),]
    vevents[[ii]]=cbind(paste(gsub('MrFast.+','',gsub('NoNs','',gsub('^New','',clustf[ii]))),'MultipleAlnsort.bam',sep=''),
           gsub('^chro:','',totestclust[,5]),as.numeric(gsub('^Inside_Start:','',totestclust[,1]))-500,
           as.numeric(gsub('^Inside_End:','',totestclust[,2]))+500,gsub('^chro:','',totestclust[,5]),
           as.numeric(gsub('^OutSide_Start:','',totestclust[,3]))-500,
           as.numeric(gsub('^Oustide_End:','',totestclust[,4]))+500,
           as.numeric(gsub('^sumProb:','',totestclust[,9])))
    colnames(vevents[[ii]])[1:8]=c('Sample','Chr1','Start1','End1','Chr2','Start2','End2','score')
    vevents[[ii]][abs(as.numeric(vevents[[ii]][,8]))==Inf,8]=as.numeric(gsub('^sup:','',totestclust[abs(as.numeric(vevents[[ii]][,8]))==Inf,7]))
}

vhevents=vector()
for(ii in seq_along(clustf)){
    vhevents=rbind(vhevents,vevents[[ii]])
}




vhrngs=GRanges(seqnames=paste(
               rep(vhevents[,'Sample'],each=2),gsub('chr','',as.vector(rbind(vhevents[,'Chr1'],
                                       vhevents[,'Chr2']))),sep='.'),
ranges=IRanges(start=as.numeric(rbind(vhevents[,'Start1'],vhevents[,'Start2'])),
end=as.numeric(rbind(vhevents[,'End1'],vhevents[,'End2']))))
elementMetadata(vhrngs)[,1]='U';
elementMetadata(vhrngs)[,2]=rep(vhevents[,'score'],each=2)

if(getTop){
    vhrngs=vhrngs[order(as.numeric(elementMetadata(vhrngs)[,2]),decreasing=T)][1:tophits]
}

return(vhrngs)

}
