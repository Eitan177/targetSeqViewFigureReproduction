gasvRngs <-
function(getTop=FALSE) {
    gevents=list()
    cutoff=2

    path <- system.file("extdata", package="targetSeqViewFigureReproduction")
    setwd(path)

    clustf=clustf=list.files(pattern="BAMToGASV_AMBIG.gasv.combined.in.clusters")
    clustf=grep('1320KB',clustf,value=T)


    for(ii in seq_along(clustf)){
    clust=read.delim(file=clustf[ii],header=T)
    totestclust=clust[which(clust[,2]>=cutoff),]
    gevents[[ii]]=cbind(paste(grep('1320KB',strsplit(clustf[ii],'\\.')[[1]],value=T),'MultipleAlnsort.bam',sep=''),
           totestclust[,6],as.numeric(gsub(',.+','',totestclust[,8]))-550,as.numeric(gsub(',.+','',totestclust[,8]))+550,totestclust[,7],
           as.numeric(gsub(',.+','',gsub('^[0-9]+, ','',totestclust[,8])))-550,as.numeric(gsub(',.+','',gsub('^[0-9]+, ','',totestclust[,8])))+550,totestclust[,2])
    gevents[[ii]][,c(2,5)]=paste('chr',gevents[[ii]][,c(2,5)],sep='')
    colnames(gevents[[ii]])[1:8]=c('Sample','Chr1','Start1','End1','Chr2','Start2','End2','score')
}

    gasvevents=vector()
    for(ii in seq_along(clustf)){
        gasvevents=rbind(gasvevents,gevents[[ii]])
    }


    gasvrngs=GRanges(seqnames=paste(
                     rep(gasvevents[,'Sample'],each=2),gsub('chr','',as.vector(rbind(gasvevents[,'Chr1'],
                                               gasvevents[,'Chr2']))),sep='.'),
    ranges=IRanges(start=as.numeric(rbind(gasvevents[,'Start1'],gasvevents[,'Start2'])),
    end=as.numeric(rbind(gasvevents[,'End1'],gasvevents[,'End2']))))
    elementMetadata(gasvrngs)[,1]='U';
    elementMetadata(gasvrngs)[,2]=rep(gasvevents[,'score'],each=2)

    g1s=gasvrngs[seq(1,length(gasvrngs),2)]
    g2s=gasvrngs[seq(2,length(gasvrngs),2)]
    tocloseSoremove=(which(width(pgap(ranges(g1s),ranges(g2s)))<1e3))
    closeremove=tocloseSoremove*2
    closerem=sort(c(closeremove-1,closeremove))
    gasvrngs=gasvrngs[-closerem]

    if(getTop){
    gasvrngs=gasvrngs[order(as.numeric(elementMetadata(gasvrngs)[,2]),decreasing=T)][1:tophits]
}

    return(gasvrngs)
}
