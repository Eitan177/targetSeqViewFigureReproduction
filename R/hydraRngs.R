hydraRngs <-
function(getTop=FALSE){
    events=list()
    cutoff=2
    path <- system.file("extdata", package="targetSeqViewFigureReproduction")
    setwd(path)

    ll=list.files(pattern='1320KB00[0-9][0-9].breaks.final$')

    for(ii in seq_along(ll)){
        rr=read.delim(file=ll[ii],header=F)
        totest=rr[which(rr[,19]>=cutoff),]
        totest=apply(totest, 2, as.numeric)
        if(length(which(is.na(totest[,1])))>0) totest=totest[-which(is.na(totest[,1])),]
        print(dim(totest))
        events[[ii]]=cbind(gsub('.breaks.final','MultipleAlnsort.bam',ll[ii]),totest[,c(1:6,19)])
        events[[ii]][,c(3,6)]=as.numeric(events[[ii]][,c(3,6)])-500
        events[[ii]][,c(4,7)]=as.numeric(events[[ii]][,c(4,7)])+500
        events[[ii]][,c(2,5)]=paste('chr',events[[ii]][,c(2,5)],sep='')
        colnames(events[[ii]])[1:8]=c('Sample','Chr1','Start1','End1','Chr2','Start2','End2','score')
    }

    hydraevents=vector()
    for(ii in seq_along(ll)){
        hydraevents=rbind(hydraevents,events[[ii]])
    }

    hydrarngs=GRanges(seqnames=paste(rep(hydraevents[,'Sample'],each=2),gsub('chr','',as.vector(rbind(hydraevents[,'Chr1'],
                                                                hydraevents[,'Chr2']))),sep='.'),
    ranges=IRanges(start=as.numeric(rbind(hydraevents[,'Start1'],hydraevents[,'Start2'])),
    end=as.numeric(rbind(hydraevents[,'End1'],hydraevents[,'End2']))))
    elementMetadata(hydrarngs)[,1]='U';
    elementMetadata(hydrarngs)[,2]=rep(hydraevents[,'score'],each=2)


    if(getTop){
        hydrarngs=hydrarngs[order(as.numeric(elementMetadata(hydrarngs)[,2]),decreasing=T)][1:tophits]
    }

    return(hydrarngs)
}
