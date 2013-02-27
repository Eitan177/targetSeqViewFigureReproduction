reproduceVenn <- function(){
getTop=T;
tophits=1000
library('Vennerable')
    gasvrngs <- gasv4venn(getTop)
    vhrngs <- variationHunter4venn(getTop)
    hydrarngs <- hydra4venn(getTop)

assign('retrem',IRanges(),envir=.GlobalEnv)
subjectreturns=list()
for(ii in 1:3){
    first=c('hydrarngs','hydrarngs','vhrngs')
    second=c('gasvrngs','vhrngs','gasvrngs')
    ovs=findOverlaps(get(first[ii]),get(second[ii]))
    rngnames=seqnames(get(first[ii])[queryHits(ovs)])
    secondnames=seqnames(get(second[ii])[subjectHits(ovs)])
    keepovs=which(as.character(rngnames)==as.character(secondnames))
    ovs=ovs[keepovs]
    sovs=IRanges::split(ovs,(queryHits(ovs)-1)%/%2)
    subjectreturns[[ii]]=IRanges::lapply(sovs,function(x){ss=split(subjectHits(x),queryHits(x));firsts=IRanges(start=ss[[1]]-(ss[[1]]-1)%%2,width=2);
                                                          if(length(ss)==1){print('variationhunter junctions not found');seconds=firsts;ret=IRanges()}else{
                                                              seconds=IRanges(start=ss[[2]]-(ss[[2]]-1)%%2,width=2);#};
                                                              ff=findOverlaps(firsts,seconds);
                                                              ret=firsts[queryHits(ff)]
                                                              if(nrow(IRanges::as.matrix(ff))!=1){
                                                                  print('more than 1 loc will be removed')
                                                              }
                                                              if(length(ret)>1){
                                                                  matchret=na.omit(IRanges::match(get("retrem",envir=.GlobalEnv),ret))
                                                                  if(length(matchret)>0){
                                                                      print('omitting duplicated hit')
                                                                      ret = ret[-matchret]
                                                                  }
                                                                  ret=ret[which.max(as.numeric(elementMetadata(vhrngs[start(ret)])[,2]))]
                                                                  assign("retrem",c(retrem,ret),envir=.GlobalEnv)

                                                              }};
                                                          ret})
    tt=table(((Rle(queryHits(ovs)))@values-1)%/%2);names(subjectreturns[[ii]])=paste(names(tt),tt,sep='.')
    assign("retrem",IRanges(),envir=.GlobalEnv)
}

varigasv=names(subjectreturns[[3]])[which(unlist(lapply(subjectreturns[[3]],length))>0)]
varigasv=varigasv[as.numeric(gsub('.+\\.','',varigasv))>1]
varigasv=gsub('\\..+','',varigasv)
hydragasv=names(subjectreturns[[1]])[which(unlist(lapply(subjectreturns[[1]],length))>0)]
hydragasv=hydragasv[as.numeric(gsub('.+\\.','',hydragasv))>1]
hydragasv=gsub('\\..+','',hydragasv)
hydravh=names(subjectreturns[[2]])[which(unlist(lapply(subjectreturns[[2]],length))>0)]
hydravh=hydravh[as.numeric(gsub('.+\\.','',hydravh))>1]
hydravh=gsub('\\..+','',hydravh)


all3=hydravh[na.omit(match(hydragasv,hydravh))]

hydravh2=hydravh[-which(!is.na(match(hydravh,all3)))]
hydragasv2=hydragasv[-which(!is.na(match(hydragasv,all3)))]


all3n=length(all3)
hgn=length(hydragasv2)
hvn=length(hydravh2)
hn=(length(hydrarngs)/2)-all3n-hgn-hvn
gvn=max(length(varigasv)-all3n,0)
gn=max((length(gasvrngs)/2)-all3n-hgn-gvn,0)
vn=max((length(vhrngs)/2)-all3n-hvn-gvn,0)


hs=c(make.unique(rep('all3',all3n)),make.unique(rep('h',hn)),make.unique(rep('hg',hgn)),make.unique(rep('hv',hvn)))
gs=c(make.unique(rep('all3',all3n)),make.unique(rep('g',gn)),make.unique(rep('hg',hgn)),make.unique(rep('gv',gvn)))
vs=c(make.unique(rep('all3',all3n)),make.unique(rep('v',vn)),make.unique(rep('gv',gvn)),make.unique(rep('hv',hvn)))


C3 <- Venn(list(hs,vs,gs),SetNames=c('HYDRA','VariationHunter','GASV'))
plot(C3, type = "circles")
}
