targetSeqView4Roc <-
function(ret){
    likelihoodScores=ret[[2]]
    ## the first 38 scores are from events that failed validation
    nonvalids=likelihoodScores[1:38]
    ## the rest are PCR validated or canonical V(D)J recombination
    valids=likelihoodScores[39:length(likelihoodScores)]

    labels=c(rep(0,length(nonvalids)),rep(1,length(valids)))
    prag=prediction(likelihoodScores, labels)
    pred=performance(prag,"tpr", "fpr")
    return(pred)
}
