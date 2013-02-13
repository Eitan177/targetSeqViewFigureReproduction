reproduceRoc <-
function(ret){

    ## get the likelihood scores from our realignments
    pred <- targetSeqView4Roc(ret)
    ## get scores from gasv results
    gprag <- gasv4Roc()
    ## gets scores from variationHunter results
    vhprag <- variationHunter4Roc()
    ## get scores from hydra results
    hprag <- hydra4Roc()


    ## I have a score for all 26 tp junctions, the others return scores for 25 or less tp junctions.
    xnumMine=xnumHydra=38;ynumMine=26
    xnumGasv=36; ynumGasv=ynumHydra=25
    xnumVH=21; ynumVH=12


    ## format the dataFrame for ggplot
    x=c(pred@x.values[[1]]*xnumMine,gprag@x.values[[1]]*xnumGasv,hprag@x.values[[1]]*xnumHydra,vhprag@x.values[[1]]*xnumVH )
    y=c(pred@y.values[[1]]*ynumMine,gprag@y.values[[1]]*ynumGasv,hprag@y.values[[1]]*ynumHydra,vhprag@y.values[[1]]*ynumVH )
    meth=c(rep('SeqView',length(pred@x.values[[1]])),rep('GASV',length(gprag@x.values[[1]])),rep('HYDRA',length(hprag@x.values[[1]])),
    rep('VariationHunter',length(vhprag@x.values[[1]])))
    forRoc=data.frame(x,y,meth)
    colnames(forRoc)[3]='group'

    groc=ggplot(forRoc,aes(x=x,y=y,colour=group))+geom_line(cex=2)+xlab('Failed Validation')+ylab('Validated and Canonical V(D)J')+
        scale_x_continuous(breaks = seq(5,35,5), labels = seq(5,35,5))+
            scale_y_continuous(breaks = seq(2,25,2), labels = seq(2,25,2))

    ## plot the roc
    dev.new()
    print(groc)

}
