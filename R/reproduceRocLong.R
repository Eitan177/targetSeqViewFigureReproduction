reproduceRocLong <-
function(nodes=1){
    library(ROCR)
    registerDoMC(nodes)

    path=system.file("extdata", package="targetSeqViewFigureReproduction")
    filename='SVJunctionsInRocPlot.txt'
    normalBam=file.path(path,'1320KB0001.bam')
    bamFilePath=path
    ## do full realignment
    ret=ViewAndScore(filename=filename,normalBam=normalBam,bamFilePath=bamFilePath,verbose=TRUE)
    reproduceRoc(ret)

}
