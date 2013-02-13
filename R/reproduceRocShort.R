reproduceRocShort <-
function(){
    path <- system.file("extdata", package="targetSeqViewFigureReproduction")
    load(file.path(path,"realignments.rda"))
    library(ROCR)
    reproduceRoc(ret)
}
