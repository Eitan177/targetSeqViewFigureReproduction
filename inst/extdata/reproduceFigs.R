#The purpose of this code is to reproduce the figures from our paper accompanying
#the release of the targetSeqView package (Halper-Stromberg, 2013).

library(targetSeqViewFigureReproduction)
library(grid)

## perform realignment for Figs 1 and 3
path <- system.file("extdata", package="targetSeqViewFigureReproduction")
filename=file.path(path,'Figure1and3regions.txt')

### If you have access to 2 cpus, set nodes to 2
nodes=1
registerDoMC(nodes)

## some alignment defaults are changed for the second event, notably the allowd mismatches and the gap opening penalty.
## The gap opening penalty is increased so that mismatching bases are aligned as mismatches rather than gaps, and the
## allowed mismatches is increased so that these mismatches are shown in the display. With default settings the second event will still
## demonstrate markedly better alignments supporting the deletion. However, we would like to show the mismatches at the edges since these bases
## represent the action of the enzyme terminal deoxynucleotidyl transferase (TdT). Extra bases at the edges are expected for this junction
## based upon a context clue: this junction juxtaposes two coding V(D)J elements (in this case two V elements), a process that involves the action
## of TdT.

retFig1and3=ViewAndScore(filename=filename,bamFilePath=path,
estimateIndelRate=FALSE,estimateMmRate=FALSE,getReadLength=FALSE,
build='hg19',verbose=TRUE,gapOpeningArg=c(-4,-100),allowedMM=c(6,25),refexpansion=c(400,0))
## The likelihood scores for the these two events
print(retFig1and3[[2]])

    ## View the alignments
for(ii in seq_along(retFig1and3[[1]])){
## There are >120 readpairs in the Fig1 regions, we set the covMaxDisplay to 50 so that our display does not look squished
p1=formatPlot(retFig1and3[[1]][[ii]][[1]][[2]],title='Alignment supporting a structural variant',mismatchcolor='red',covMaxDisplay=50)
p2=formatPlot(retFig1and3[[1]][[ii]][[2]][[2]],title='Alignment supporting no structural variant',mismatchcolor='red',covMaxDisplay=50)
p3=formatPlot(retFig1and3[[1]][[ii]][[3]][[2]],title='Alignment supporting no structural variant',mismatchcolor='red',covMaxDisplay=50)
dev.new()
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
print(p1,vp = viewport(layout.pos.row = 1, layout.pos.col=1))
pushViewport(viewport(layout = grid.layout(2, 2)))
print(p2,vp = viewport(layout.pos.row = 2, layout.pos.col=1))
print(p3,vp = viewport(layout.pos.row = 2, layout.pos.col=2))

}

## Reproduce the ROC

## To do full realignment call reproduceRocLong. For full realignment, speed can be greatly increased if multiple cpus are available.
## To take advantage of multiple cpus, pass the number of available cpus into the function reproduceRocLong using the 'nodes' argument.
## The reproduceRocShort function loads an R object containing the realignments precomputed.

reproduceRocShort()

## Reproduce Repetitive DNA bar chart
forRepeatFig=reproduceMappability()

g1=ggplot(forRepeatFig,aes(x=map,fill=val))+geom_bar(colour="black")+ scale_fill_manual(values=c('#FFFF00','red','green',"#0099FF","#0066FF","#66CCFF"))+
ylab('Base Pair Count')+xlab('Alignability of Baits and Surrounding Sequence')+facet_wrap(~ Type,scales="free_y")

forRepeatFig2=forRepeatFig[-which(forRepeatFig$map==1),]
newname=rep('',nrow(forRepeatFig2))
newname[which(forRepeatFig2$Type=='All Baits')]="All Baits < 100% Alignable"
newname[which(forRepeatFig2$Type=='Reported Events')]="Reported Events < 100% Alignable"
forRepeatFig2[,4]=newname
g2=ggplot(forRepeatFig2,aes(x=map,fill=val))+geom_bar(colour="black")+ scale_fill_manual(values=c('#FFFF00','red','green',"#0099FF","#0066FF",
          "#66CCFF"))+ylab('Base Pair Count')+xlab('Alignability of Baits and Surrounding Sequence')+facet_wrap(~ Type,scales="free_y")

dev.new()
print(g1)
dev.new()
print(g2)


