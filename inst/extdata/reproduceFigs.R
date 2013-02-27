#The purpose of this code is to reproduce the figures from our paper accompanying
#the release of the targetSeqView package (Halper-Stromberg, 2013).

library(targetSeqViewFigureReproduction)
library(grid)


path <- system.file("extdata", package="targetSeqViewFigureReproduction")


### If you have access to 2 cpus, set nodes to 2
nodes=1
registerDoMC(nodes)
filename=file.path(path,'FigureRegs.txt')
substitutionMatSplits=nucleotideSubstitutionMatrix(match = 5, mismatch = -3)
substitutionMat=nucleotideSubstitutionMatrix(match = 1, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)]
retFigs=ViewAndScore(filename=filename,bamFilePath=path,estimateIndelRate=FALSE,estimateMmRate=FALSE,getReadLength=FALSE,verbose=T,
gapOpeningArg=c(-4,-100,-150,-4),substitutionMat=list(substitutionMat,substitutionMat,substitutionMatSplits,substitutionMat),
gapExtensionArg=c(-1,-1,0,-1),allowedMM=c(6,25,15,6),refexpansion=c(400,0,400,400),
findSplitReads=T,dedup=F)

### This will plot (in order) Figure 3, supplemental Figure 3, Figure 4B, supplemental Figure 1
for(ii in 1:4){
p1=formatPlot(retFigs[[1]][[ii]][[1]][[2]],title='Alignment supporting a structural variant',mismatchcolor='red',covMaxDisplay=50,cropGaps =T,
flipLeftandRight=ifelse(ii > 2,T,F))
p2=formatPlot(retFigs[[1]][[ii]][[2]][[2]],title='Alignment supporting no structural variant',mismatchcolor='red',covMaxDisplay=50,cropGaps =T)
p3=formatPlot(retFigs[[1]][[ii]][[3]][[2]],title='Alignment supporting no structural variant',mismatchcolor='red',covMaxDisplay=50,cropGaps =T)
dev.new()
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
print(p1,vp = viewport(layout.pos.row = 1, layout.pos.col=1))
pushViewport(viewport(layout = grid.layout(2, 2)))
print(p2,vp = viewport(layout.pos.row = 2, layout.pos.col=1))
print(p3,vp = viewport(layout.pos.row = 2, layout.pos.col=2))
}
### This will plot Fig 4A
path <- system.file("extdata", package="targetSeqViewFigureReproduction")
load(file.path(path,"realignments.rda"))
p1=formatPlot(ret[[1]][[11]][[1]][[2]],title='Alignment supporting a structural variant',mismatchcolor='red',covMaxDisplay=50,cropGaps =T)
p2=formatPlot(ret[[1]][[11]][[2]][[2]],title='Alignment supporting no structural variant',mismatchcolor='red',covMaxDisplay=50,cropGaps =T)
p3=formatPlot(ret[[1]][[11]][[3]][[2]],title='Alignment supporting no structural variant',mismatchcolor='red',covMaxDisplay=50,cropGaps =T)
dev.new()
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
print(p1,vp = viewport(layout.pos.row = 1, layout.pos.col=1))
pushViewport(viewport(layout = grid.layout(2, 2)))
print(p2,vp = viewport(layout.pos.row = 2, layout.pos.col=1))
print(p3,vp = viewport(layout.pos.row = 2, layout.pos.col=2))


## Reporduce VennDiagram (Figure 1)
reproduceVenn()

## Reproduce the ROC (Figure 2)

## To do full realignment call reproduceRocLong. For full realignment, speed can be greatly increased if multiple cpus are available.
## To take advantage of multiple cpus, pass the number of available cpus into the function reproduceRocLong using the 'nodes' argument.
## The reproduceRocShort function loads an R object containing the realignments precomputed.

reproduceRocShort()

## Reproduce Repetitive DNA bar chart (supplemental Figure 2)
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


