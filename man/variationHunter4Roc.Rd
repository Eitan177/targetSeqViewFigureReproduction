\name{variationHunter4Roc}
\alias{variationHunter4Roc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
reproduce figures from manuscript
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
variationHunter4Roc(getTop = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{getTop}{
%%     ~~Describe \code{getTop} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (getTop = FALSE) 
{
    vhrngs <- variationHunterRngs(getTop)
    eventsknown = rngs = list()
    for (inds in 1:5) {
        eventnames = c("NotValidatedR1.txt", "NotValidatedR2.txt", 
            "ValidatedR.txt", "ValidatedR2.txt", "1211canonicalNotValidated.txt")[inds]
        exampleEvents = as.matrix(read.delim(file = eventnames, 
            header = T))
        if (inds == 5) {
            exampleEvents[, 2] = gsub(" ", "", paste("chr", exampleEvents[, 
                2], sep = ""))
            exampleEvents = cbind(exampleEvents[seq(1, nrow(exampleEvents), 
                2), ], exampleEvents[seq(2, nrow(exampleEvents), 
                2), 2:4])
            colnames(exampleEvents)[5:7] = c("Chr2", "Start2", 
                "End2")
        }
        bamFile = vector()
        for (ii in 1:nrow(exampleEvents)) {
            bamFile[ii] = exampleEvents[ii, "Sample"]
            if (length(grep("^1320KB", bamFile[ii])) == 0) {
                numsamp = as.numeric(gsub("[^0-9]", "", bamFile[ii]))
                prefi = ifelse(numsamp < 10, "1320KB000", "1320KB00")
                bamFile[ii] = paste(prefi, numsamp, sep = "")
            }
            bamFile[ii] = paste(as.character(bamFile[ii]), "MultipleAlnsort.bam", 
                sep = "")
        }
        exampleEvents[, "Sample"] = bamFile
        eventsknown[[inds]] = exampleEvents
        rngs[[inds]] = GRanges(seqnames = paste(rep(exampleEvents[, 
            "Sample"], each = 2), gsub("chr", "", as.vector(rbind(exampleEvents[, 
            "Chr1"], exampleEvents[, "Chr2"]))), sep = "."), 
            ranges = IRanges(start = as.numeric(rbind(exampleEvents[, 
                "Start1"], exampleEvents[, "Start2"])), end = as.numeric(rbind(exampleEvents[, 
                "End1"], exampleEvents[, "End2"])))) + 400
        elementMetadata(rngs[[inds]]) = rep(exampleEvents[, "Sample"], 
            each = 2)
    }
    subjectreturns = list()
    retrem = IRanges()
    for (ii in 1:5) {
        ovs = findOverlaps(rngs[[ii]], vhrngs)
        rngnames = seqnames(rngs[[ii]][queryHits(ovs)])
        vhnames = seqnames(vhrngs[subjectHits(ovs)])
        keepovs = which(as.character(rngnames) == as.character(vhnames))
        ovs = ovs[keepovs]
        sovs = split(ovs, (queryHits(ovs) - 1)\%/\%2)
        subjectreturns[[ii]] = lapply(sovs, function(x) {
            ss = split(subjectHits(x), queryHits(x))
            firsts = IRanges(start = ss[[1]] - (ss[[1]] - 1)\%\%2, 
                width = 2)
            if (length(ss) == 1) {
                print("variationhunter junctions not found")
                seconds = firsts
            }
            else {
                seconds = IRanges(start = ss[[2]] - (ss[[2]] - 
                  1)\%\%2, width = 2)
            }
            ff = findOverlaps(firsts, seconds)
            ret = firsts[queryHits(ff)]
            if (nrow(as.matrix(ff)) != 1) {
                print("more than 1 loc will be removed")
            }
            if (length(ret) > 1) {
                matchret = na.omit(match(get("retrem"), ret))
                if (length(matchret) > 0) {
                  print("omitting duplicated hit")
                  ret = ret[-matchret]
                }
                ret = ret[which.max(as.numeric(elementMetadata(vhrngs[start(ret)])[, 
                  2]))]
                assign("retrem", c(retrem, ret), envir = .GlobalEnv)
            }
            ret
        })
        print(length(rngs[[ii]]))
        print(length(unique(queryHits(ovs))))
    }
    ovsN = do.call(c, (unlist(subjectreturns[1:2])))
    elementMetadata(vhrngs[unique(sort(c(start(ovsN), end(ovsN))))])[, 
        1] = "N"
    ovsV = do.call(c, (unlist(subjectreturns[3:4])))
    elementMetadata(vhrngs[unique(sort(c(start(ovsV), end(ovsV))))])[, 
        1] = "V"
    ovsP = do.call(c, (unlist(subjectreturns[5])))
    elementMetadata(vhrngs[unique(sort(c(start(ovsP), end(ovsP))))])[, 
        1] = "P"
    vtprep = c(as.numeric(elementMetadata(vhrngs[elementMetadata(vhrngs)[, 
        1] == "V"])[, 2]), as.numeric(elementMetadata(vhrngs[elementMetadata(vhrngs)[, 
        1] == "P"])[, 2]))
    vfprep = as.numeric(elementMetadata(vhrngs[elementMetadata(vhrngs)[, 
        1] == "N"])[, 2])
    vtprep = vtprep[seq(1, length(vtprep), 2)]
    vfprep = vfprep[seq(1, length(vfprep), 2)]
    vfprep = c(vfprep, rep(0, 17))
    vtprep = c(vtprep, rep(0, 6))
    vtprep = c(vtprep, rep(0, 7))
    vhp = prediction(c(vtprep, vfprep), c(rep(1, length(vtprep)), 
        rep(0, length(vfprep))))
    vhprag = performance(vhp, "tpr", "fpr")
    torem = which(vhprag@y.values[[1]] == 1)
    vhprag@x.values[[1]] = vhprag@x.values[[1]][-torem]
    vhprag@y.values[[1]] = vhprag@y.values[[1]][-torem]
    return(vhprag)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
