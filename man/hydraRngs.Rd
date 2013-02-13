\name{hydraRngs}
\alias{hydraRngs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
reproduce figures from manuscript
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hydraRngs(getTop = FALSE)
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
    events = list()
    cutoff = 2
    path <- system.file("extdata", package = "targetSeqViewReproduceFigures")
    setwd(path)
    ll = list.files(pattern = "1320KB00[0-9][0-9].breaks.final$", 
        full.names = TRUE)
    for (ii in seq_along(ll)) {
        rr = read.delim(file = ll[ii], header = F)
        totest = rr[which(rr[, 19] >= cutoff), ]
        totest = apply(totest, 2, as.numeric)
        if (length(which(is.na(totest[, 1]))) > 0) 
            totest = totest[-which(is.na(totest[, 1])), ]
        print(dim(totest))
        events[[ii]] = cbind(gsub(".breaks.final", "MultipleAlnsort.bam", 
            ll[ii]), totest[, c(1:6, 19)])
        events[[ii]][, c(3, 6)] = as.numeric(events[[ii]][, c(3, 
            6)]) - 500
        events[[ii]][, c(4, 7)] = as.numeric(events[[ii]][, c(4, 
            7)]) + 500
        events[[ii]][, c(2, 5)] = paste("chr", events[[ii]][, 
            c(2, 5)], sep = "")
        colnames(events[[ii]])[1:8] = c("Sample", "Chr1", "Start1", 
            "End1", "Chr2", "Start2", "End2", "score")
    }
    hydraevents = vector()
    for (ii in seq_along(ll)) {
        hydraevents = rbind(hydraevents, events[[ii]])
    }
    hydrarngs = GRanges(seqnames = paste(rep(hydraevents[, "Sample"], 
        each = 2), gsub("chr", "", as.vector(rbind(hydraevents[, 
        "Chr1"], hydraevents[, "Chr2"]))), sep = "."), ranges = IRanges(start = as.numeric(rbind(hydraevents[, 
        "Start1"], hydraevents[, "Start2"])), end = as.numeric(rbind(hydraevents[, 
        "End1"], hydraevents[, "End2"]))))
    elementMetadata(hydrarngs)[, 1] = "U"
    elementMetadata(hydrarngs)[, 2] = rep(hydraevents[, "score"], 
        each = 2)
    if (getTop) {
        hydrarngs = hydrarngs[order(as.numeric(elementMetadata(hydrarngs)[, 
            2]), decreasing = T)][1:tophits]
    }
    return(hydrarngs)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
