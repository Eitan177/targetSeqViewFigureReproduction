\name{gasvRngs}
\alias{gasvRngs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
reproduce figures from manuscript

}
\description{
  reproduce figures from manuscript
}
\usage{
gasvRngs(getTop = FALSE)
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
    gevents = list()
    cutoff = 2
    path <- system.file("extdata", package = "targetSeqView")
    setwd(path)
    clustf = system("find . -name \"BAMToGASV_AMBIG.gasv.combined.in.clusters\"", 
        intern = T)
    clustf = grep("1320KB", clustf, value = T)
    for (ii in seq_along(clustf)) {
        clust = read.delim(file = clustf[ii], header = T)
        totestclust = clust[which(clust[, 2] >= cutoff), ]
        gevents[[ii]] = cbind(paste(grep("1320KB", strsplit(clustf[ii], 
            "/")[[1]], value = T), "MultipleAlnsort.bam", sep = ""), 
            totestclust[, 6], as.numeric(gsub(",.+", "", totestclust[, 
                8])) - 550, as.numeric(gsub(",.+", "", totestclust[, 
                8])) + 550, totestclust[, 7], as.numeric(gsub(",.+", 
                "", gsub("^[0-9]+, ", "", totestclust[, 8]))) - 
                550, as.numeric(gsub(",.+", "", gsub("^[0-9]+, ", 
                "", totestclust[, 8]))) + 550, totestclust[, 
                2])
        gevents[[ii]][, c(2, 5)] = paste("chr", gevents[[ii]][, 
            c(2, 5)], sep = "")
        colnames(gevents[[ii]])[1:8] = c("Sample", "Chr1", "Start1", 
            "End1", "Chr2", "Start2", "End2", "score")
    }
    gasvevents = vector()
    for (ii in seq_along(clustf)) {
        gasvevents = rbind(gasvevents, gevents[[ii]])
    }
    gasvrngs = GRanges(seqnames = paste(paste("/amber1/archive/sgseq/data/KathyBurns/hydraBreakFiles/", 
        rep(gasvevents[, "Sample"], each = 2), sep = ""), gsub("chr", 
        "", as.vector(rbind(gasvevents[, "Chr1"], gasvevents[, 
            "Chr2"]))), sep = "."), ranges = IRanges(start = as.numeric(rbind(gasvevents[, 
        "Start1"], gasvevents[, "Start2"])), end = as.numeric(rbind(gasvevents[, 
        "End1"], gasvevents[, "End2"]))))
    elementMetadata(gasvrngs)[, 1] = "U"
    elementMetadata(gasvrngs)[, 2] = rep(gasvevents[, "score"], 
        each = 2)
    g1s = gasvrngs[seq(1, length(gasvrngs), 2)]
    g2s = gasvrngs[seq(2, length(gasvrngs), 2)]
    tocloseSoremove = (which(width(pgap(ranges(g1s), ranges(g2s))) < 
        1000))
    closeremove = tocloseSoremove * 2
    closerem = sort(c(closeremove - 1, closeremove))
    gasvrngs = gasvrngs[-closerem]
    if (getTop) {
        gasvrngs = gasvrngs[order(as.numeric(elementMetadata(gasvrngs)[, 
            2]), decreasing = T)][1:tophits]
    }
    return(gasvrngs)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
