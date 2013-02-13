\name{reproduceMappability}
\alias{reproduceMappability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
reproduce figures from manuscript
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
reproduceMappability()
}
%- maybe also 'usage' for other objects documented here.
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
function () 
{
    sureselect = read.delim("sureselectformat.bed", header = FALSE)
    rr = read.delim(file = "repeatmaskerbaits.txt", header = F)
    segd = read.delim(file = "genomicSuperDups.txt", header = F)
    repm = GRanges(seqnames = rr[, 5], ranges = IRanges(start = as.numeric(as.character(rr[, 
        6])), end = as.numeric(as.character(gsub("[^0-9]", "", 
        rr[, 7])))), elementMetadata = as.character(rr[, 11]))
    segdGR = GRanges(seqnames = as.character(segd[, 2]), ranges = IRanges(start = as.numeric(as.character(segd[, 
        3])), end = as.numeric(as.character(segd[, 4]))))
    baits2 = cbind(as.character(sureselect[, 1]), sureselect[, 
        2] + 1000, sureselect[, 3] - 1000)
    baits = cbind(as.character(sureselect[, 1]), sureselect[, 
        2] + 100, sureselect[, 3] - 100)
    baitR = IRanges(start = as.numeric(baits[, 2]), end = as.numeric(baits[, 
        3]), names = baits[, 1])
    baitR2 = IRanges(start = as.numeric(baits2[, 2]), end = as.numeric(baits2[, 
        3]), names = baits2[, 1])
    baitGR = GRanges(seqnames = names(baitR), IRanges(baitR))
    baitGR = resize(baitGR, fix = "center", width = 2000)
    baitGR2 = GRanges(seqnames = names(baitR2), IRanges(baitR2))
    allbaits = processbaits(baitR, baitGR)
    baitRcaptured = sort(c(mareGR, mairGR))
    names(baitRcaptured) = seqnames(baitRcaptured)
    baitGRcaptured = baitRcaptured
    capturedevents = processbaits(baitRcaptured, baitGRcaptured)
    newrept = rbind(allbaits, capturedevents)
    newrept2 = cbind(newrept, c(rep("All Baits", nrow(allbaits)), 
        rep("Reported Events", nrow(capturedevents))))
    colnames(newrept2)[4] = "Type"
    return(newrept2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
