\name{targetSeqView4Roc}
\alias{targetSeqView4Roc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
reproduce figures from manuscript
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
targetSeqView4Roc(ret)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ret}{
%%     ~~Describe \code{ret} here~~
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
function (ret) 
{
    likelihoodScores = ret[[2]]
    nonvalids = likelihoodScores[1:38]
    valids = likelihoodScores[39:length(likelihoodScores)]
    labels = c(rep(0, length(nonvalids)), rep(1, length(valids)))
    prag = prediction(likelihoodScores, labels)
    pred = performance(prag, "tpr", "fpr")
    return(pred)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
