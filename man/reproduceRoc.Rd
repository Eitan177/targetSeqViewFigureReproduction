\name{reproduceRoc}
\alias{reproduceRoc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
reproduce figures from manuscript
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
reproduceRoc(ret)
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
    pred <- targetSeqView4Roc(ret)
    gprag <- gasv4Roc()
    vhprag <- variationHunter4Roc()
    hydraprag <- hydra4Roc()
    xnum = xnumMine = 38
    ynum = 25
    ynumMine = 26
    x = c(pred@x.values[[1]] * xnumMine, hprag@x.values[[1]] * 
        xnum, hydrahprag@x.values[[1]] * xnum, vhprag@x.values[[1]] * 
        xnum)
    y = c(pred@y.values[[1]] * ynumMine, hprag@y.values[[1]] * 
        ynum, hydrahprag@y.values[[1]] * ynum, vhprag@y.values[[1]] * 
        ynum)
    meth = c(rep("SeqView", length(pred@x.values[[1]])), rep("GASV", 
        length(hprag@x.values[[1]])), rep("HYDRA", length(hydrahprag@x.values[[1]])), 
        rep("VariationHunter", length(vhprag@x.values[[1]])))
    forRoc = data.frame(x, y, meth)
    colnames(forRoc)[3] = "group"
    groc = ggplot(forRoc, aes(x = x, y = y, colour = group)) + 
        geom_line(cex = 2) + xlab("Failed Validation") + ylab("Validated and Canonical V(D)J") + 
        scale_x_continuous(breaks = seq(5, 35, 5), labels = seq(5, 
            35, 5)) + scale_y_continuous(breaks = seq(2, 25, 
        2), labels = seq(2, 25, 2))
    print(groc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
