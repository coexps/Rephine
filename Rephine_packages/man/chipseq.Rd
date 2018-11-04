\name{chipseq}
\alias{chipseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ChiPSeq caculation
%%  ~~function to do ... ~~
}
\description{
Use elastic net and permuation to caculate TF erichment.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
chipseq(chip = matrix(0, 0, 0), stage_all_partial = matrix(0, 0, 0), genes = NULL, samples = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chip}{
  a ChiPseq \code{matrix}
%%     ~~Describe \code{chip} here~~
}
  \item{stage_all_partial}{
  partial correlation genertated by \code{\link[Rephine:parcol]{parcol}}
%%     ~~Describe \code{stage_all_partial} here~~
}
  \item{genes}{
  \code{vector} of genes
%%     ~~Describe \code{genes} here~~
}
  \item{samples}{
  \code{vector} of samples
%%     ~~Describe \code{samples} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a list contains pvalue from elastic net and permutaion.
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
First created by Shiyi Liu
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
ptable<-chipseq(chip=chip,stage_all_partial=stage_all_partial)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
