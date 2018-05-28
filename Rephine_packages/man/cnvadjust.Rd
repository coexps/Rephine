\name{cnvadjust}
\alias{cnvadjust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CNV adjust expression
%%  ~~function to do ... ~~
}
\description{
Use linear regression to adjust expression with CNV data.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cnvadjust(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A \code{\link[Rephine:Rephine]{Rephine}} object
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a \code{matrix} contains adjusted expression
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
reg_factor<-cnvadjust(Rephine)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
