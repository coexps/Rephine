\name{snvselector}
\alias{snvselector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SNV selector
%%  ~~function to do ... ~~
}
\description{
select snv based on drugs and mutaion information by using lasso and lrtest.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
snvselector(object, p=0.1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A \code{\link[Rephine:Rephine]{Rephine}} objected created by Rephine
%%     ~~Describe \code{object} here~~
}
  \item{p}{
  P value used as cutoff to filter out significant features
  }
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a feature \code{matrix} with logical numbers
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
Frist created by Shiyi Liu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
feature_matrix<-snvselector(Rephine)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
