\name{R phine}
\alias{Rephine}
\title{
Raphine Constructor
}
\description{
Create a \code{Raphine} object from drug activity, expression, snv, cnv and tissue matrix.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Raphine(drugs = matrix(0, 0, 0), exp = matrix(0, 0, 0), snv = matrix(0, 0, 0), cnv = matrix(0, 0, 0), tissue = matrix(0, 0, 0), cellline = NULL, genes = NULL, drug = NULL, remove.zeros = FALSE)
}
  \arguments{
  \item{drugs}{
  \code{matrix} of drug activity, cols are drugs, rows are ...
}
  \item{exp}{
  \code{matrix} of gene expression though tissues/cell lines
}
  \item{snv}{
  \code{matrix} of snv information
%%     ~~Describe \code{snv} here~~
}
  \item{cnv}{
  cnv \code{matrix}
}
  \item{tissue}{
  tissue \code{matrix}
%%     ~~Describe \code{tissue} here~~
}
  \item{cellline}{
  celline \code{vector}
%%     ~~Describe \code{cellline} here~~
}
  \item{genes}{
  gene \code{vector}
%%     ~~Describe \code{genes} here~~
}
  \item{drug}{
  drug \code{vector}
%%     ~~Describe \code{drug} here~~
}
  \item{remove.zeros}{
  ogical, whether to remove rows with zero counts
%%     ~~Describe \code{remove.zeros} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a \code{Rephine} object
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
First created by Shiyi Liu.
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
data(Rephine)
Rep<-Rephine(drugs=drugs,exp=exp,snv=snv,cnv=cnv,tissue=tissue)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
