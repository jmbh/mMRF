\name{mMRFfit}
\alias{mMRFfit}

\title{
Estimating mixed Markov fields
}
\description{
Estimates a mixed Markov random field using a L1-constraint.}
\usage{
mMRFfit(data, type, levels, d, rule.w = "OR", rule.b = "OR")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
n x p data matrix
}
  \item{type}{
p-vector containing the types of variable ("g" Gaussian, "e" Exponential, "p" Poisson, "c" Categorical)
}
  \item{levels}{
p-vector containing the number of levels for each variable p (level=1 for all continuous variables)
}
  \item{d}{
Degrees of augmented interactions (See Loh \& Wainwright, 2013)
}
  \item{rule.w}{
Combining parameters within regression.
}
  \item{rule.b}{
Combining parameters between regression}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list containing:
\item{wadj}{Weighted adjacency matrix}
\item{adj}{Adjacency matrix}
\item{wpar.matrix}{Weighted model parameter matrix}
\item{par.matrix}{Binarized model parameter matrix}
\item{lambda}{p x 2 matrix containing the L1 lambda parameter (selected with 10fold CV) and the tau threshold}
}
\references{
Loh, P. L., & Wainwright, M. J. (2013). Structure estimation for discrete graphical models: Generalized covariance matrices and their inverses. The Annals of Statistics, 41(6), 3022-3049.

Yang, E., Baker, Y., Ravikumar, P., Allen, G., & Liu, Z. (2014). Mixed graphical models via exponential families. In Proceedings of the Seventeenth International Conference on Artificial Intelligence and Statistics (pp. 1042-1050).
}
\author{
Jonas Haslbeck}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

head(data_mixed) #example data

type <- c("c", "g", "p", "g") #data contains a categorical, a gaussian, a poisson and a gaussian variable in that order
levels <- c(3, 1, 1, 1) # the categorical variable has three categories

mMRFfit(data=data_mixed,  
        type=type,
        levels=levels, 
        d=1, 
        rule.w="OR",
        rule.b="OR")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
