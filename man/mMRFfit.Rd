\name{mMRFfit}
\alias{mMRFfit}

\title{
Estimating mixed Markov fields
}
\description{
Estimates a mixed Markov random field using L1-constrained neighborhood regression.}
\usage{
mMRFfit(data, type, lev, lambda.sel, folds, gam, d, rule.reg = "AND", rule.cat = "OR")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
n x p data matrix
}
  \item{type}{
p-vector containing the types of variable ("g" Gaussian, "e" Exponential, "p" Poisson, "c" Categorical)
}
  \item{lev}{
p-vector containing the number of levels for each variable p (lev=1 for all continuous variables)
}

  \item{lambda.sel}{
Procedure to select the lambda-parameter for L1-penalized regressions. The two options are cross validation "CV" (default) and "EBIC"; While cross validation is used in the paper (see below), the Extended Bayesian Information Criterion (EBIC, see Barber et al., 2015) is a useful alternative in cases where the cross validation-requirement that all categories are present in each fold is not met, because some categories are only present in few cases.
}
  \item{folds}{
The number of folds to be used for cross-validation in case lambda.sel="CV". Defaults to folds=10.}
  \item{gam}{
Gamma hyperparameter when lambda.sel="EBIC". Defaults to gam=.25 (Barber et al., 2015)
}
  \item{d}{
Degrees of augmented interactions. The degree of augmented interactions reflects our belief about the maximal degree in the true graph. (see Loh & Wainwright, 2013)
}
  \item{rule.cat}{
Rule for combining parameters of interactions including categorical variables. The "OR"-rule determines conditional dependence if at least one parameter is non-zero. The "AND"-rule determines conditional dependence if all parameters are non-zero.
}
  \item{rule.reg}{
Rule for combining the two parameters obtained for each edge due to the neighborhood-regression approach. The "OR"-rule determines conditional dependence if at least one of the two parameters is non-zero, the "AND"-rule determines conditional dependence if both parameters are non-zero.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list containing:
\item{adj}{Adjacency matrix}
\item{wadj}{Weighted adjacency matrix; Please note that the parameters of interactions including categorical variables are not interpretable; for the actual interaction parameters inspect wpar.matrix}
\item{wpar.matrix}{Weighted model parameter matrix; This matrix is has the dimension of an overcomplete covariance matrix. By deleting empty columns/rows one obtains the minimal representation. It is given to enable the user to to inspect the exact interaction parameters involving categorical variables, which are (arbitrarily) aggregated in wadj by "rule.cat" (see above)}
\item{lambda}{p x 2 matrix containing the L1 lambda parameter (selected with CV or EBIC) and the tau threshold (see Loh & Wainwright, 2013)}
}
\references{
Barber, R. F., & Drton, M. (2015). High-dimensional Ising model selection with Bayesian information criteria. Electronic Journal of Statistics, 9, 567-607.

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
lev <- c(3, 1, 1, 1) # the categorical variable has three categories

mMRFfit(data=data_mixed,  
        type=type,
        lev=lev,
        lambda.sel="CV",
        d=1, 
        rule.cat="OR",
        rule.reg="AND")
        
        }

