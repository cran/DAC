\name{DAC_Uniform}
\alias{DAC_Uniform}
\title{DAC For Uniform Benchmark Prior And (Skewed) Normal Expert Priors}
\usage{
DAC_Uniform(data,bench.min,bench.max,experts,lb.area,ub.area,step)
}
\description{
Calculates DAC values for multiple experts that provide prections in the form of a (skewed) normal prior. The benchmark is set to be a uniform distribution with user specified boundries.
}

\arguments{
\item{data}{A vector of data needs to be entered here.}
\item{bench.min}{Define the lowerbound for the uniform benchmark prior.}
\item{bench.max}{Define the upperbound for the uniform benchmark prior.}
\item{experts}{A matrix with 1 row per expert. Columns based on input for the dsnorm function of the fGarch package. Column 1; mean, column 2; sd and column 3; skewness parameters}
\item{lb.area}{Define the lowerbound for the area of the parameter space that is to be evaluated.}
\item{ub.area}{Define the upperbound for the area of the parameter space that is to be evaluated.}
\item{step}{The step size or precision by which the area of the parameterspace is evaluated. Smaller step sizes lead to more accurate but slower results.}
}

\references{
Bousquet, N. (2008) Diagnostics of prior-data agreement in applied bayesian analysis.
\emph{Journal of Applied Statistics}, \bold{35(9)}, 1011-1029.

Veen, D., Stoel, D., Schalken, N., van de Schoot, R. (2017). \emph{Using the Data Agreement Criterion to Rank Experts' Beliefs}. Manuscript submitted for publication.
}


\examples{
out <- DAC_Uniform(data=Turnover$actual.data,bench.min = 0,bench.max = 5,
                   experts = Turnover$expert.priors,lb.area = 0,ub.area = 5,step=.01)
out
}
