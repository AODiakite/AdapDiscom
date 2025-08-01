% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{lambda_max}
\alias{lambda_max}
\title{Compute Lambda Max for L1 Regularization using KKT Conditions}
\usage{
lambda_max(X, y, Methode = "lasso", robust = 0)
}
\arguments{
\item{X}{Matrix, design matrix}

\item{y}{Vector, response vector}

\item{Methode}{Character, method for computation}

\item{robust}{Integer, 0 for classical, 1 for robust}
}
\value{
Maximum lambda value
}
\description{
Compute Lambda Max for L1 Regularization using KKT Conditions
}
\examples{
# Generate sample data
set.seed(123)
n <- 50; p <- 20
X <- matrix(rnorm(n*p), n, p)
y <- rnorm(n)

# Different methods for lambda_max computation
lmax_lasso <- lambda_max(X, y, Methode = "lasso")
lmax_discom <- lambda_max(X, y, Methode = "discom")

print(paste("Lambda max (lasso):", round(lmax_lasso, 4)))
print(paste("Lambda max (discom):", round(lmax_discom, 4)))
}
