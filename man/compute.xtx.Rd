% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{compute.xtx}
\alias{compute.xtx}
\title{Compute X'X Matrix}
\usage{
compute.xtx(x, robust = 0, k_value = 1.5)
}
\arguments{
\item{x}{Matrix, input data matrix}

\item{robust}{Integer, 0 for classical estimate, 1 for Huber robust estimate}

\item{k_value}{Numeric, tuning parameter for Huber function}
}
\value{
Covariance matrix
}
\description{
Compute X'X Matrix
}
\examples{
# Create sample data with missing values
set.seed(123)
x <- matrix(rnorm(100), 20, 5)
x[sample(100, 10)] <- NA  # Introduce missing values

# Classical covariance estimation
xtx_classical <- compute.xtx(x, robust = 0)
print(round(xtx_classical, 3))

# Robust covariance estimation
xtx_robust <- compute.xtx(x, robust = 1, k_value = 1.5)
print(round(xtx_robust, 3))
}
