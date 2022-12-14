% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_param_confirmed_bor.R
\name{aval_resp}
\alias{aval_resp}
\title{Map Character Response Values to Numeric Values}
\usage{
aval_resp(arg)
}
\arguments{
\item{arg}{Character vector}
}
\value{
\itemize{
\item \code{1} if \code{arg} equals \code{"CR"},
\item \code{2} if \code{arg} equals \code{"PR"},
\item \code{3} if \code{arg} equals \code{"SD"},
\item \code{4} if \code{arg} equals \code{"NON-CR/NON-PD"},
\item \code{5} if \code{arg} equals \code{"PD"},
\item \code{6} if \code{arg} equals \code{"NE"},
\item \code{7} if \code{arg} equals \code{"MISSING"},
\item \code{NA_real_} otherwise
}
}
\description{
Map character response values like \code{"PR"} or \code{"SD"} to numeric values.
}
\examples{

aval_resp(c("CR", "PR", "SD", "NON-CR/NON-PD", "PD", "NE", "MISSING", "ND", NA_character_))
}
\author{
Stefan Bundfuss
}
\concept{utils_fmt}
\keyword{utils_fmt}
