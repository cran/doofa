% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doofa.R
\name{gen.design2}
\alias{gen.design2}
\title{Repeat the process of design generation using doofa.pwo several times and return the best design}
\usage{
gen.design2(n, m, num.repeat = 10)
}
\arguments{
\item{n}{number of runs, a positive integer}

\item{m}{number of components, a positive integer}
\item{num.repeat}{number of repeats, a positive integer}
}
\value{
a design with D-efficiency
}
\description{
Repeat the process of design generation using doofa.pwo several times and return the best design
}
\examples{
gen.design2(n = 5, m = 3, num.repeat = 10);
}
