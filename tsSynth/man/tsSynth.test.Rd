% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tsSynth.test.R
\name{tsSynth.test}
\alias{tsSynth.test}
\title{Apply a two-step hypothesis test to determine the appropriate method based on the testing result.}
\usage{
tsSynth.test(
  tsSynth.prep.obj = NULL,
  m = NULL,
  nb = 1000,
  Y1 = NULL,
  X1 = NULL
)
}
\arguments{
\item{m}{The value of subsample size. Noted that the value of m should not be larger than the number of pre-treatment periods.}

\item{nb}{The number of subsampling-bootstrap replication.}

\item{Y1}{A matrix that contains values of the treatment unit over the pre-intervention periods.}

\item{X1}{A matrix that contains values of the control units over the pre-intervention periods.}

\item{synth.prep.obj}{The object that comes from running \code{tsSynth.prep}. This objects contains matrices necessary for running this function. Therefore, once users import this object, there is no need to specify Y1, X1 manually.}
}
\value{
The decision on which method is the most appropriate to apply.
}
\description{
This function provides a formal testing procedure for the parallel trends assumption in the first step, and the application of an appropriate synthetic control method used for estimating accurate ATT of treated unit in the second step.
}
\details{
Users are recommended to supply data object obtained by running the \code{tsSynth.prep}, or import data of both treated and control units over the pretreatment periods manually. Moreover, a subsampling method serves to yield consistent estimators of the test statistic, and thus users need to choose appropriate subsample size and times of subsampling repetition.
}
\examples{
  ##In order to run tsSynth.test() properly, we recommended users
  ##to run tsSynth.prep() for component extraction first.

  #load data from the package.
  data(synth.data)

  #extract relevant components necessary for running tsSynth.test()
  #from wide-format panel data.
  tsSynth.prep.out =
    tsSynth.prep(
     synth.data,
     treatment.identifier = 9,
     controls.identifier = c(2:8),
     time.variable = 1,
     time.prior = c(1984:1989),
     time.post = c(1990:1998)
     )

  ##apply the two-step hypothesis test to determine the most appropriate
  ##method based on the results of the test.
  tsSynth.test.out =
    tsSynth.test(tsSynth.prep.out, m = 4)

}
\references{
Li K, Shankar V (2020) Estimating the causal effect of a digitally native retailer opening a new store: a new two-step synthetic control method. Working Paper, University of Texas at Austin, Austin, TX.

  Li K (2019) Statistical inference for average treatment effects estimated by synthetic control methods. \emph{Journal of American Statistical Association}.
}
