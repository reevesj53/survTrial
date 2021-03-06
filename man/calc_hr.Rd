% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_hr.R
\name{calc_hr}
\alias{calc_hr}
\title{Generate statistics from simulation runs on hazard ratios}
\usage{
calc_hr(sim, ci.range = 0.95)
}
\arguments{
\item{sim}{A \code{trialsim} class object generated by \code{\link[=trial_sim]{trial_sim()}} function.}

\item{ci.range}{Prediction interval for simulated hazard ratios (defaults to 0.95).}
}
\value{
A \code{trialsim.hr} object that contains a data frame (tibble) for the hazard ratios of each simulation
with the following columns:
\itemize{\item \strong{rep}: ID for simulation runs
\item \strong{HR}: hazard ratio for each simulation
}
And a data frame containing the summary statistics for hazard ratios across simulations:
\itemize{\item \strong{description}: Labels for hazard ratios: lower limit \ median \ upper limit
\item \strong{HR}: Statistic for hazard ratios
\item \strong{quantile}: Quantile of simulations (corresponding to \code{ci.range} parameter)
}
}
\description{
Generate statistics from simulation runs on hazard ratios
}
\details{
Extracts median and prediction interval by treatment group across all repeated simulations for hazard ratio.
}
\examples{
# Generate some trial simulation data
enrol <- c(seq(2,10,length.out=5),rep(10,times=3))
schedule <- seq(0,100,4)
sim <- trial_sim(schedule, enrol, c(12,10), 40, adjust=TRUE, trt=c("Treatment","Placebo"),
death.prop=0.1, censor.prop=0.1, n.rep=1000)
# Extract hazard ratios from simulated data
\donttest{sim.hr <- calc_hr(sim)}
}
