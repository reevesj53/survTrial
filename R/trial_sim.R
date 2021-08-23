#' Simulation of trial data
#'
#' The main function to generate two-arm clinical trial PFS data assuming an exponential distribution.
#' @export
#' @param schedule A required vector describing the visit schedule in weeks.  Simulated events occurring
#' after the last scheduled visit are censored.
#' @param enrol A required vector that displays the number of subjects projected to be enrolled each
#'   month.  It incorporates information on the ramp-up period (see example below).
#' @param rxrate A required vector providing the median PFS for each of two groups in months
#'   (treatment / placebo).
#' @param nevent A required integer providing the event number that defines the data cut-off.
#' @param adjust A logical to specify whether progression events are moved forward to next
#' scheduled visit, and censored events remain at last available visit (defaults to TRUE).
#' @param trt An optional character vector specifying names of 2 treatment arms (order corresponds to `rxrate`).
#' @param death.prop Specifies proportion of simulated PFS events that are deaths (i.e. simulated
#' times of death are not adjusted to `schedule` and `adjust` parameter has no effect).
#' @param censor.prop Proportion of subjects censored in simulation.
#' @param n.rep An integer defining number of simulations.
#' @return A `trialsim` object that contains a data frame (tibble) for the simulated survival profiles with the
#'   following columns:
#'   \itemize{\item \strong{rep}: ID for simulation runs
#'            \item \strong{rx}: Treatment group
#'            \item \strong{subjid}: ID for subjects in simulations
#'            \item \strong{pfst}: Simulated PFS time (weeks)
#'            \item \strong{status}: Event status, 0=censored, 1=event
#'            \item \strong{death}: Death status, 0=alive, 1=dead
#'            \item \strong{eventt}: Simulated adjusted PFS Time (weeks) (output if `adjust=TRUE`)
#'            \item \strong{enrol}: Enrollment time (weeks)
#'            \item \strong{totalt}: Total time: enrollment time plus PFS time (weeks)
#'            }
#'   And a data frame containing the subject with the event that represents the data cut-off based on
#'   the `nevent` parameter (for all simulations).
#'
#' @details
#' [trial_sim()] returns simulated PFS data from a hypothetical clinical trial using the ramp-up enrollment provided
#' in `schedule`. Censoring proportion and proportion of PFS events that result in death can be specified.  The option to
#' move progression events forward to the next scheduled visit (and censored events back to previous visit) is also
#' provided.
#'
#' @examples
#' # Generate enrollment schedule, with ramp-up period of 5 months, 10 subjects per month thereafter,
#' # up to a total of 60 subjects with 1:1 randomization
#' enrol <- c(seq(2,10,length.out=5),rep(10,times=3))
#' # Input visit schedule for tumor assessments in weeks
#' schedule <- seq(0,100,4)
#' # Set Median PFS to 10 months for Placebo, and 12 months for Sip-T
#' rxrate <- c(12,10)
#' # Data cut-off is set at event number 40, proportion of events that are deaths is set at 10% and
#' # censoring rate is 10%.  Accept default method of moving progression events forward to next
#' # scheduled visit.
#' # Number of trial simulations = 1000.
#' nevent <- 40
#' sim <- trial_sim(schedule, enrol, rxrate, nevent, adjust=TRUE, trt=c("Treatment","Placebo"),
#' death.prop=0.1, censor.prop=0.1, n.rep=1000)

trial_sim <- function(schedule, enrol, rxrate, nevent, adjust=TRUE, trt=c("treatment","placebo"), death.prop=0.1,
                      censor.prop=0.1, n.rep = 1000){

  if(missing(schedule)) stop("`schedule` needs to be provided")
  if(missing(enrol)) stop("`enrol` needs to be provided")
  if(missing(rxrate)) stop("`rxrate` needs to be provided") else if (length(rxrate) !=2) {
    stop("`rxrate` has to be a vector of length two")
  }
  if(missing(nevent)) stop("`nevent` needs to be provided")

  week <- (365.25/12)/7
  rxrate <- log(2)/(rxrate*week)
  lamcens <- rxrate*censor.prop / (1-censor.prop)
  nsub <- sum(enrol)
  tot <- nsub*n.rep

  tenrol <- vector("list",length(enrol))
  for (i in seq_along(enrol)) tenrol[[i]] <- rep((i-1),times=enrol[i])
  tenrol <- purrr::flatten_dbl(tenrol)

    ## set randomization variable
    rx <- sample(c(0,1),tot,replace=T)

    ## simulate event times for 2 groups
    event <- stats::rexp(tot,rx*rxrate[1]+(1-rx)*rxrate[2])

    ## simulate censoring times for 2 groups
    cens <- stats::rexp(tot,rx*lamcens[1]+(1-rx)*lamcens[2])

    ## status - event (1) or censoring (0)
    status <- 1*(event<=cens)
    pfst <- pmin(event,cens)
    status <- dplyr::if_else(pfst>max(schedule),0,status)
    pfst <- pmin(pfst,max(schedule))

    sim <- tibble::tibble(pfst,status,rx=trt[2-rx])
    sim <- sim %>% dplyr::mutate(rep=rep(1:n.rep,each=nsub),subjid=rep(1:nsub,times=n.rep),
                          death=stats::rbinom(tot,1,death.prop) * status,
                          enrol=(tenrol+stats::runif(tot))*week,
                          totalt=enrol+pfst)

    cutoff <- sim %>% dplyr::group_by(rep) %>% dplyr::filter(status==1) %>% dplyr::arrange(rep,totalt) %>%
    dplyr::slice(nevent)
    cut1 <- cutoff %>% dplyr::select(rep, totalt) %>% dplyr::rename(totalmax=totalt)

    sim <- sim %>% dplyr::left_join(cut1) %>% dplyr::filter(is.na(totalmax) | enrol<=totalmax) %>%
      dplyr::mutate(diff=dplyr::if_else(is.na(totalmax),0,pmax(totalt-totalmax,0)),status = dplyr::if_else(diff>0,0,status),
                    dplyr::across(c(pfst,totalt),~.x-diff)) %>% dplyr::relocate(rep,rx,subjid) %>%
                    dplyr::select(!c(totalmax,diff))
    if (adjust) {
      sim <- sim %>% dplyr::mutate(eventt=status*death*pfst + status*(1-death)*forward(pfst,schedule)
                            + (1-status)*backward(pfst,schedule),totalt=enrol+eventt)
    }
    cutoff <- sim %>% dplyr::semi_join(cutoff,by=c("rep","subjid"))

    out <- list()
    out$sim <- sim
    out$cutoff <- cutoff
    structure(out, class = "trialsim")
}

forward <- Vectorize(function(a, vec) vec[which(vec>=a)[1]], "a")
backward <- Vectorize(function(a, vec) vec[rev(which(vec<=a))[1]], "a")
