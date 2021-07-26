#' Generate Kaplan-Meier curves with prediction intervals from simulation routine
#'
#' @export

calc_km <- function(sim){

t.last <- survival::survfit(Surv(pfst, status)~1, data = sim$sim)$time %>% max()
t.out <- seq(0, t.last, length.out = 100)

approx_km <- function(x){
  surv <- stats::approx(c(0,x$time), c(1,x$surv), xout=t.out, method="constant", rule=2)$y
  tibble(time = t.out, surv = surv)

  sim.grouped <-
    sim$sim %>%
    dplyr::group_by(rep, rx) %>% nest()


}

}
