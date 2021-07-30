#' Generate hazard ratio with prediction intervals from trial simulation
#'
#' @export
calc_hr <- function(sim){
  sim.nested <-
    sim$sim %>%
    dplyr::group_by(rep) %>%
    nest()
  calc_hr_each_sim <- function(x){
    cfit <- survival::coxph(Surv(pfst, status)~rx, data=x) %>%
      broom::tidy(exponentiate = TRUE) %>%
      dplyr::select(HR=estimate)
    return(cfit)
  }
  sim.hr <-
    sim.nested %>%
    dplyr::mutate(HR = purrr::map(data, calc_hr_each_sim)) %>%
    tidyr::unnest(HR) %>% dplyr::select(-data) %>% dplyr::ungroup()

  ci.range <- 0.95
  quantiles <-
    tibble::tibble(description = c("low", "median", "high"),
                   quantile = c(0.5 - ci.range/2, 0.5, 0.5 + ci.range/2))
  sim.hr.quantile <-
    sim.hr %>%
    dplyr::summarize(low = as.numeric(stats::quantile(HR, probs = 0.5 - ci.range/2, na.rm = TRUE)),
                     median = as.numeric(stats::quantile(HR, probs = 0.5, na.rm = TRUE)),
                     high= as.numeric(stats::quantile(HR, probs = 0.5 + ci.range/2, na.rm = TRUE))) %>%
    pivot_longer(low:high, names_to = "description", values_to = "HR") %>%
    dplyr::inner_join(quantiles,by="description")
  # Output
  out <- list()
  out$sim.hr.quantile <- sim.hr.quantile
  out$sim.hr <- sim.hr
  structure(out, class = c("trialsim.hr"))
}

plot_hr <- function(sim,ci = FALSE){

  # Plot
  ## Generate ggplot object with aes specified using simulated data
  hrplot <- sim$sim.hr
  hrquantile <- sim$sim.hr.quantile

  g <-
    ggplot2::ggplot(hrplot, ggplot2::aes(HR)) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.5, color = "black") +
    ggplot2::geom_vline(data = dplyr::filter(hrquantile, description %in% c("low", "high")),
                        aes(xintercept = HR), lty="dashed", color="red", lwd=1)
  tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                       base_size = 10,
                       padding = unit(c(2, 4), "mm"))
  tbl <- tableGrob(hrquantile, rows=NULL, theme=tt)
  grid.arrange(g, tbl,
               nrow = 2, heights = c(2, 0.5))
  return(g)
}

