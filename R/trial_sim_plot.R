#' Generic graph function
#'
#' Plots simulation results for both Kaplan-Meier medians by treatment group, and hazard ratios generated
#' from simulated trials.
#' @export
#' @param sim An object from [calc_km()] or [calc_hr()]
#' @param ci A logical specifying whether confidence bounds are shown on the histogram plot (default=FALSE)
#' @param summary A logical specifying whether a summary table is presented on the graph (default=FALSE)
plot_sim <- function(sim, ci, summary) UseMethod("plot_sim")

#' @export
#' @describeIn plot_sim Simulation results presented for Kaplan-Meier medians
plot_sim.trialsim.km <- function(sim, ci = FALSE, summary = FALSE){
  kmplot <- sim$sim.km
  kmquantile <- sim$median.quantile
  upperlim <- stats::quantile(kmplot$median, 0.995, na.rm=TRUE)
  g <- kmplot %>%
    ggplot2::ggplot(ggplot2::aes(x=median, fill=rx)) +
    ggplot2::geom_histogram(color="#e9ecef", alpha=0.6, na.rm=TRUE, bins=30) +
    ggplot2::coord_cartesian(xlim=c(0,upperlim))+
    ggplot2::facet_wrap(~rx)+
    ggplot2::labs(x="Simulated KM medians (weeks)",y="N simulated trials")+
    ggplot2::theme(legend.position = "none")
  if(ci) {
  g <- g +
    ggplot2::geom_vline(data = dplyr::filter(kmquantile, description %in% c("sim_low", "sim_high")),
                        ggplot2::aes(xintercept = KM_median, color=rx), lty="dashed", lwd=1)
  }
  if(!summary) {
    g <- g + ggplot2::theme(plot.margin=grid::unit(c(0.5,0.5,2,0.5),"cm"))
    return(g)
  } else
    {
    g <- g + ggplot2::theme(plot.margin=grid::unit(c(0.5,0.5,-0.5,0.5),"cm"))
    tt <- gridExtra::ttheme_default(base_size = 9)
    kmquantile <- kmquantile %>% dplyr::mutate(KM_median=sprintf("%.1f",KM_median)) %>%
    dplyr::select(-description) %>% dplyr::rename(`KM median`=KM_median,treatment=rx) %>%
    dplyr::group_split(treatment)
    tbl1 <- gridExtra::tableGrob(kmquantile[[1]], rows=NULL, theme=tt)
    tbl1 <- gtable::gtable_add_padding(tbl1, grid::unit(c(0,0,0,2), "cm"))
    tbl2 <- gridExtra::tableGrob(kmquantile[[2]], rows=NULL, theme=tt)
    hlay <- rbind(c(1,1),c(2,3))
    return(gridExtra::grid.arrange(g, tbl1, tbl2, layout_matrix=hlay,
                      heights = c(2, 1)))
    }
}

#' @export
#' @describeIn plot_sim Simulation results presented for hazard ratios
plot_sim.trialsim.hr <- function(sim, ci = FALSE, summary = FALSE){
  hrplot <- sim$sim.hr
  hrquantile <- sim$sim.hr.quantile
  graphlim <- stats::quantile(hrplot$HR, c(0.005,0.995))
  g <- hrplot %>%
    ggplot2::ggplot(ggplot2::aes(x=HR)) +
    ggplot2::geom_histogram(color="#e9ecef", alpha=0.6, na.rm=TRUE, bins=30) +
    ggplot2::coord_cartesian(xlim=graphlim)+
    ggplot2::labs(x="Simulated HR",y="N simulated trials")+
    ggplot2::theme(legend.position = "none")+
    ggplot2::scale_x_continuous(trans="log2")+
    ggplot2::geom_vline(xintercept=1,lwd=1)
  if(ci) {
    g <- g +
      ggplot2::geom_vline(data = dplyr::filter(hrquantile, description %in% c("sim_low", "sim_high")),
                          ggplot2::aes(xintercept = HR), lty="dashed", color = "red", lwd=1)
  }
  if(!summary) {
    g <- g + ggplot2::theme(plot.margin=grid::unit(c(0.5,0.5,2,0.5),"cm"))
    return(g)
  } else
  {
    g <- g + ggplot2::theme(plot.margin=grid::unit(c(0.5,0.5,-0.5,0.5),"cm"))
    tt <- gridExtra::ttheme_default(base_size = 9)
    hrquantile <- hrquantile %>% dplyr::mutate(HR=sprintf("%.3f",HR)) %>%
      dplyr::select(-description)
    tbl1 <- gridExtra::tableGrob(hrquantile, rows=NULL, theme=tt)
    return(gridExtra::grid.arrange(g, tbl1, nrow=2,
                        heights = c(2, 1)))
  }
}

