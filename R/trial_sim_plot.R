#' @export
plot_sim <- function(sim) UseMethod("plot_sim")
#' @export
plot_sim.trialsim.km <- function(sim, ci = FALSE){

  # Plot
  ## Generate ggplot object with aes specified using simulated data
  #sim <- sim.km
  kmplot <- sim$sim.km
  kmquantile <- sim$median.quantile
  km1 <- kmquantile %>% filter(rx==0) %>% select(-rx)
  km2 <- kmquantile %>% filter(rx==1) %>% select(-rx)

  g <-
    ggplot2::ggplot(kmplot, ggplot2::aes(median)) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.5, color = "black") +
    ggplot2::geom_vline(data = dplyr::filter(kmquantile, description %in% c("int_low", "int_high")),
                        aes(xintercept = median), lty="dashed", color="red", lwd=1) +
    ggplot2::facet_wrap(~rx)
  tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                       base_size = 10,
                       padding = unit(c(2, 4), "mm"))
  tbl1 <- tableGrob(km1, rows=NULL, theme=tt)
  tbl2 <- tableGrob(km2, rows=NULL, theme=tt)
  hlay <- rbind(c(1,1),
                c(2,3))
  g <- grid.arrange(g, tbl1, tbl2, layout_matrix=hlay,
                    heights = c(2, 1))
  return(g)
}

