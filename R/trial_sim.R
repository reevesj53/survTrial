#' Simulation of trial data for Sip-T project
#'
#' @export
#' @param schedule A required vector describing the visit schedule in weeks.
#' @param enrol A required vector that displays the number of subjects projected to be enrolled each    
#'   month.  It incorporates information on the ramp-up period (see example below).  
#' @param rxrate A required vector providing the median PFS for each of two groups 
#'   (reference / treatment).  
#' @param censor.prop Proportion of subjects censored in simulation (e.g. 0.1 for 10%)
#' @param n.rep An integer defining number of simulations
#' @return A `survparamsim` object that contains a data frame (tibble) for the simulated survival profiles with the
#'   following columns:
#'   \itemize{\item \strong{rep}: ID for simulation runs
#'            \item \strong{subjid}: ID for subjects in simulations  
#'            \item \strong{pfst}: PFS time
#'            \item \strong{status}: event status, 0=censored, 1=event
#'   }

trial_sim <- function(schedule, enrol, rxrate, censor.prop=0.1, n.rep = 1000){
  
  if(missing(schedule)) stop("`schedule` needs to be provided")
  if(missing(enrol)) stop("`enrol` needs to be provided")
  if(missing(rxrate)) stop("`enrol` needs to be provided") else if (length(rxrate) !=2) {
    stop("`rxrate` has to be a vector of length two")
  }
    
  lamcens <- rxrate*pcens / (1-pcens)
  nsub <- sum(enrol)
  tot <- nsub*n.rep
  
  tenrol <- vector("list",length(enrol))
  for (i in seq_along(enrol)) tenrol[[i]] <- rep((i-1),times=enrol[i])
  tenrol <- flatten_dbl(tenrol)
  
    ## set randomization variable
    rx <- sample(c(0,1),tot,replace=T)
    
    ## simulate event times for 2 groups
    event <- rexp(tot,rx*rxrate[1]+(1-rx)*rxrate[2])
    
    ## simulate censoring times for 2 groups
    cens <- rexp(tot,rx*lamcens[1]+(1-rx)*lamcens[2])
    
    ## status - event (1) or censoring (0)
    status <- 1*(event<=cens)
    pfst <- pmin(event,cens)
    
    sim <- tibble(pfst,status)
    sim <- sim %>% mutate(rep=rep(1:n.rep,each=nsub),subjid=rep(1:nsub,times=n.rep),
                          death=rbinom(tot,1,0.05),eventt=
                          status*death*pfst + status*(1-death)*
                          forward(pfst,schedule) + (1-status)*backward(pfst,schedule),
                          enrol=tenrol+runif(tot),
                          totalt=enrol+eventt)  
    
    cutoff <- sim %>% group_by(rep) %>% select(rep, subjid, status, totalt) %>% 
      filter(status==1) %>% arrange(rep, totalt) %>% summarize(cut=nth(totalt,226,default=last(totalt)))
   
    sim <- sim %>% inner_join(cutoff,by="rep") %>% arrange(rep,totalt) %>% filter(enrol<=cut) %>% 
      mutate(diff=pmax(totalt-cut,0),status = ifelse(diff>0,0,status), 
             pfst=ifelse(diff>0,pfst-diff,pfst), totalt=ifelse(totalt>cut,cut,totalt)) %>% 
      select(-diff,-cut)
    out <- list()
    out$sim <- sim
    out$cutoff <- cutoff
    structure(out, class = "trialsim")
}

forward <- Vectorize(function(a, vec) vec[which(vec>a)[1]], "a")
backward <- Vectorize(function(a, vec) vec[rev(which(vec<a))[1]], "a")

