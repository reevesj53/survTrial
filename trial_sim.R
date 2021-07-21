#' Simulation of trial data for Sip-T project
#'

trial_sim <- function(schedule, enrol, censor.prop=0.1, n.rep = 1000){
  sched <- seq(0,600,12)
  if(missing(schedule)) stop("`schedule` needs to be provided")
  if(missing(enrol)) stop("`enrol` needs to be provided")
  lam <- log(2)/c(7.7,5)
  pcens <- 0.1
  lamcens <- lam*pcens / (1-pcens)
  n.rep <- 6
  nsub <- sum(enrol)
  tot <- nsub*n.rep
  
  enrol <- rep(10,each=28)
  tenrol <- vector("list",length(enrol))
  for (i in seq_along(enrol)) tenrol[[i]] <- rep((i-1),times=enrol[i])
  tenrol <- flatten_dbl(tenrol)
  
  
    ## set randomization variable
    rx <- sample(c(0,1),tot,replace=T)
    
    ## simulate event times for 2 groups
    event <- rexp(tot,rx*lam[1]+(1-rx)*lam[2])
    
    ## simulate censoring times for 2 groups
    cens <- rexp(tot,rx*lamcens[1]+(1-rx)*lamcens[2])
    
    ## status - event (1) or censoring (0)
    status <- 1*(event<=cens)
    pfst <- pmin(event,cens)
    
    sim <- tibble(pfst,status)
    sim <- sim %>% mutate(rep=rep(1:n.rep,each=nsub),subjid=rep(1:nsub,times=n.rep),
                          death=rbinom(tot,1,0.05),eventt=
                          status*death*pfst + status*(1-death)*
                          forward(pfst,sched) + (1-status)*backward(pfst,sched),
                          enrol=tenrol+runif(tot),
                          totalt=enrol+eventt)  
    
    cutoff <- sim %>% group_by(rep) %>% select(rep, subjid, status, totalt) %>% 
      filter(status==1) %>% arrange(rep, totalt) %>% summarize(cut=nth(totalt,226,default=last(totalt)))
   
    sim <- sim %>% inner_join(cutoff,by="rep") %>% arrange(rep,totalt) %>% filter(enrol<=cut) %>% 
      mutate(diff=pmax(totalt-cut,0),status = ifelse(diff>0,0,status), 
             pfst=ifelse(diff>0,pfst-diff,pfst), totalt=ifelse(totalt>cut,cut,totalt)) %>% 
      select(-diff)
     
}

forward <- Vectorize(function(a, vec) vec[which(vec>a)[1]], "a")
backward <- Vectorize(function(a, vec) vec[rev(which(vec<a))[1]], "a")
