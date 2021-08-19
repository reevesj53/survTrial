
enrol <- c(seq(2,10,length.out=5),rep(10,times=3))
schedule <- seq(0,100,4)
rxrate <- c(12,10)
nevent <- 40
sim <- trial_sim(schedule, enrol, rxrate, nevent, adjust=TRUE, trt=c("Sip-T","placebo"),death.prop=0.1,
                 censor.prop=0.1,n.rep=1000)
newdata <- sim$sim %>% dplyr::select(pfst,eventt,status,death) %>% dplyr::filter(death==0) %>%
  dplyr::mutate(adjsched=(eventt<=pfst))

test_that("all eventt>pfst if event, all eventt<pfst if censored (equal if visit maximum)", {
  expect_false(newdata %>% dplyr::filter(status==1) %>% dplyr::pull() %>% any())
  expect_true(newdata %>% dplyr::filter(status==0) %>% dplyr::pull() %>% all())
})

test_that("all events at time zero are censored", {
  expect_false(sim$sim %>% dplyr::transmute(time0=(eventt==0 & status==1)) %>% any())
})

test_that("no PFS times beyond maximum scheduled visit", {
  expect_false(sim$sim %>% dplyr::transmute(timemax=(pmax(eventt,pfst)>max(schedule))) %>% any())
})


