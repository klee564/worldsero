VEdata2 <- preprocVE2()
bmodel <- fitstan(VEdata2)
mydat <- df_to_standat2(VEdata2)


function(mydat){
  #simuldata로 수정하기
  mu0 <- runif(1,-4,0)
  dmu0 <- runif(1,0,2)
  tau <- runif(1,0,1)
  kappa <- runif(1,0,1)
  pkappa <- runif(1,0,1)

  head(VEdata2)

  K <- mydat$K
  J1 <- mydat$J1
  J2 <- mydat$J2
  g1 <- mydat$g1
  g2 <- mydat$g2
  v1 <- mydat$v1
  v2 <- mydat$v2

  mu <- rep(0,K)
  pmu <- rep(0,K)
  y1 <- rep(0,J1)
  y2 <- rep(0,J2)
  for(k in 1:K){
    mu[k] <- rnorm(1,mu0,kappa);
    pmu[k] <- rnorm(1,mu0+dmu0,pkappa);
  }
  for(j in 1:J1){
    y1[j] <- rnorm(1,pmu[g1[j]],sqrt(tau^2+v1[j]));
    v1[j] <- runif(1,0,1)
  }
  for(j in 1:J2){
    y2[j] <- rnorm(1,mu[g2[j]],sqrt(tau^2+v2[j]));
    v2[j] <- runif(1,0,1)
  }


  updated_mydat <- mydat
  updated_mydat$y1 <- y1
  updated_mydat$y2 <- y2
  fit <- stan(file = "hierarchicalmeta2.stan",
              data = updated_mydat, iter = 8000, warmup = 5000, chains = 4)


  rnames <- rownames(summary(fit)$summary)

  mucred <- summary(fit)$summary[stringr::str_detect(rnames,"^mu\\["),c(4,8)]
  1:length(mu) %>% purrr::map(~between(mu[.x],mucred[.x,1],mucred[.x,2])) %>% unlist()

  pmucred <- summary(fit)$summary[stringr::str_detect(rnames,"^pmu\\["),c(4,8)]
  1:length(pmu) %>% purrr::map(~between(pmu[.x],pmucred[.x,1],pmucred[.x,2])) %>% unlist()


  tbetween(mu0,summary(fit)$summary[rnames == "mu0",c(4,8)])
  tbetween(dmu0,summary(fit)$summary[rnames == "dmu0",c(4,8)])

  tbetween(kappa,summary(fit)$summary[rnames == "kappa",c(4,8)])
  tbetween(pkappa,summary(fit)$summary[rnames == "pkappa",c(4,8)])
  tbetween(tau,summary(fit)$summary[rnames == "tau",c(4,8)])


  list(truparam =
         list(mu0=mu0,dmu0=dmu0,kappa=kappa,pkappa=pkappa,tau=tau,mu=mu,pmu=pmu),
       fit=fit)

}


traceplot(fit, pars = c("mu0", "dmu0","kappa","pkappa","tau"), inc_warmup = TRUE)

summary(fit)$summary[1:7,c(4,8)]
summary(fit)$summary[1:7,4] <= mu


tbetween <- function(value,vec){
  between(value,vec[1],vec[2])
}
