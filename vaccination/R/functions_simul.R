#'
#' @examples
#' train_list <- train_input(1,del=21)
#' inputlist <- input_multinom()
#' mydat <- preprocVE2()  %>% df_to_standat2
#'
#'
runsimul <- function(train_list,inputlist,mydat){
  library(rstan)
  library(posterior)
  library(dplyr)

  pois_alpha <- runif(1,-1,0)
  pois_beta <- runif(1,0,1) # 이전엔 0.5
  multi_beta <- runif(1,0.5,1.5)

  simulres_pois <- pois_simul(train_list,pois_alpha,pois_beta)
  simulres_multi <- multi_simul(inputlist,beta=multi_beta)

  VEres <- VE_simul(mydat)
  rbind(eval_simul1(simulres_multi,simulres_pois),
        eval_VEsimul(VEres$fit,VEres$truparam))
}



#'
#' @examples
#' param <- gen_paramXY()
#' train_list <- train_input(1,del=21)
#' mcmc.out <- pois_simul(train_list,param$pois_coef$pois_alpha,param$pois_coef$pois_beta)
#'
#'
pois_simul <- function(train_list,alpha=-1,beta=0.1){
  library(nimble)

  traindata <- train_list$datalist
  trainresp <- 1:train_list$constants$n %>%
    purrr::map(~rpois(1,with(traindata,Xadj[.x]+pop[.x]*(
      weightedV[.x,1]*exp(alpha+beta*log(Z[.x,1]/pop[.x]*100))+
        weightedV[.x,2]*exp(alpha+beta*log(Z[.x,2]/pop[.x]*100))+
        weightedV[.x,3]*exp(alpha+beta*log(Z[.x,3]/pop[.x]*100)))/100))) %>%
    unlist


  updated_trainlist <- train_list
  updated_trainlist$datalist$response <- trainresp


  nc <- nimble::nimbleCode({
    for(i in 1:n){
      response[i] ~
        dpois(Xadj[i]+pop[i]*(
          weightedV[i,1]*exp(alpha+beta*log(Z[i,1]/pop[i]*100))+
            weightedV[i,2]*exp(alpha+beta*log(Z[i,2]/pop[i]*100))+
            weightedV[i,3]*exp(alpha+beta*log(Z[i,3]/pop[i]*100)))/100);
    }
    beta ~ dflat();
    alpha ~ dflat();

  })
  constants <- updated_trainlist$constants
  datalist <- updated_trainlist$datalist
  mcmc.out <- nimble::nimbleMCMC(nc,constants = constants,data = datalist,
                                 monitors = c("beta","alpha"),
                                 inits = list(beta=beta,alpha=alpha),
                                 nburnin=5000,niter=10000,summary = TRUE,
                                 nchains = 4,samplesAsCodaMCMC = TRUE)

  list(truparam = c(alpha=alpha,beta=beta), fit = mcmc.out)
}




#'
#'
#' @examples
#'
#' inputlist <- input_multinom()
#' beta<-0.9
#' simulres_multi <- multi_simul(inputlist,beta=0.9)
#'
multi_simul <- function(inputlist,beta=0.9){
  n <- inputlist$constants$n
  nvec <- inputlist$constants$nvec
  N <- inputlist$constants$N
  X <- inputlist$constants$X

  trainY <- inputlist$datalist$Y
  for(i in 1:n){
    pp <- exp(beta*X[i,1:nvec[i]])
    trainY[i,1:nvec[i]] <- rmulti(1,prob=pp,size=N[i]);
  }

  library(nimble)
  nc_multinom <- nimbleCode({
    for(i in 1:n){
      Y[i,1:nvec[i]] ~ dmulti(pp[i,1:nvec[i]],N[i]);
    }
    for(i in 1:n){
      for(j in 1:nvec[i]){
        pp[i,j] <- exp(beta * X[i,j]);
      }
    }
    beta ~ dflat();


  })
  updated_inputlist <- inputlist
  updated_inputlist$datalist$Y <- trainY
  result_multinom <-
    nimbleMCMC(nc_multinom,constants = updated_inputlist$constants,
               data = updated_inputlist$datalist,monitors = c("beta"),
               inits = list(beta=1),nburnin=5000,niter=10000,
               summary = TRUE,WAIC=TRUE,nchains = 4,
               samplesAsCodaMCMC = TRUE)

  list(truparam = beta, fit = result_multinom)
}

eval_simul1 <- function(simulres_multi,simulres_pois){
  data.frame(param=c("beta_multi","alpha","beta"),
             coverage=c(tbetween(simulres_multi$truparam,simulres_multi$fit$summary$all.chains[4:5]),
  tbetween(simulres_pois$truparam[1],simulres_pois$fit$summary$all.chains[1,4:5]),
  tbetween(simulres_pois$truparam[2],simulres_pois$fit$summary$all.chains[2,4:5])),
  truval = c(simulres_multi$truparam,simulres_pois$truparam),
  err = c((simulres_multi$truparam-simulres_multi$fit$summary$all.chains[1])^2,
          (simulres_pois$truparam[1]-simulres_pois$fit$summary$all.chains[1,1])^2,
          (simulres_pois$truparam[2]-simulres_pois$fit$summary$all.chains[2,1])^2))
}


#'
#' @examples
#' VEdata2 <- preprocVE2()
#' bmodel <- fitstan(VEdata2)
#' mydat <- df_to_standat2(VEdata2)
#'
VE_simul <- function(mydat){
  #simuldata로 수정하기
  mu0 <- runif(1,-2,-1)
  dmu0 <- runif(1,0,1)
  tau <- runif(1,0,1)
  kappa <- runif(1,0,1)
  pkappa <- runif(1,0,1)

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

    pmu[k] <- truncnorm::rtruncnorm(1,mean=mu0+dmu0,sd = pkappa,
                                            a=mu[k]);
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
  list(truparam =
         list(mu0=mu0,dmu0=dmu0,kappa=kappa,pkappa=pkappa,tau=tau,mu=mu,pmu=pmu),
       fit=fit)

}

#'
#' @examples
#' VEdata2 <- preprocVE2()
#' bmodel <- fitstan(VEdata2)
#' mydat <- df_to_standat2(VEdata2)
#' VEres <- VE_simul(mydat)
#' eval_VEsimul(VEres$fit,VEres$truparam)
#'
eval_VEsimul <- function(fit,truparam){
  library(dplyr)
  library(rstan)
  fsummary <- rstan::summary(fit)$summary
  rnames <- rownames(fsummary)

  mucred <- fsummary[stringr::str_detect(rnames,"^mu\\["),c(4,8)]

  pmucred <- fsummary[stringr::str_detect(rnames,"^pmu\\["),c(4,8)]

  dfmu <- 1:length(truparam$mu) %>%
    purrr::map(~data.frame(param=paste0("mu",.x),coverage =
                             between(truparam$mu[.x],mucred[.x,1],mucred[.x,2]),
                           truval = truparam$mu[.x],
                           err= (truparam$mu[.x]-
                                   fsummary[stringr::str_detect(rnames,"^mu\\["),1][.x])^2 )) %>%
    do.call("rbind",.)
  dfpmu <- 1:length(truparam$pmu) %>%
    purrr::map(~data.frame(param=paste0("pmu",.x),coverage =
                             between(truparam$pmu[.x],pmucred[.x,1],pmucred[.x,2]),
                           truval = truparam$pmu[.x],
                           err= (truparam$pmu[.x]-
                                   fsummary[stringr::str_detect(rnames,"^pmu\\["),1][.x])^2 )) %>%
    do.call("rbind",.)

  rbind(dfmu,dfpmu, data.frame(
    param = c("mu0","dmu0","kappa","pkappa","tau"),
    coverage = c(tbetween(truparam$mu0,fsummary[rnames == "mu0",c(4,8)]),
    tbetween(truparam$dmu0,fsummary[rnames == "dmu0",c(4,8)]),
    tbetween(truparam$kappa,fsummary[rnames == "kappa",c(4,8)]),
    tbetween(truparam$pkappa,fsummary[rnames == "pkappa",c(4,8)]),
    tbetween(truparam$tau,fsummary[rnames == "tau",c(4,8)])),
    truval = c(truparam$mu0,truparam$dmu0,truparam$kappa,truparam$pkappa,truparam$tau),
    err = c((truparam$mu0-fsummary[rnames == "mu0",1])^2,
            (truparam$dmu0-fsummary[rnames == "dmu0",1])^2,
            (truparam$kappa-fsummary[rnames == "kappa",1])^2,
            (truparam$pkappa-fsummary[rnames == "pkappa",1])^2,
            (truparam$tau-fsummary[rnames == "tau",1])^2))
  )
}

tbetween <- function(value,vec){
  between(value,vec[1],vec[2])
}
