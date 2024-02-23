

#'
#'
#' @export
#'
#'
#'
#'
argmax_marginal <- function(){
  vaccine_data <- make_vacclinical()
  library(nimble)
  nimblecode <- nimble::nimbleCode({
    for(i in 1:N){
      nv[i] ~ dpois(sqrt(lambda[i])*(Nv[i]*(1-ve[i])/(Nv[i]*(1-ve[i])+Nc[i])))
      nc[i] ~ dpois(sqrt(lambda[i])*(Nc[i]/(Nv[i]*(1-ve[i])+Nc[i])))
      lambda[i] ~ dhalfflat()
      ve[i] ~ dbeta(alpha,beta)

    }
    alpha ~ dhalfflat()
    beta ~ dhalfflat()
  })

  vaccine_data %>% split(vaccine_data$type) %>%
    purrr::map(function(df){
    constants <- df %>% dplyr::select(Nv,Nc) %>%
      as.list() %>% c(N=dim(df)[1])
    datalist <- df %>% dplyr::select(nv,nc) %>% as.list()
    initlist <- list(lambda=rep(1,dim(df)[1]),ve=rep(0.6,dim(df)[1]),
                     alpha=1,beta=1)

    nimblemodel <- nimbleModel(code = nimblecode, name = "vac",
                               constants = constants,data = datalist,
                               inits = initlist)
    box <- list( list(c("alpha","beta"), c(0, Inf))) ## Constraints for the parameters
    MCEMfit <- buildMCEM(model = nimblemodel, boxConstraints = box,
                         latentNodes = c(paste0("ve[1:",dim(df),"]"),
                             paste0("lambda[1:",dim(df),"]")))
    MCEMfit$run()
    }
    )

}


#'
#'
#' @export
#'
#' @examples
#' hyperparam <- list(dose1 = c(alpha=10.5,beta=4.497525),
#' dose2 = c(alpha=10.5,beta=4.497525))
#' hyperparam <- c(alpha=10.5,beta=4.497525)
#' vaccinedata <- make_vacclinical() %>% dplyr::filter(doses==2)
#' res_hier <- hierarchical_train(hyperparam)
#'
#'
hierarchical_train <- function(vaccinedata,hyperparam,vacnames,samplenum=100){
  library(nimble)
  nimblecode_vac <- nimble::nimbleCode({
    for(i in 1:N){
      nv[i] ~ dpois(sqrt(lambda[i])*(Nv[i]*(1-ve[i])/(Nv[i]*(1-ve[i])+Nc[i])))
      nc[i] ~ dpois(sqrt(lambda[i])*(Nc[i]/(Nv[i]*(1-ve[i])+Nc[i])))
      lambda[i] ~ dhalfflat()
      ve[i] ~ dbeta(alpha,beta)

    }
  })
  N = dim(vaccinedata)[1]
  constants <- c(vaccinedata %>% dplyr::select(Nv,Nc) %>% as.list(),
                 list(N=N,
                    alpha=hyperparam[["alpha"]],
                    beta=hyperparam[["beta"]]))
  datalist <-  vaccinedata %>% dplyr::select(nv,nc) %>% as.list()
  initlist <- list(lambda=rep(1,N),ve=rep(0.6,N))
  mcmc.out <- nimbleMCMC(nimblecode_vac,constants = constants,data = datalist,
             monitors = c("ve","lambda"),inits = initlist,nburnin=1000,niter=2000,
             summary = TRUE,WAIC=TRUE,nchains = 4,samplesAsCodaMCMC = TRUE)
  samplemat <- mcmc.out$samples %>% do.call("rbind",.)
  veind <- grep(pattern = "ve",colnames(samplemat))

  res <- rbeta(samplenum,hyperparam[["alpha"]],hyperparam[["beta"]]) %*%
    t(rep(1,length(vacnames)))
  tind <- vaccinedata$vaccine %>% purrr::map(~which(.x==vacnames)) %>%
    unlist()

  for(i in 1:length(tind)){
    res[,tind[i]] <- sample(samplemat[,veind[i]],samplenum)
  }
  return(res)
}


function(){
  reslist <- res_hier$res_MCEM %>%
    purrr::map(~rbeta(samplenum,.x[["alpha"]],.x[["beta"]])) %>%
    purrr::map(~.x %*% t(rep(1,length(vacnames)))) %>%
    setNames(c("VEhalf","VEfull"))


  #hierarchical_train()
  #위 사용된 것의 이름과 vacnames의 위치를 대조해서 덮어씌운다.

}


#'
#' @export
#'
#' @examples
#' res_hier <- readRDS("res_hier.rds")
#' doselist <- load_doses()
#' vacnames <- doselist$vacnames
#'
#'
vaceff_sample <- function(res_hier,vacnames,samplenum=100){
  vaccine_data <- make_vacclinical()
  vacinterval <- make_vacinterval(vacnames)

  #1. 데이터가 있는 경우
  1:2 %>% purrr::map(function(targetdose){
    samplemat <- res_hier$res_EB[[targetdose]]$samples %>% do.call("rbind",.)

    veind <- grep(pattern = "ve",colnames(samplemat))


    localeff_df <- data.frame(samplemat[sample(dim(samplemat)[1],samplenum),veind]) %>%
      magrittr::set_colnames((vaccine_data %>%
                                dplyr::filter(doses==targetdose))$vaccine) %>%
      dplyr::mutate(postind = 1:samplenum) %>%
      reshape::melt(id="postind",variable_name = "vaccine")

    globaleff_df <- data.frame(postind = 1:samplenum,
                               value = rbeta(samplenum,res_hier$res_MCEM[[targetdose]][["alpha"]],
                                             res_hier$res_MCEM[[targetdose]][["beta"]]))

    setdiff((vacinterval %>% dplyr::filter(doses>=targetdose))$vaccine,
            localeff_df$vaccine) %>%
      purrr::map(~data.frame(vaccine=.x,globaleff_df)) %>%
      bind_rows(localeff_df)
  })
}

