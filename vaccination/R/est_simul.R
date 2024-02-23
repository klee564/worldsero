#'
#' @export
#'
#' @examples
#'
#'
#'
simultrain_input <- function(datalist_doses,betahat=1,del=21){
  delivery_table <- make_delivery()

  vacnames <- datalist_doses$vacnames
  reservemat <- make_reserves(delivery_table)

  popmat <- load_population() %>%
    dplyr::filter(country %in% datalist_doses$country_names) %>%
    dplyr::mutate(country = factor(country,levels = datalist_doses$country_names)) %>%
    dplyr::arrange(country)

  #supply를 국가별로
  tres <- 1:length(datalist_doses$country_names) %>%
    purrr::map(~poisson_covariate(datalist_doses$datasets[[.x]],
                                  reservemat %>% dplyr::filter(country == datalist_doses$country_names[.x]) %>%
                                    dplyr::select(-c("country")) %>% unlist %>% as.vector,
                                  vacnames,betahat,popmat$pop[.x],del=del))

  tres2 <- list(
    response = tres %>% purrr::map(~.x$response) %>% do.call("c",.),
    pop = tres %>% purrr::map(~.x$pop) %>% do.call("c",.),
    Xadj = tres %>% purrr::map(~.x$Xadj) %>% do.call("c",.),
    Zmat = tres %>% purrr::map(~.x$Zmat) %>% do.call("rbind",.),
    weightedV = tres %>% purrr::map(~.x$weightedV) %>% do.call("rbind",.)
  )

  library(nimble)
  nc <- nimble::nimbleCode({
    for(i in 1:n){
      response[i] ~
        dpois(Xadj[i]+pop[i]*(
          weightedV[i,1]*exp(alpha+beta*log(Z[i,1]/pop[i]*100))+
            weightedV[i,2]*exp(alpha+beta*log(Z[i,2]/pop[i]*100)))/100);
    }
    beta ~ dflat();
    alpha ~ dflat();

  })

  tind <- !is.na(tres2$response)
  constants <- list(n=sum(tind))

  datalist <- list(pop = tres2$pop[tind],
                   Xadj = tres2$Xadj[tind],
                   Z = tres2$Zmat[tind,],
                   weightedV = tres2$weightedV[tind,],
                   response=tres2$response[tind])
  return(list(constants=constants,datalist=datalist))
}



#'
#' @export
#'
#' @examples
#'
#' datalist_doses <- dataXY$dosedf_list
#'
simultrain_pois <- function(datalist_doses,betahat=1,del=21){

  #train_input(betahat,del=del)
  library(nimble)
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

  train_list <- simultrain_input(datalist_doses,betahat,del=21)
  constants <- train_list$constants
  datalist <- train_list$datalist
  mcmc.out <- nimble::nimbleMCMC(nc,constants = constants,data = datalist,
                                 monitors = c("beta","alpha"),
                                 inits = list(beta=rep(1),alpha=rep(0)),
                                 nburnin=2000,niter=4000,summary = TRUE,
                                 nchains = 4,samplesAsCodaMCMC = TRUE)

  mcmc.out
}

