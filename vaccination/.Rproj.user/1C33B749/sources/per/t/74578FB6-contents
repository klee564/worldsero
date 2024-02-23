#regression.R



#'
#'
#' @export
#'
#' @examples
#' datalist_doses <- load_doses(FALSE)
#' doses_df <- datalist_doses$datasets[[148]]
#' vacnames <- datalist_doses$vacnames
#' supply <- 1:length(vacnames)
#'
multinom_prob <- function(doses_df,supply,vacnames,beta=1){
  W <- make_dw(doses_df,vacnames,supply) %>% make_w
  data.frame(date=W$date,exp(log(W$Wmat) * beta) %>% rowwise_standard)
}

#'
#' @export
#'
#' @examples
#' datalist_doses <- load_doses(FALSE)
#' doses_df <- datalist_doses$datasets[[148]]
#' vacnames <- datalist_doses$vacnames
#' supply <- 1:length(vacnames)
#'
calcZ <- function(doses_df,supply,vacnames,del=21){
  library(lubridate)
  jind <- pastdate_index(doses_df$date,del)

  dt <- lubridate::date(doses_df$date)-lubridate::date(doses_df$date)[jind]
  dX <- doses_df$total_vaccinations - doses_df$total_vaccinations[jind]

  dW <- make_dw(doses_df,vacnames,supply)
  Wmat <- make_w(dW,del) %>% rowwise_standard

  interval_df <- make_vacinterval(vacnames)

  sort(unique(interval_df$doses)) %>%
    purrr::map(~which(interval_df$doses==.x)) %>%
    purrr::map(~(Wmat[,.x,drop=F] %*% interval_df$intervals[.x])*
                   dX/as.integer(dt)) %>%
    purrr::map(~ifelse(is.na(.x),0,.x)) %>%
    as.data.frame %>%
    magrittr::set_colnames(paste0("dose", sort(unique(interval_df$doses))))


}

#'
#' @export
#'
#' @examples
#' datalist_doses <- load_doses(FALSE)
#' doses_df <- datalist_doses$datasets[[148]]
#' vacnames <- datalist_doses$vacnames
#' supply <- 1:length(vacnames)
#'
calcv <- function(doses_df,supply,vacnames,betahat=1){

  weightmat <- multinom_prob(doses_df,supply,vacnames,beta=betahat) %>%
    dplyr::select(-c("date"))
  interval_df <- make_vacinterval(vacnames)

  sort(unique(interval_df$doses)) %>%
    purrr::map(~which(interval_df$doses==.x)) %>%
    purrr::map(~rowSums(weightmat[,.x,drop=FALSE])) %>%
    as.data.frame %>%
    magrittr::set_colnames(paste0("dose", sort(unique(interval_df$doses))))


}


#'
#'
#' @export
#'
#' @examples
#'
#' datalist_doses <- load_doses(onerow = FALSE)
#' doses_df <- datalist_doses$datasets[[51]]
#' supply <- 1:11
#' vacnames <- datalist_doses$vacnames
#' betahat <- 1
#' pop <- 1000000
#' poisson_covariate(doses_df,supply,vacnames,betahat,pop)
#'
poisson_covariate <- function(doses_df,supply,vacnames,betahat,pop,del=21){
  library(lubridate)
  jind <- pastdate_index(doses_df$date,del)

  #dt <- lubridate::date(doses_df$date)-lubridate::date(doses_df$date)[jind]
  #dX <- doses_df$total_vaccinations - doses_df$total_vaccinations[jind]
  #dW <- make_dw(doses_df,vacnames,supply)

  #vac_interval <- make_vacinterval(vacnames = vacnames)
  X <- doses_df$total_vaccinations

  Zmat <- calcZ(doses_df,supply,vacnames,del)
  Vmat <- calcv(doses_df,supply,vacnames,betahat)
  dosevec <- as.integer(gsub(names(Vmat),pattern = "dose",replacement = ""))
  #설명변수 계산


  weightedV <- Vmat*(rep(1,dim(Vmat)[1]) %*% t(1-1/dosevec))
  list(Zmat=Zmat,weightedV=weightedV,Xadj = rowSums(weightedV)*X,pop=rep(pop,length(X)),
       response = doses_df$total_vaccinations-doses_df$people_fully_vaccinated)


  #weightedT <- (make_w(dW,del) %>% rowwise_standard) %*%
  #  vac_interval$intervals


  #Xadj <- rowSums((multinom_prob(doses_df,supply,vacnames,
  #                       beta=betahat) %>%
  #           dplyr::select(-c("date","postind")))[,vac_interval$doses==2]) *
  #  doses_df$total_vaccinations

  #data.frame(Z=as.vector(weightedT*dX/as.integer(dt)),Xadj=Xadj,pop=pop,
  #           X=doses_df$total_vaccinations,
  #           Y=doses_df$people_fully_vaccinated) %>%
  #  dplyr::mutate(response = 2*(X-Y),
  #                scaledresp = 2*(X-Y)/pop*100,
  #                scaledXadj = Xadj/pop*100,
  #                scaledZ = Z/pop*100)

}




#'
#' @export
#'
#' @examples
#'
#'
#'
train_input <- function(betahat,del=21){
  datalist_doses <- load_doses(onerow = FALSE)
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
#' df2 <- df %>% dplyr::mutate(resp = (response-Xadj)/pop*100, covariate = Z/pop*100)  %>%
#' dplyr::filter(resp>0 & covariate >0)
#'
#'
train_pois <- function(betahat,del=21){

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

  train_list <- train_input(betahat,del=21)
  constants <- train_list$constants
  datalist <- train_list$datalist
  mcmc.out <- nimble::nimbleMCMC(nc,constants = constants,data = datalist,
                                 monitors = c("beta","alpha"),
                                 inits = list(beta=rep(1),alpha=rep(0)),
                                 nburnin=2000,niter=4000,summary = TRUE,
                                 nchains = 4,samplesAsCodaMCMC = TRUE)

  mcmc.out
}

