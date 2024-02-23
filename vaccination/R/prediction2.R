#'
#' @export
#'
#' @examples
#' datalist_doses <- load_doses(onerow=FALSE)
#' res_multinom <- readRDS("res_multinom.rds")
#' res_pois <- readRDS("res_pois.rds")
#' res_hier <- readRDS("res_hier.rds")
#' vacnames <- datalist_doses$vacnames
#' mcmc_df <- make_mcmcdf(res_multinom,res_pois,res_hier,vacnames,samplenum=3)
#'
make_mcmcdf <- function(res_multinom,res_pois,res_eff,vacnames,samplenum=100){

  coef <- res_pois$samples %>% do.call("rbind",.) %>%
    as.data.frame %>% purrr::map(~sample(.x,size=samplenum)) %>%
    as.data.frame %>% magrittr::set_colnames(c("pois_alpha","pois_beta")) %>%
    cbind(data.frame(multi_beta=res_multinom$samples %>%
                       do.call("rbind",.) %>% as.vector %>%
                       sample(size=samplenum)))

  tind <- sample(1:dim(res_eff$VEhalf)[1],replace = TRUE,samplenum)
  c(list(coef=coef),res_eff %>% purrr::map(~.x[tind,]))
}


#'
#' @export
#'
#' @examples
#' datalist_doses <- load_doses(onerow=FALSE)
#' Xbymanu <- load_bymanu(datalist_doses$vacnames)
#' delivery_table <- make_delivery()
#' agreement_table <- make_agreement()
#' reservemat <- make_reserves(delivery_table,agreement_table)
#' popmat <- load_population() %>%
#' dplyr::filter(country %in% datalist_doses$country_names) %>%
#'   dplyr::mutate(country = factor(country,levels = datalist_doses$country_names)) %>%
#'  dplyr::arrange(country)
#' pop <- popmat$pop[130]
#' doses_df <- datalist_doses$datasets[[130]]
#' supply <- reservemat %>% dplyr::filter(country==datalist_doses$country_names[130]) %>%
#' dplyr::select(-c("country")) %>% unlist %>% as.vector
#' vacnames <- datalist_doses$vacnames
#'
predict_effpop <- function(doses_df,pop,supply,vacnames,thisXbymanu=NULL,
                           mcmc_df,del=21){
  countryname <- unique(doses_df$location)
  vac_interval <- make_vacinterval(vacnames)
  coefY <- rep(1,dim(doses_df)[1]) %*% t(vac_interval$doses)

  1:dim(mcmc_df$coef)[1] %>% purrr::map(function(i){
    Xbyvac <- predict_Xbyvac(doses_df,supply,vacnames,
                             coef=mcmc_df$coef$multi_beta[i],
                             thisXbymanu=thisXbymanu)
    if(dim(doses_df)[1]==1){
      Ybyvac <- rep(0,length(vacnames))
    }else{
      fullY <- predict_Y(doses_df,pop,supply,vacnames,
                         mcmc_df$coef[i,],
                         betahat=mean(mcmc_df$coef$multi_beta),
                         Xbyvac=Xbyvac,del=del)
      Ybyvac <- predict_Ybyvac(doses_df,supply,Xbyvac,fullY,vacnames)
    }

    halfpop <- (Xbyvac-coefY*Ybyvac) %>%
      purrr::map(~ifelse(.x<0,0,.x)) %>% as.data.frame
    VEhalf <- mcmc_df$VEhalf[i,]
    VEfull <- mcmc_df$VEfull[i,]

    efffull <- purrr::map2(as.data.frame(Ybyvac),VEfull,
                           function(x,y){
                             unlist(purrr::map(x,~rbinom(1,.x,y)) )
                           }) %>% as.data.frame() %>% rowSums
    effhalf <- purrr::map2(as.data.frame(halfpop),VEhalf,
                           function(x,y){
                             unlist(purrr::map(x,~rbinom(1,.x,y)))
                           })%>% as.data.frame() %>% rowSums
    data.frame(date=doses_df$date,postind = i, effpop = cummax(efffull+effhalf))
  }) %>% do.call("rbind",.)

}




#'
#' @export
#'
#' @examples
#'
#' res_multinom <- readRDS("res_multinom.rds")
#' res_pois <- readRDS("res_pois.rds")
#' res_hier <- readRDS("res_hier.rds")
#' pois_coef <- res_pois$samples$chain1[1,]
#'
#'
predict_effvac <- function(multinom_coef,pois_coef,res_eff,betahat,
                           samplenum=2,del=21){

  datalist_doses <- load_doses(onerow=TRUE)
  Xbymanu <- load_bymanu(datalist_doses$vacnames,nareplace = TRUE)
  delivery_table <- make_delivery()
  vacnames <- datalist_doses$vacnames
  reservemat <- make_reserves(delivery_table)
  popmat <- load_population() %>%
    dplyr::filter(country %in% datalist_doses$country_names) %>%
    dplyr::mutate(country = factor(country,levels = datalist_doses$country_names)) %>%
    dplyr::arrange(country)

  mcmc_df <- make_mcmcdf(multinom_coef,pois_coef,res_eff,
                         vacnames,samplenum=samplenum)

  #library(parallel)
  #detectCores()

  library("doFuture")
  library(foreach)
  library(doRNG)
  registerDoFuture()
  plan(multisession,workers=14)

  length(datalist_doses$datasets)
  ##134,141,72
  #for(i in 151:182){
  #  print(i)
  indvec <- 1:length(datalist_doses$datasets)
  res <- foreach(i = indvec) %dorng% {
    countryname <- datalist_doses$country_names[i]
    doses_df <- datalist_doses$datasets[[i]]
    supply <- reservemat %>% dplyr::filter(country == countryname) %>%
      dplyr::select(-c("country")) %>% unlist %>% as.vector

    if(is.element(countryname,Xbymanu$location)){
     # thisXbymanu <- Xbymanu %>% dplyr::filter(location==countryname)
    }else{thisXbymanu <- NULL}
    thisXbymanu <- NULL
    pop <- popmat %>% dplyr::filter(country ==countryname) %>%
      dplyr::pull(pop)

    predict_effpop(doses_df,pop,supply,vacnames,
                   thisXbymanu=thisXbymanu,mcmc_df,del=del)


  }
  names(res) <- datalist_doses$country_names[indvec]
  return(res)

}
