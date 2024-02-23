#'
#'
gen_paramXY <- function(){
  pois_coef <- data.frame(pois_alpha = runif(1,-2,0),pois_beta=runif(1,0,2))
  #multi_coef <- runif(1,-0.5,1.5)
  multi_coef <- runif(1,0.5,1.5)
  list(pois_coef=pois_coef,multi_coef=multi_coef)
}


#'
#' @examples
#' params <- gen_paramXY()
#' trueXY <- make_trueXY(params)
#'
make_trueXY <- function(params){
  miss_dflist <- load_doses(onerow=TRUE)
  miss_dflist$datasets <- miss_dflist$datasets %>%
    purrr::map(~.x %>% dplyr::mutate(people_fully_vaccinated=NA))

  vacnames <- miss_dflist$vacnames

  delivery_table <- make_delivery()
  reservemat <- make_reserves(delivery_table)
  popmat <- load_population() %>%
    dplyr::filter(country %in% miss_dflist$country_names) %>%
    dplyr::mutate(country = factor(country,levels = miss_dflist$country_names)) %>%
    dplyr::arrange(country)

  #fulldatalist_doses <-load_doses(onerow=TRUE)

  generated_XY <- 1:length( miss_dflist$datasets)%>%
    purrr::map(function(cind){
      miss_df <- miss_dflist$datasets[[cind]]
      pop <- popmat$pop[cind]
      countryname <- miss_dflist$country_names[cind]
      supply <- reservemat %>% dplyr::filter(country == countryname) %>%
        dplyr::select(-c("country")) %>% unlist %>% as.vector

      generate_XY(miss_df,pop,supply,vacnames,params)

  })
  trueXbyvac <- generated_XY %>% purrr::map(~.x$Xbyvac)
  trueYbyvac <- generated_XY %>% purrr::map(~.x$Ybyvac)
  trueY <- generated_XY %>% purrr::map(~.x$fullY)

  #for(cind in 1:length( miss_dflist$datasets)){
  #if(cind%%10==0) print(cind)
  #for(cind in 1:20){
  #  miss_df <- miss_dflist$datasets[[cind]]
  #  pop <- popmat$pop[cind]
  #  countryname <- miss_dflist$country_names[cind]
  #  supply <- reservemat %>% dplyr::filter(country == countryname) %>%
  #    dplyr::select(-c("country")) %>% unlist %>% as.vector
  #  generated <- generate_XY(miss_df,pop,supply,vacnames,params)

  #  trueY[[cind]] <- generated$fullY

  #  trueXbyvac[[cind]] <- generated$Xbyvac
  #  trueYbyvac[[cind]] <- generated$Ybyvac
  #}
  list(Y=trueY,Ybyvac=trueYbyvac,Xbyvac=trueXbyvac)
  #Xtotal을 제외하고 나머지를 missing으로 만든후에,
  #Xbyvac과 Ybyvac을 채운다.
}

#'
#' @export
#'
#' @examples
#'
#' cind <- 1
#' miss_dflist <- load_missdoses()
#' miss_df <- miss_dflist$datasets[[cind]]
#' delivery_table <- make_delivery()
#' vacnames <- datalist_doses$vacnames
#' reservemat <- make_reserves(delivery_table)
#' popmat <- load_population() %>%
#' dplyr::filter(country %in% datalist_doses$country_names) %>%
#'   dplyr::mutate(country = factor(country,levels = datalist_doses$country_names)) %>%
#'  dplyr::arrange(country)
#' pop <- popmat$pop[cind]
#' countryname <- miss_dflist$country_names[cind]
#' supply <- reservemat %>% dplyr::filter(country == countryname) %>%
#' dplyr::select(-c("country")) %>% unlist %>% as.vector
#' vacnames <- datalist_doses$vacnames
#'
#'
generate_XY <- function(doses_df,pop,supply,vacnames,params,thisXbymanu=NULL,del=21){
  pois_coef <- params$pois_coef
  multi_coef <- params$multi_coef

  countryname <- unique(doses_df$location)
  vac_interval <- make_vacinterval(vacnames)
  coefY <- rep(1,dim(doses_df)[1]) %*% t(vac_interval$doses)

  Xbyvac <- predict_Xbyvac(doses_df,supply,vacnames,
                           coef=multi_coef)

  #Ybyvac 생성
  if(dim(doses_df)[1]==1){
    fullY <- 0
    Ybyvac <- rep(0,length(vacnames))
  }else{
    fullY <- predict_Y(doses_df,pop,supply,vacnames,
                       pois_coef,
                       betahat=1,
                       Xbyvac=Xbyvac,del=del)
    Ybyvac <- predict_Ybyvac(doses_df,supply,Xbyvac,fullY,vacnames)
  }

  return(list(fullY=fullY,Ybyvac=Ybyvac,Xbyvac=Xbyvac))
}



#'
#' @examples
#' params <- gen_paramXY()
#' trueXY <- make_trueXY(params)
#' dataXY <- make_simulXY(trueXY)
#'
make_simulXY <- function(trueXY){
  #trueXY를 input으로 받고,

  dosedf_list <- load_doses(onerow=TRUE)

  for(ind in 1:length(dosedf_list$datasets)){
    nonnaind <- !is.na(dosedf_list$datasets[[ind]]$people_fully_vaccinated)
    dosedf_list$datasets[[ind]]$people_fully_vaccinated[nonnaind] <-
      trueXY$Y[[ind]][nonnaind]
  }
  list(Xbyvac =trueXY$Xbyvac,dosedf_list = dosedf_list)
}


#'
#' @examples
#' params <- gen_paramXY()
#' trueXY <- make_trueXY(params)
#' dataXY <- make_simulXY(trueXY)
#' estimatedXY <- estimate_XY(dataXY,params)
#'
#'
estimate_XY <- function(dataXY,params=NULL){
  if(is.null(params)){
    params <- list()
    params$pois_coef <- c(pois_alpha=-1.9,pois_beta=0.2)
    params$multi_coef <- 1
  }

  #coefficient 만들기
  multicoef <- simultrain_multi(dataXY$dosedf_list,dataXY$Xbyvac)
  poiscoef <- simultrain_pois(dataXY$dosedf_list,betahat = multicoef$summary$all.chains[1])
  #poiscoef3 <- simultrain_pois(dataXY$dosedf_list)


  #between(params$multi_coef,multicoef$summary$all.chains[4],multicoef$summary$all.chains[5])
  #params에 바꿔넣기
  estparams <- params
  estparams$multi_coef <- multicoef$summary$all.chains[1]
  estparams$pois_coef[1]

  estparams$multi_coef <- multicoef$summary$all.chains[,1]
  estparams$pois_coef <- poiscoef$summary$all.chains[,1]

  poiscoef$summary$all.chains[1,1]
  poiscoef$summary$all.chains[2,1]

  names(estparams$pois_coef) <- names(params$pois_coef)

  miss_dflist <- dataXY$dosedf_list

  vacnames <- miss_dflist$vacnames

  delivery_table <- make_delivery()
  reservemat <- make_reserves(delivery_table)
  popmat <- load_population() %>%
    dplyr::filter(country %in% miss_dflist$country_names) %>%
    dplyr::mutate(country = factor(country,levels = miss_dflist$country_names)) %>%
    dplyr::arrange(country)

  #Ybyvac생성하기
  #Ybyvac 생성
  dataXbyvac <- list()
  dataYbyvac <- list()
  dataY <- list()
  for(cind in 1:length( miss_dflist$datasets)){
    if(cind%%10==0) print(cind)
    #for(cind in 1:20){
    miss_df <- miss_dflist$datasets[[cind]]
    pop <- popmat$pop[cind]
    countryname <- miss_dflist$country_names[cind]
    supply <- reservemat %>% dplyr::filter(country == countryname) %>%
      dplyr::select(-c("country")) %>% unlist %>% as.vector

    generated <- generate_XY(miss_df,pop,supply,vacnames,estparams)

    dataY[[cind]] <- generated$fullY

    dataXbyvac[[cind]] <- generated$Xbyvac
    dataYbyvac[[cind]] <- generated$Ybyvac
  }
  list(Y=dataY,Ybyvac=dataYbyvac,Xbyvac=dataXbyvac,estparams=estparams)
}


