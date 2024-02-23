#predictive.R

#'
#' @export
#'
#' @examples
#'
#' datalist_doses <- load_doses(FALSE)
#' beta_sample <- runif(4)
#' supply <- 1:10
#' Xbymanu <- load_bymanu(datalist_doses$vacnames,TRUE)
#' thisXbymanu <- Xbymanu %>% dplyr::filter(location == "US")
#' doses_df <- datalist_doses$datasets[[148]]
#' Xtotal <- doses_df$total_vaccinations
#' datevec <- doses_df$date
#'
#' vacnames <- datalist_doses$vacnames
#' ubmat <- upper_Xbyvac(date,thisXbymanu)
#'
predict_Xbyvac <- function(doses_df,supply,vacnames,coef=1,
                           thisXbymanu=NULL){
  library(zoo)
  Xtotal <- doses_df$total_vaccinations
  datevec <- doses_df$date
  prob_df <- multinom_prob(doses_df,supply,vacnames,beta=coef)
  probmat <- prob_df %>% dplyr::select(-c("date"))
  vac_interval <- make_vacinterval(vacnames)

  if(!is.null(thisXbymanu)){
    ubmat <- upper_Xbyvac(datevec,thisXbymanu)
  } else{ubmat<-NULL}

  if(length(Xtotal)==1){
    return(as.vector(rmultinom(1,Xtotal,probmat[,1])))
  }

  dose1_ind <- which(vac_interval$doses==1)
  dose2_ind <- which(vac_interval$doses>1)
  tlbvec <- pmax(2*tidyr::replace_na(na.locf(doses_df$people_fully_vaccinated,
                    na.rm=FALSE),replace=0) -Xtotal,0)
  tlbvec <- lbfun(tlbvec,Xtotal)
  ubvec <- na.locf(doses_df$people_fully_vaccinated,
                  na.rm=FALSE,fromLast = TRUE)
  prob1 <- rowSums(probmat[,dose1_ind])
  lbval <- 0
  X1vec <- Xtotal
  preXval <- 0
  for(i in 1:length(Xtotal)){
    X1vec[i] <- rbinom_bound(size=Xtotal[i],prob1[i],
                             lb=max(lbval,tlbvec[i]),
                             ub=min(ubvec[i],preXval+Xtotal[i]))
    lbval <- X1vec[i]
    preXval <- X1vec[i] - Xtotal[i]
  }

  resmat <- probmat
  res1mat <- probmat[,dose1_ind]
  lbvec <- rep(0,length(vacnames))
  lbvec1 <- rep(0,length(dose1_ind))
  for(i in 1:length(Xtotal)){
    #if(is.element(datevec[i],thisXbymanu$date)){
    #    resmat[i,] <- thisXbymanu %>%
    #      dplyr::filter(date %in% datevec[i]) %>%
    #      dplyr::select(-c("date","location")) %>% as.vector
    #}else{
    sampleind <- base::intersect(which(probmat[i,]>0),dose2_ind)


    #resmat[i,sampleind] <-
    #      rmultinom_bound(size=Xtotal[i]-X1vec[i],
    #                      prob=as.vector(probmat[i,sampleind]),
    #                      lb=as.vector(lbvec)[sampleind],
    #                      ub=as.vector(ubmat[i,sampleind]))
    if(length(sampleind)==0) {
      resmat[i,] <- 0
    }else{
      resmat[i,sampleind] <- rmultinom(1,size=Xtotal[i]-X1vec[i],
                                       prob=pmax(as.vector(probmat[i,sampleind]),0.00001))
      resmat[i,-sampleind] <- 0
      lbvec <- resmat[i,]

    }

    #res1mat[i,] <-
    #  rmultinom_bound(size=X1vec[i],prob=as.vector(probmat[i,dose1_ind]),
    #                lb=lbvec1,ub=as.vector(ubmat[i,dose1_ind]))
    res1mat[i,] <-
      rmultinom(1,size=X1vec[i],prob=pmax(as.vector(probmat[i,dose1_ind]),0.00001))


    lbvec1 <- res1mat[i,]
    #}
  }
  resmat[dose1_ind] <- res1mat



  return(resmat)
}

#'
#'
#' @export
#'
#' @examples
#' thisXbymanu <- Xbymanu %>% dplyr::filter(location == "US")
#' doses_df <- datalist_doses$datasets[[148]]
#' date <- doses_df$date
#'
#'
upper_Xbyvac <- function(date,thisXbymanu){
  data.frame(date=date) %>%
    dplyr::left_join(thisXbymanu,by="date") %>%
    dplyr::select(-c("date","location")) %>%
    purrr::map(~c(zoo::na.locf(.x,fromLast=TRUE,na.rm=FALSE)[-1],Inf)) %>%
    purrr::map(~ifelse(is.na(.x),Inf,.x)) %>%
    as.data.frame
}


#'
#' @export
#'
#' @examples
#' pois_coef <- data.frame(alpha = -0.26,beta = 0.78)
#' fullY <- predict_pois(doses_df,supply,vacnames,betas_sample,betahat=1)
#'
#'
#'
predict_Y <- function(doses_df,pop,supply,vacnames,
                      pois_coef,Xbyvac,betahat=1,del=21){
  vac_interval <- make_vacinterval(vacnames)
  Ytotal <- doses_df$people_fully_vaccinated
  X <- doses_df$total_vaccinations
  resvec <- Ytotal

  if(sum(is.na(resvec))==0) return(resvec)

  Zmat <- calcZ(doses_df,supply,vacnames,del)
  Vmat <- calcv(doses_df,supply,vacnames,betahat)
  dosevec <- as.integer(gsub(names(Vmat),pattern = "dose",replacement = ""))
  #설명변수 계산

  weightedV <- Vmat*(rep(1,dim(Vmat)[1]) %*% t(1-1/dosevec))
  poismean <-
    ((pop*exp(pois_coef[["pois_alpha"]] + pois_coef[["pois_beta"]] *log(Zmat/pop*100))/100) *
      weightedV) %>% rowSums + rowSums(weightedV)*X
  #upperbound 계산


  upper1 <- ((rep(1,dim(Xbyvac)[1]) %*%
      t(1/vac_interval$doses)) *
      as.data.frame(Xbyvac))%>% apply(1,function(x){sum(as.integer(x))})

  ubvec <-
    pmin(c(zoo::na.locf(Ytotal,na.rm=FALSE,fromLast=TRUE)[-1],Inf),
         upper1,na.rm = T) #####이부분이 실제 관측값 기준으로 바꾸기

  #addlbvec <- c(0,diff(rowSums(Xbyvac[,which(vac_interval$doses==1)])))
  tlbvec <- rowSums(Xbyvac[,which(vac_interval$doses==1)])
  lbval = 0
  for(i in 1:length(Ytotal)){
    if(is.na(Ytotal[i])){
      resvec[i] <- X[i] -
        rpois_bound(lambda = poismean[i],ub=X[i] -lbval,lb=X[i] -ubvec[i])
    }
    lbval <- max(resvec[i],tlbvec[i])
  }
  return(resvec)
}

#'
#' @export
#'
#' @examples
#' datalist_doses <- load_doses(onerow=FALSE)
#' delivery_table <- make_delivery()
#' agreement_table <- make_agreement()
#' vacnames <- datalist_doses$vacnames
#' reservemat <- make_reserves(delivery_table,agreement_table)
#' popmat <- load_population() %>%
#'  dplyr::filter(country %in% datalist_doses$country_names) %>%
#'  dplyr::mutate(country = factor(country,levels = datalist_doses$country_names)) %>%
#'  dplyr::arrange(country)
#'  i<-148
#'  countryname <- datalist_doses$country_names[i]
#' doses_df <- datalist_doses$datasets[[i]]
#' supply <- reservemat %>% dplyr::filter(country == countryname) %>%
#' dplyr::select(-c("country")) %>% unlist %>% as.vector
#' pop <- popmat %>% dplyr::filter(country ==countryname) %>% dplyr::pull(pop)
#' Xbyvac <- predict_Xbyvac(doses_df,supply,vacnames,coef=1,thisXbymanu=thisXbymanu)
#' fullY <- predict_Y(doses_df,pop,supply,vacnames,pois_coef,betahat=1,del=21)
#' datevec <- doses_df$date
#'
predict_Ybyvac <- function(doses_df,supply,Xbyvac,fullY,vacnames){

  vac_interval <- make_vacinterval(vacnames)
  prob_df <- multinom_prob(doses_df,supply,vacnames,beta=1) %>%
    dplyr::select(-c("date"))
  index_past <- vac_interval$intervals %>%
    purrr::map(~pastdate_index(doses_df$date,.x))

  weightmat <-
    purrr::map2(index_past[which(vac_interval$doses>1)],
              as.data.frame(Xbyvac)[which(vac_interval$doses>1)],
              ~.y[.x]) %>%
    purrr::map2(vac_interval$doses[vac_interval$doses>1],
                ~.x/.y) %>% do.call("cbind",.) %>%
    rowwise_standard() *0.99 +
    0.01 *prob_df[,which(vac_interval$doses>1)]
  ubmat <- as.data.frame((rep(1,length(fullY)) %*%
    t(1/vac_interval$doses[which(vac_interval$doses>1)])) *
    Xbyvac[,which(vac_interval$doses>1)]) %>% purrr::map(~ceiling(.x)) %>%
    data.frame
  lbvec <-rep(0,dim(ubmat)[2])

  resmat <- Xbyvac
  for(i in 1:length(fullY)){
    #if(sum(lbvec)>(fullY[i]-sum(resmat[i,which(vac_interval$doses==1)]))){
    #  resmat[i,which(vac_interval$doses==1)] <-
    #    as.integer(resmat[i,which(vac_interval$doses==1)] *
    #    (fullY[i]-sum(lbvec))/sum(resmat[i,which(vac_interval$doses==1)]))
    #}
    resmat[i,which(vac_interval$doses>1)] <-
      rmultinom(1,
          size=max(0,fullY[i]-sum(resmat[i,which(vac_interval$doses==1)])),
          prob = pmax(weightmat[i,],0.00001))
      #rmultinom_bound(
      #  min(fullY[i]-sum(resmat[i,which(vac_interval$doses==1)]),sum(ubmat[i,])),
      #              weightmat[i,],ub=ubmat[i,],lb=lbvec)
    lbvec <- resmat[i,which(vac_interval$doses>1)]
  }

  return(resmat)

}




