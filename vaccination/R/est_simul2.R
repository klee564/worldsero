

#'
#'
#' @examples
#' tdata <- simultrain_inputmulti(dataXY$dosedf_list)
#'
simultrain_inputmulti <- function(datalist_doses,Xbyvac=NULL){

  if(is.null(Xbyvac)){
    bymanu <- load_bymanu(datalist_doses$vacnames)
  } else{
    tbymanu <- load_bymanu(datalist_doses$vacnames)
    ttbymanu <- purrr::map2(datalist_doses$datasets,Xbyvac,~.x %>%
      dplyr::select(location,date) %>% data.frame(.,.y)) %>%
      do.call("rbind",.)
    names(ttbymanu) <- names(tbymanu)
    bymanu <- tbymanu %>% dplyr::select(location,date) %>%
      dplyr::left_join(ttbymanu,by=c("location"="location","date"="date"))

  }





  delivery_table <- make_delivery()
  vacnames <- datalist_doses$vacnames
  reservemat <- make_reserves(delivery_table)
  reservemat <- reservemat[,-1] %>% purrr::map(~tidyr::replace_na(.x,0)) %>%
    data.frame(country=reservemat$country,.)

  #국가의 목록


  resmat_list <- unique(bymanu$location) %>%
    purrr::map(function(countryname){
      reservevec <- reservemat %>% dplyr::filter(country ==countryname) %>%
        dplyr::select(-c("country")) %>% unlist
      doses_df <- datalist_doses$datasets[[which(datalist_doses$country_names ==countryname)]]
      W <- make_dw(doses_df,datalist_doses$vacnames,reservevec) %>% make_w
      tmp <- make_XWmatrix(W,bymanu %>% dplyr::filter(location == countryname))

      ttt <- 1:dim(tmp[[1]])[1] %>%
        purrr::map(~proc_NAcell(tmp[[1]][.x,],tmp[[2]][.x,]))

      list(response=ttt %>% purrr::map(~.x[[1]]) %>% do.call("rbind",.),
           covariate = ttt %>% purrr::map(~.x[[2]]) %>% do.call("rbind",.))
    }
    )
  response <- resmat_list %>% purrr::map(~.x[[1]]) %>% do.call("rbind",.)
  covariate <- resmat_list %>% purrr::map(~.x[[2]]) %>% do.call("rbind",.)

  resmat <- response
  covmat <- log(covariate)
  nvec <- apply(response,1,function(x){sum(!is.na(x))})
  constants <- list(n=dim(covmat)[1],N=rowSums(resmat,na.rm = T),
                    nvec = nvec,X=covmat)
  datalist <- list(Y=resmat)
  list(constants=constants,datalist=datalist)

}

#'
#' @examples
#'
#' datalist_doses <- dataXY$dosedf_list
#'
simultrain_multi <- function(datalist_doses,Xbyvac=NULL){
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

  inputlist <- simultrain_inputmulti(datalist_doses,Xbyvac)
  result_multinom <-
    nimbleMCMC(nc_multinom,constants = inputlist$constants,
               data = inputlist$datalist,monitors = c("beta"),
               inits = list(beta=1),nburnin=5000,niter=10000,
               summary = TRUE,WAIC=TRUE,nchains = 4,
               samplesAsCodaMCMC = TRUE)

  result_multinom

}
