

#'
#' @export
#'
#' @examples
#'
#' datalist_doses <- load_doses()
#' bymanu <- load_bymanu(datalist_doses$vacnames)
#' doses_df <- datalist_doses$datasets[[which(datalist_doses$country_names=="US")]]
#' X_manu <- bymanu %>% dplyr::filter(location=="US")
#' reservevec <- make_reserves(delivery_table,vacnames) %>%
#' dplyr::filter(country=="US") %>% dplyr::select(-c("country")) %>% unlist
#'
#' W <- make_dw(doses_df,datalist_doses$vacnames,reservevec) %>% make_w
#' res <- covariate_multinom(W,X_manu)
#'
#'
make_XWmatrix <- function(W,X_manu){
  #1. 겹치는 날짜만 찾기
  interdate1 <- base::intersect(X_manu$date,W$date)

  #2. X_manu에 관측된 값이 W에 있는 날짜 찾기
  Xmanu_mat <- X_manu %>% dplyr::filter(date %in% interdate1) %>%
    dplyr::select(-c("location","date"))
  Wmat2 <- W$Wmat[is.element(W$date,interdate1),]

  valid_ind <- 1:dim(Xmanu_mat)[1] %>%
    purrr::map(~is.element(which(Xmanu_mat[.x,]>0),which(Wmat2[.x,]>0)) %>% mean) %>%
    unlist
  Xmanu_mat <- Xmanu_mat[valid_ind==1,]
  Wmat2 <- Wmat2[valid_ind==1,]
  return(list(Xmat=Xmanu_mat,Wmat = Wmat2))
}


proc_NAcell <- function(Xvec,Wvec){
  len <- length(Xvec)
  if(length(which(Xvec>0))<2) return(NULL)

  list(c(unlist(Xvec[which(Xvec>0)]),rep(NA,len))[1:len],
       c(Wvec[which(Xvec>0)],rep(NA,len))[1:len])
}




#'
#' @export
#'
#' @examples
#'
multinom_datamat <- function(){
  datalist_doses <- load_doses(onerow = FALSE)
  bymanu <- load_bymanu(datalist_doses$vacnames)
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
  list(response=response,covariate=covariate)
}

#'
#' @export
#'
#' @examples
#'
input_multinom <- function(){
  multidata <- multinom_datamat()
  resmat <- multidata$response
  covmat <- log(multidata$covariate)
  nvec <- apply(multidata$response,1,function(x){sum(!is.na(x))})
  constants <- list(n=dim(covmat)[1],N=rowSums(resmat,na.rm = T),
                    nvec = nvec,X=covmat)
  datalist <- list(Y=resmat)
  list(constants=constants,datalist=datalist)
}

#'
#' @export
#'
#' @examples
#'
train_multinom <- function(){
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
  inputlist <- input_multinom()
  result_multinom <-
    nimbleMCMC(nc_multinom,constants = inputlist$constants,
               data = inputlist$datalist,monitors = c("beta"),
               inits = list(beta=1),nburnin=5000,niter=10000,
               summary = TRUE,WAIC=TRUE,nchains = 4,
               samplesAsCodaMCMC = TRUE)

  result_multinom
}
#for visualize
#plot(log(resmat[,2]/resmat[,1]),covmat[,2]-covmat[,1])

