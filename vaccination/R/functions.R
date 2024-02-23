

#'
#' @export
#'
#' @examples
#' doses_df <- datalist_doses$datasets[[12]]
#' datevec <- doses_df$date
#'
#'
pastdate_index <- function(datevec,del=21){
  library(lubridate)
  outer(lubridate::date(datevec),lubridate::date(datevec),`-`) %>%
    apply(2,function(x){which.min(abs(x[x<0]+del))[1]}) %>%
    tidyr::replace_na(1)

}

#'
#' @export
#'
rowwise_standard <- function(mtx){
  if(!is.matrix(mtx)){
    return(tidyr::replace_na(mtx/sum(mtx),replace=0))
  }
  data.frame(mtx/outer(rowSums(mtx),rep(1,dim(mtx)[2]))) %>%
    purrr::map(~tidyr::replace_na(.x,replace=0)) %>% do.call("cbind",.)
}





lbfun <- function(lbvec,Xtotal){
  res <- lbvec
  for(ind in length(lbvec):2){
    tvalue <- lbvec[ind]
    tvector <- diff(Xtotal[1:ind])

    res[1:(ind-1)] <-
      pmax(res[1:(ind-1)],
           pmax(tvalue - rev(cumsum(rev(tvector))),0))
  }
  return(res)
}
