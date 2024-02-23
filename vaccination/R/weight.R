
#'
#' @export
#'
#' @examples
#' datalist_doses <- load_doses()
#' vacnames <- datalist_doses$vacnames
#' doses_df <- datalist_doses$datasets[[2]]
#' supply <- 1:10
#'
#'
make_dw <- function(doses_df,vacnames,supply){

  vaccine_list <- doses_df$vaccine %>%
    purrr::map(~unlist(strsplit(.x,split = ",")) %>%
                 stringr::str_trim(.) ) %>%
    purrr::map(~is.element(vacnames,.x)) %>%
    do.call("rbind",.)

  dXsupply <- diff(c(0,doses_df$total_vaccinations)) %*%
    t(supply)
  list(dWmat = vaccine_list * dXsupply, date=doses_df$date)
}

#'
#' @export
#'
#' @examples
#' datalist_doses <- load_doses()
#' vacnames <- datalist_doses$vacnames
#' doses_df <- datalist_doses$datasets[[2]]
#' supply <- 1:10
#' dW <- make_dw(doses_df,vacnames,supply)
#'
#'
make_w <- function(dW,duration=NULL){
  if(is.null(duration)){
    return( list(Wmat = apply(dW$dWmat,2,cumsum) ,date=dW$date))
  }
  indexvec <- pastdate_index(dW$date,duration)
  1:length(indexvec) %>% purrr::map(~colSums(matrix(dW$dWmat[indexvec[.x]:.x,],
                                                    ncol=dim(dW$dWmat)[2]))) %>%
    do.call("rbind",.)

}


