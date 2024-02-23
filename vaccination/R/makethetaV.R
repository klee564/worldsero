
#'
#' @examples
#' Xbyvac <- trueXY$trueXbyvac[[1]]
#' Ybyvac <- trueXY$trueYbyvac[[1]]
#' datalist_doses <- load_doses(onerow=FALSE)
#' vac_interval <- make_vacinterval(datalist_doses$vacnames)
#' ndoses <- vac_interval$doses
#' VEdf <- trueVE/100
#'
#'
makethetaV <- function(Xbyvac,Ybyvac,ndoses,VEdf){
  coefY <- rep(1,nrow(Xbyvac)) %*% t(ndoses)

  halfpop <- (Xbyvac-coefY*Ybyvac) %>%
    purrr::map(~ifelse(.x<0,0,.x)) %>% as.data.frame

  efffull <- purrr::map2(as.data.frame(Ybyvac),VEdf$final,
                         function(x,y){
                           unlist(purrr::map(x,~rbinom(1,.x,y)) )
                         }) %>% as.data.frame() %>% rowSums
  effhalf <- purrr::map2(as.data.frame(halfpop),VEdf$`1`,
                         function(x,y){
                           unlist(purrr::map(x,~rbinom(1,.x,y)))
                         })%>% as.data.frame() %>% rowSums
  cummax(efffull+effhalf)
}
