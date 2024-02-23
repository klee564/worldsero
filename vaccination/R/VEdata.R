
#'
#' @examples
#' VEdata2 <- preprocVE2()
#'
preprocVE2 <- function(){
  library(magrittr)

  VEdata2 <- system.file("extdata","VE3.csv",
                         package = "vaccination") %>% read.csv %>%
    dplyr::filter(vaccine!="AstraZeneca" | outcome!="asymptomatic")  %>%
    dplyr::filter(vaccine!="AstraZeneca" | author!="Alali") %>% VE_to_effect

  VEdata2$vaccine <-
    dplyr::recode(VEdata2$vaccine,AstraZeneca="Oxford/AstraZeneca",
                  Gamaleya="Sputnik V",Janssen = "Johnson&Johnson",
                  Pfizer = "Pfizer/BioNTech",Sinopharm = "Sinopharm/Beijing")

  return(VEdata2)
}

VE_to_effect <- function(df){

  df %>% dplyr::filter(!is.na(LCI) & !is.na(UCI)) %>%
    dplyr::mutate(UCI= ifelse(UCI==100,99.99,UCI),
                  VE = ifelse(VE==100,99.99,VE)) %>%
    dplyr::mutate(yi = log(1-VE/100),upper = log(1-LCI/100),lower = log(1-UCI/100)) %>%
    dplyr::mutate(vi = ((upper-lower)/(2*1.96))^2) %>%
    dplyr::select(vaccine,dose_number,dose,yi,vi,upper,lower,VE,LCI,UCI)
}



RR_to_VE <- function(RR){
  pmax((1-exp(RR))*100,0)
  #(1-RR)*100
}



make_VEdata <- function(vacdf){

  genrow <- function(onerow){
    vis <- runif(onerow$numdat,0,.5)
    yis <- purrr::map(vis,~rnorm(1,log(1-onerow$VE /100),sd=sqrt(.x) )) %>%unlist
    data.frame(vaccine = onerow$vacnames,dose=onerow$doses,
               yi=yis ,vi=vis  )
  }

  tdf <- vacdf %>% dplyr::filter(numdat>0)
  1:nrow(tdf) %>% purrr::map(~genrow(tdf[.x,])) %>%
    do.call("rbind",.) %>% dplyr::mutate(vaccine=as.character(vaccine))


}
