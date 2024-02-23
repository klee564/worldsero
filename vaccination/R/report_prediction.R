

pred_total <- function(countrynames=NULL){
  library(dplyr)
  confirmed_list <- load_confirmed()
  datalist_doses <- load_doses(onerow=FALSE)

  if(is.null(countrynames)){
    countrynames <- datalist_doses$country_names
  }
  pred_result <- predict_overall(countrynames)


  res <- purrr::map2(countrynames,pred_result,~make_preddf(datalist_doses$datasets
                            [[which(datalist_doses$country_names==.x)]],
                            confirmed_list[[which(names(confirmed_list)==.x)]],.y))

  names(res) <- countrynames

}

#'
#'
#'
#'
#'
make_preddf <- function(doses_df,confirmed_df,effpop_df){
  summ_effpop <- effpop_df %>% dplyr::group_by(date) %>%
    dplyr::mutate(date=lubridate::date(date)) %>%
    summarise(postmean=mean(effpop),
              lower = quantile(effpop,0.025),
              upper = quantile(effpop,0.975))

  confirmed_df %>%
    dplyr::left_join(doses_df %>% dplyr::mutate(date=lubridate::date(date)) %>%
                       dplyr::select(c("date","total_vaccinations")),by="date") %>%
    dplyr::left_join(summ_effpop,by="date")
}





make_ts <- function(preddf,pop){
  library(ggplot2)
  library(patchwork) # To display 2 charts together
  library(hrbrthemes)
  popmat <- load_population()

  #국가별로 필터하고, 아래에 넣는다.
  poplist <- names(preddf) %>%
    purrr::map(~(popmat %>% dplyr::filter(country==.x))$pop)

  coeff <- 10
  # A few constants
  effpopColor <- "#69b3a2"
  confirmedColor <- rgb(0.2, 0.6, 0.9, 1)


  glist <- purrr::map2(preddf,poplist,
              function(x,y){
                tscale <- 10*max(x$confirmed)/y
                ggplot(x,aes(x=lubridate::date(date),y=postmean/y*100))+
                  geom_line(linetype="dashed",color=effpopColor) +
                  geom_ribbon(aes(ymin=lower/y*100,ymax=upper/y*100), linetype=2, alpha=0.1) +
                  geom_line(aes(y=(pmax(confirmed,0)/y)*1000/tscale),color=confirmedColor)+
                  scale_y_continuous(
                    name = "Effective vaccinated (%)",
                    # Add a second axis and specify its features
                    sec.axis = sec_axis(~.*tscale  , name="Daily confirmed per 1000")
                  )  +
                  theme_ipsum() +
                  theme(
                    axis.title.y = element_text(color = effpopColor, size=13),
                    axis.title.y.right = element_text(color = confirmedColor, size=13)
                  )
              }
              )

}

#'
#' @examples
#' countrylist <- c("US","United Kingdom","China","India","Korea, South","Israel")
#'
visual1 <- function(countrylist,glist){


  res2 <- countrylist %>% purrr::map(~
    glist[[which(names(glist)==.x)]]+ggtitle(.x)
  )
}

visual2 <- function(preddf){

  #국가별 weight를 계산한다.

  tres2 <- preddf %>% purrr::map(~summ_bymon(.x)) %>%
    purrr::map2(poplist,~.x %>% dplyr::mutate(effpop_mean = effpop_mean/.y,
                                              effpop_first = effpop_first/.y,
                                              confirmed_ratio = confirmed/.y))

  tres3 <- do.call("rbind",tres2 ) %>%
    dplyr::mutate( effpop_mean = tidyr::replace_na(effpop_mean,0)*100)

  hist(tres3$effpop_mean)
  library(ggridges)
  ggplot(tres3 ,aes(x=effpop_mean,y=mon,group=mon)) +
    geom_density_ridges(fill = "#00AFBB") +
    xlab("effective vaccinated(%)") +
    ylab("month")


}





vis_overall <- function(df){
  popmat <- load_population()

  #국가별로 필터하고, 아래에 넣는다.
  poplist <- names(df) %>% purrr::map(~(popmat %>% dplyr::filter(country==.x))$pop)


  tres2 <- df %>% purrr::map(~summ_bymon(.x)) %>%
    purrr::map2(poplist,~.x %>% dplyr::mutate(effpop_mean = effpop_mean/.y,
                                              effpop_first = effpop_first/.y,
                                              confirmed_ratio = confirmed/.y)) %>%
    purrr::map(~c(ratio=last(.x$confirmed)/first(.x$confirmed),
                  lasteff = last(.x$effpop_mean),lastconfirmed = last(.x$confirmed_ratio))) %>%
    do.call("rbind",.) %>% as.data.frame

  tres2 <- df %>% purrr::map(~summ_bymon(.x)) %>%
    purrr::map2(poplist,~.x %>% dplyr::mutate(effpop_mean = effpop_mean/.y,
                                              effpop_first = effpop_first/.y)) %>% purrr::map(~c(ratio=last(.x$confirmed)/first(.x$confirmed),
                                                                                                 lasteff = last(.x$effpop_mean))) %>%
    do.call("rbind",.) %>% as.data.frame
  return(tres2)
  #plot(tres2$lasteff,log(tres2$ratio))
}



#' 국가,달별 요약통계(관심있는 달을 넣으면 그 시점에서 effpop의 평균과 확진자수의 평균을 계산한다.)
#' @export
#'
#'
#'
summ_bymon <- function(preddf,monvec=NULL){
  lubridate::month(preddf$date)
  preddf %>% dplyr::mutate(mon = lubridate::month(date)) %>%
    dplyr::group_by(mon) %>%
    dplyr::summarise(confirmed = mean(confirmed,na.rm = T),
                      effpop_mean = mean(postmean,na.rm = T),
                     effpop_first = min(postmean,na.rm = T))
}

