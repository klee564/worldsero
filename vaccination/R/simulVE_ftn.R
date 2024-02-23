


#'
#'
#'
gen_paramVE <- function(numvac=12){
  purrr::rerun(numvac,sort(runif(2,50,100))) %>%
    do.call("rbind",.) %>% as.data.frame %>%
    magrittr::set_colnames(c("1","final"))
}

#'
#' yi,vi를 생성한다.
#' @examples
#' trueVE <- gen_paramVE(numvac=12)
#' VEdata <- make_dataVE(trueVE)
#'
#'
make_dataVE <- function(trueVE){
  VEdata2 <- preprocVE2()
  datalist_doses <- load_doses(onerow=FALSE)


  #reshape
  numdats <- table(VEdata2$vaccine,VEdata2$dose) %>% as.data.frame %>%
    magrittr::set_colnames(c("vacnames","doses","numdat"))

  #true VE를 붙인다.
  vacdf <- expand.grid(vacnames=datalist_doses$vacnames,doses=c("1","final")) %>%
    dplyr::left_join(numdats,by=c("vacnames"="vacnames","doses"="doses")) %>%
    dplyr::mutate(numdat = tidyr::replace_na(numdat,0)) %>%
    data.frame(VE=reshape::melt(trueVE) %>%
                 magrittr::set_colnames(c("doses","VE")) %>% dplyr::pull(VE))
  data.frame()

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




#'
#'
#' @examples
#' trueVE <- gen_paramVE(numvac=12)
#' VEdata <- make_dataVE(trueVE)
#' VEest <- estimate_VE(VEdata)
#'
#'
estimate_VE <- function(VEdata){
  library(rstan)
  library(posterior)
  bmodel <- fitstan(VEdata)
  datalist_doses <- load_doses(onerow=FALSE)

  postsample2 <- bmodel$postsamples

  dfbox <- rbind(
    1:7 %>%
      purrr::map(~data.frame(vind=.x, type= "fully vaccinated",
                             VE=RR_to_VE(postsample2[[paste0("mu[",.x,"]")]]))) %>%
      do.call("rbind",.),
    1:7 %>%
      purrr::map(~data.frame(vind=.x, type= "partially vaccinated",
                             VE=RR_to_VE(postsample2[[paste0("pmu[",.x,"]")]]))) %>%
      do.call("rbind",.)) %>%
    dplyr::left_join(data.frame(vind=1:length(bmodel$vacnames),vaccine=bmodel$vacname),
                     by="vind")#%>%dplyr::filter(!(type=="partially vaccinated" & vaccine=="Janssen" ))






  overallVE <- rbind(data.frame(type="fully vaccinated",
                                VE=purrr::map2(postsample2$mu0,postsample2$kappa,
                                               ~rnorm(1,.x,.y)) %>% unlist %>% RR_to_VE() ),
                     data.frame(type="partially vaccinated",
                                VE=purrr::map2(postsample2$mu0+postsample2$dmu0,postsample2$pkappa,
                                               ~rnorm(1,.x,.y)) %>% unlist %>% RR_to_VE() ))

  mcmc_df <- list(VEhalf=c(),VEhalf=c())
  for(i in 1:length(datalist_doses$vacnames)){
    if(is.element(datalist_doses$vacnames[i],dfbox$vaccine)){
      mcmc_df$VEhalf <- dfbox %>%
        dplyr::filter(vaccine==datalist_doses$vacnames[i]&
                        type=="partially vaccinated") %>%
        dplyr::pull(VE) %>% sample(size=100) %>% cbind(mcmc_df$VEhalf,.)
      mcmc_df$VEfull <- dfbox %>%
        dplyr::filter(vaccine==datalist_doses$vacnames[i]&
                        type=="fully vaccinated") %>%
        dplyr::pull(VE) %>% sample(size=100) %>% cbind(mcmc_df$VEfull,.)
    } else{
      mcmc_df$VEhalf <- overallVE %>%
        dplyr::filter(type=="partially vaccinated") %>%
        dplyr::pull(VE) %>% sample(size=100) %>% cbind(mcmc_df$VEhalf,.)
      mcmc_df$VEfull <- overallVE %>%
        dplyr::filter(type=="fully vaccinated") %>%
        dplyr::pull(VE) %>% sample(size=100) %>% cbind(mcmc_df$VEfull,.)
    }

  }
  return(mcmc_df)

  #output은 백신별로 VE sample.
}


fitstan <- function(data){
  mydat <- df_to_standat2(data)
  fit <- stan(file = "hierarchicalmeta2.stan",
              data = mydat, iter = 8000, warmup = 5000, chains = 2)
  list(vacnames = mydat$vacname,postsamples =as_draws_df(fit), fit=fit)
}

#'
#' @examples
#' df <- preprocVE()
#' df_to_standat2(df)
#'
#'
df_to_standat2 <- function(df){
  spliteddf <- df %>%
    dplyr::mutate(vaccine= as.factor(vaccine)) %>%
    split(.$dose)

  df1 <- spliteddf[["1"]]
  df2 <- spliteddf[["final"]]
  list(K=length(unique(df$vaccine)), J1=nrow(df1),
       y1=df1$yi,v1=df1$vi,g1=as.integer(df1$vaccine),
       J2=nrow(df2),y2=df2$yi,v2=df2$vi,g2=as.integer(df2$vaccine),
       vacname = levels(df2$vaccine))

}


