
VEdata2 <- preprocVE2()
VEest <- estimate_VE(VEdata2)
VEestlist <- 1:nrow(VEest$VEhalf) %>%
  purrr::map(~data.frame(`1`=VEest$VEhalf[.x,],final=VEest$VEfull[.x,]) %>%
               magrittr::set_colnames(c("1","final")))

dosedf_list = load_doses(onerow=TRUE)
Xbyvac <- load_bymanu(dosedf_list$vacnames)
dataXY <- list(dosedf_list = dosedf_list)
estimatedXY <- estimate_XY(dataXY) ## nimble package required

#sub data
datalist_doses <- load_doses(onerow=TRUE)
popmat <- load_population() %>%
  dplyr::filter(country %in% datalist_doses$country_names) %>%
  dplyr::mutate(country = factor(country,levels = datalist_doses$country_names)) %>%
  dplyr::arrange(country)
vac_interval <- make_vacinterval(datalist_doses$vacnames)
ndoses <- vac_interval$doses


ThetaV <- list()
for(cind in 1:length(datalist_doses$datasets)){
  eps <- 1e-8
  numMCMC <- 100
  betaparam <- VEestlist %>%
    purrr::map(~makethetaV(Xbyvac=estimatedXY$Xbyvac[[cind]],
                           Ybyvac=estimatedXY$Ybyvac[[cind]],
                           ndoses=ndoses,VEdf = .x/100)/popmat$pop[cind]) %>%
    do.call("rbind",.) %>% as.data.frame %>%
    purrr::map(~ifelse(.x>1, runif(numMCMC,0.99,0.999),.x)) %>%
    purrr::map(~ifelse(.x<=0, runif(numMCMC,0.00001,0.0001),.x)) %>%
    #purrr::map(~EnvStats::ebeta((.x+runif(numMCMC,0,eps))[.x<1 & .x>eps],method="mme")$parameter) %>%
    purrr::map(~EnvStats::ebeta(.x,method="mle")$parameter) %>%
    do.call("rbind",.) %>% as.data.frame %>%
    dplyr::mutate(thetaVhat = shape1/(shape1+shape2))

  ThetaV[[cind]]<-
    data.frame(country=datalist_doses$datasets[[cind]]$location,
               date= datalist_doses$datasets[[cind]]$date) %>%
    cbind(.,betaparam)
}

result <- list(ThetaV=ThetaV,VEestlist=VEestlist,estimatedXY=estimatedXY)
saveRDS(result,"thetaVreal.rds")
