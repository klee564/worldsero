#dataload_preproc.R


#'
#' 타겟 country와 그것들의 인구수를 load한다.
#' @export
#' load_population()
#'
load_population <- function(){
  library(dplyr)
  system.file("extdata","world_data_v2.csv",
                          package = "vaccination") %>%
    read.csv() %>% na.omit() %>% dplyr::select(country,pop) %>%
    dplyr::mutate(pop = pop*1000)
}

#del_country <- c("Anguilla","Scotland","England","Wales","Bermuda",
#                 "Cayman Islands","Falkland Islands","Gibraltar",
#                 "Guernsey","Isle of Man","Jersey","Montserrat",
#                 "Northern Ireland","Saint Helena","Turks and Caicos Islands")
#del2 <- c("Aruba","Curacao","Faeroe Islands","Greenland","Hong Kong","Macao",
#          "Nauru","Palau","Tonga")

#'
#' @export
#'
#'
#'
load_countryname <- function(omitvec=NULL){
  if(is.null(omitvec)){
    omitvec <- c("Cuba","Summer Olympics 2020",
                 "Diamond Princess","MS Zaandam")
  }
  popdf <- load_population()
  library(dplyr)
  system.file("extdata","country_mapping.csv",
              package = "vaccination") %>%
    read.csv() %>% dplyr::filter(!is.element(country,omitvec)) %>%
    dplyr::filter(country %in% popdf$country)
}



#'
#'
#' @export
#'
#' @examples
#'
#' datalist_doses <- load_doses(TRUE)
#'
load_doses <- function(onerow=TRUE){
  library(dplyr)
  folder <- system.file("extdata","country_data",package = "vaccination")
  files_country <- folder %>% list.files()
  country_names <- strsplit(files_country,split=".csv") %>% unlist
  datasets <- files_country %>% purrr::map(~read.csv(paste0(folder,"/",.x))) %>%
    purrr::map(~.x %>% dplyr::filter(total_vaccinations>0) %>%
                 dplyr::mutate(vaccine = gsub(vaccine,pattern = "Wuhan|HayatVax",
                                                replacement = "Beijing")))
  country_namemap <- load_countryname()

  if(!onerow){
    oneind <- datasets %>% purrr::map(~dim(.x)[1]) %>% unlist  %>%
      `==`(.,1)
    datasets <- datasets[!oneind]
    country_names <- country_names[!oneind]
  }

  level_key <- country_namemap$country[country_namemap$vaccination!=""]
  names(level_key) <-
    country_namemap$vaccination[country_namemap$vaccination!=""]

  country_names <- recode(country_names, !!!level_key)
  ind <-which(is.element(country_names,country_namemap$country))

  setdiff(country_names,country_namemap$country)

  datasets <- datasets %>%
    purrr::map(~.x %>% dplyr::mutate(
      location = dplyr::recode(stringr::str_trim(location),!!!level_key)))

  vacnames <- datasets[ind] %>%
    purrr::map(~strsplit(unique(.x$vaccine),split = ", ") %>%
                 unlist ) %>%
    unlist %>% unique
  result <- list(country_names=country_names[ind],
       datasets = datasets[ind],
       vacnames = vacnames)
  result$datasets <-
    result$datasets %>%
    purrr::map(~.x[c(0,diff(.x$total_vaccinations))>=0,])
  return(result)
}


#datalist_doses$datasets %>%
#  purrr::map(~length(grep(.x$vaccine,pattern = "EpiVacCorona"))) %>% unlist

#'
#'
#' @export
#'
#' @examples
#'
#' datalist_doses <- load_doses()
#' bymanu <- load_bymanu(datalist_doses$vacnames)
#'
load_bymanu <- function(vacnames,nareplace=FALSE,omitvec=NULL){
  if(is.null(omitvec)){
    omitvec <- c("Netherlands")
  }

  country_namemap <- load_countryname()
  level_key <- country_namemap$country[country_namemap$vaccination!=""]
  names(level_key) <-
    country_namemap$vaccination[country_namemap$vaccination!=""]

  res <- system.file("extdata","vaccinations-by-manufacturer.csv",
              package = "vaccination") %>%
    read.csv()  %>%
    dplyr::mutate(total_vaccinations =
                    ifelse(total_vaccinations==0,NA,total_vaccinations)) %>%
    mutate(vaccine = factor(vaccine,levels = vacnames)) %>%
    split(.,.$location) %>%
    purrr::map(~.x %>%
                 tidyr::spread(key=vaccine,value=total_vaccinations,drop=FALSE)) %>%
    do.call("rbind",.) %>% dplyr::mutate(location = dplyr::recode(location,!!!level_key)) %>%
    dplyr::filter(location %in% country_namemap$country) %>%
    dplyr::filter(!is.element(location,omitvec))

  if(nareplace){
    res2 <- purrr::map(res,~ifelse(is.na(.x),0,.x)) %>%
      as.data.frame()
    return(res2 %>% split(res2$location) %>%
      purrr::map(~cbind(.x%>% dplyr::select(c("location","date")),
                        .x%>% dplyr::select(-c("location","date")) %>%
                   purrr::map(~cummax(.x)) %>% as.data.frame())) %>%
        do.call("rbind",.))

  }
  return(res)
}


#'
#'
#' @export
#'
#' @examples
#' delivery_table <- make_delivery()
#'
#'
make_delivery <- function(){


  file_delivery <- system.file("extdata","vaccine_delivery.csv",
                               package = "vaccination")

  country_namemap <- load_countryname()

  level_key <- country_namemap$country[country_namemap$delivery!=""]
  names(level_key) <-country_namemap$delivery[country_namemap$delivery!=""]

  data_df <-  read.csv(file_delivery) %>%
    mutate(vaccine = recode(vaccine,AstraZeneca="Oxford/AstraZeneca",
                            Janssen="Johnson&Johnson",
                            Sinopharm=
                              "Sinopharm/Beijing",
                            `Vaxzevria` = "Oxford/AstraZeneca",
                            `Anhui ZL` = "RBD-Dimer",
                            Covishield = "Oxford/AstraZeneca")) %>%
    dplyr::filter(vaccine != "Unknown") %>% dplyr::mutate(country = dplyr::recode(stringr::str_trim(country),
                          !!!level_key)) %>%
    dplyr::filter(country %in% country_namemap$country) %>%
    tibble %>% dplyr::group_by(country,vaccine) %>%
    dplyr::summarise(doses=sum(doses,na.rm = T))
}

#'
#' @export
#'
#' @examples
#'
#' agreement_table <- make_agreement()
#'
make_agreement <- function(){
  file_agreement <- system.file("extdata","Agreements_table.xlsx",
                               package = "vaccination")
  country_namemap <- load_countryname()
  level_key <- country_namemap$country[country_namemap$agreement!=""]
  names(level_key) <-
    country_namemap$agreement[country_namemap$agreement!=""]

  readxl::read_excel(file_agreement) %>%
    dplyr::rename(doses=`Secured doses`,
           vaccine = `Vaccine developer`,
           country = Recipient) %>%
    mutate(doses=as.integer(stringr::str_remove_all(doses,","))) %>%
    dplyr::group_by(country,vaccine) %>%
    dplyr::summarise(doses = sum(doses,na.rm = TRUE)) %>% tibble %>%
    tidyr::drop_na() %>%
    mutate(vaccine = recode(vaccine,AstraZeneca="Oxford/AstraZeneca",
                            Janssen="Johnson&Johnson",
                            `Gamaleya Research Institute` = "Sputnik V",
                            `Beijing Institute of Biological Products (CNBG)`=
                              "Sinopharm/Beijing",
                            `CanSino Biologicals` = "CanSino",
                            `Bharat Biotech` = "Covaxin"),
           country = dplyr::recode(stringr::str_trim(country),
                                   !!!level_key)) %>%
    dplyr::group_by(country,vaccine) %>%
    dplyr::summarise(doses = sum(doses,na.rm = TRUE))


}




#'
#' @export
#'
#' @examples
#'
#' vac_interval <- make_vacinterval(datalist_doses$vacnames)
#'
make_vacinterval <- function(vacnames=NULL){

  res <- system.file("extdata","vac_interval.csv",
                        package = "vaccination") %>%
    read.csv()
  if(!is.null(vacnames)){
    return(res %>% dplyr::filter(vaccine %in% vacnames) %>%
             dplyr::mutate(vaccine=factor(vaccine,levels=vacnames),
                    intervals=tidyr::replace_na(intervals,0)) %>%
      dplyr::arrange(vaccine))
  }
  return(res)

}

#'
#'
#' @export
#'
#' @examples
#' data <- make_vacclinical()
#'
#'
make_vacclinical <- function(){
  rbind(
    data.frame(vaccine="Pfizer/BioNTech",
               doses=c(1,2), Nv = 21669, Nc = 21686,
               nv = c(39,11), nc = c(82,193),type=c("half","full")) ,
    data.frame(vaccine="Sputnik V",
               doses=c(1,2), Nv = c(14999,14094), Nc = c(4950,4601),
               nv = c(30,13), nc = c(79,47),type=c("half","full")) ,
    data.frame(vaccine="Moderna",
               doses=c(1,2), Nv = c(996,13934), Nc = c(1079,13883),
               nv = c(7,5), nc = c(39,90),type=c("half","full")) ,
    data.frame(vaccine="Oxford/AstraZeneca",
               doses=c(1,2), Nv = c(9257,8597), Nc = c(9237,8581),
               nv = c(32,84), nc = c(89,248),type=c("half","full")),
    data.frame(vaccine="Johnson&Johnson",
               doses=1,Nv = 19630,Nc=19691, nv=116,nc=348,type="full")

  )
}


#'
#'
#' @export
#'
#' @examples
#' confirmed_data <- load_confirmed()
#'
load_confirmed <- function(){
  library(dplyr)
  # github에서 바로 실시간 파일을 불러올 때는 다음의 URI구조를 사용한다 : https://raw.githubusercontent.com/
  Confirmed_raw = read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
  Confirmed <- Confirmed_raw %>% dplyr::select(-c("Lat", "Long")) %>%
    reshape::melt(id.vars = c("Province.State", "Country.Region" )) %>%
    dplyr::rename(date=variable, cumconfirmed = value, country = Country.Region) %>%
    dplyr::mutate(date = as.Date(date,format = "X%m.%d.%y" ))

  conf1 <- Confirmed %>% dplyr::filter(Province.State =="") %>%
    dplyr::select(-c("Province.State"))

  Confirmed %>%
    dplyr::filter(!is.element(country,conf1$country )) %>%
    dplyr::group_by(country,date) %>%
    dplyr::summarise(cumconfirmed=sum(cumconfirmed)) %>%
    dplyr::bind_rows(conf1) %>%
    dplyr::arrange(date,country) %>%
    dplyr::mutate(country=as.factor(country)) %>%
    split(.$country) %>%
    purrr::map(~.x %>% dplyr::mutate(confirmed=c(0,diff(.x$cumconfirmed)))) %>%
    do.call("rbind",.) %>%
    dplyr::mutate(confirmed=ifelse(confirmed<0,0,confirmed))


}


