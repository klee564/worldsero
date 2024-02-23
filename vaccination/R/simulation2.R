
#'
#'
#' @examples
#'
#' params <- generate_param()
#'
generate_simuldata <- function(params){
  miss_dflist <- load_missdoses()
  delivery_table <- make_delivery()
  vacnames <- miss_dflist$vacnames
  reservemat <- make_reserves(delivery_table)
  popmat <- load_population() %>%
    dplyr::filter(country %in% miss_dflist$country_names) %>%
    dplyr::mutate(country = factor(country,levels = miss_dflist$country_names)) %>%
    dplyr::arrange(country)

  datalist_doses <-load_doses(onerow=TRUE)
  fulldatalist_doses <-load_doses(onerow=TRUE)
  #for(cind in 1:length( miss_dflist$datasets)){
    #if(cind%%10==0) print(cind)
  for(cind in 1:20){
    miss_df <- miss_dflist$datasets[[cind]]
    pop <- popmat$pop[cind]
    countryname <- miss_dflist$country_names[cind]
    supply <- reservemat %>% dplyr::filter(country == countryname) %>%
      dplyr::select(-c("country")) %>% unlist %>% as.vector

    generated <- generate_XY(miss_df,pop,supply,vacnames,params)
    naind <- is.na(datalist_doses$datasets[[cind]]$people_fully_vaccinated)
    fulldatalist_doses$datasets[[cind]]$people_fully_vaccinated <-
      generated$fullY
    datalist_doses$datasets[[cind]]$people_fully_vaccinated[!naind] <-
      generated$fullY[!naind]

  }

  list(vacnames=vacnames)

}



