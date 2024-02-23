#reporting.R


check_vaccination <- function(){
  datalist_doses <- load_doses()

  #length(datalist_doses$country_names)
  naratio <- datalist_doses$datasets %>%
    purrr::map(~mean(is.na(.x$people_fully_vaccinated))) %>% unlist

  c(complete = sum(naratio==0),
    partial = sum(naratio<1 &naratio>0),
    allna = sum(naratio==1))
}


check_reserve <- function(){

  datalist_doses <- load_doses()
  delivery_table <- make_delivery()
  agreement_table <- make_agreement()


  delivery_check <- purrr::map2(datalist_doses$datasets %>%
                purrr::map(~strsplit(unique(.x$vaccine),split = ", ") %>%
                             unlist %>% unique),
              datalist_doses$country_names,
              ~(delivery_table %>% dplyr::filter(country==.y))$vaccine %>%
                unique %>% is.element(.x,.) %>% mean ) %>% unlist

  agreement_check <- purrr::map2(datalist_doses$datasets %>%
                                  purrr::map(~strsplit(unique(.x$vaccine),split = ", ") %>%
                                               unlist %>% unique),
                                datalist_doses$country_names,
                                ~(agreement_table %>% dplyr::filter(country==.y))$vaccine %>%
                                  unique %>% is.element(.x,.) %>% mean ) %>% unlist


  c(delivery = sum(delivery_check ==1),
    agreement =sum(delivery_check !=1 & agreement_check==1),
    nonavail = sum(delivery_check !=1 & agreement_check!=1))
}



