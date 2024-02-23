
#'
#'
#' @export
#'
#' @examples
#' datalist_doses <- load_doses()
#' delivery_table <- make_delivery()
#' agreement_table <- make_agreement()
#' df_doses <- datalist_doses[[2]][[131]]
#' check_supply(df_doses,delivery_table,agreement_table)
#' check_res <- purrr::map(datalist_doses[[2]],
#' ~check_supply(.x,delivery_table,agreement_table)) %>% unlist
#'
check_supply <- function(df_doses,delivery_table){


  vac_types <- strsplit(unique(df_doses$vaccine),split = ", ") %>%
    unlist %>% unique
  if(length(vac_types)==1){
    return("onetype")
  }
  delivery_table <- delivery_table %>% dplyr::filter(
    country==df_doses$location[1]
  )
  if(length(setdiff(vac_types,delivery_table$vaccine))==0){
    return("delivery")
  }

  return("else")

  # type : onetype, delivery, agreement, etc
}


#'
#' 국가별 reserve matrix 만들기
#'
#' @export
#'
#' @examples
#'
#' delivery_table <- make_delivery()
#' reservemat <- make_reserves(delivery_table)
#'
make_reserves <- function(delivery_table){

  #1. check_supply
  #2. type별로 reserve mat 만들기
  #3. else를 위한 reserve mat 만들기
  datalist_doses <- load_doses(onerow = TRUE)
  vacnames <- datalist_doses$vacnames

  types <- datalist_doses$datasets %>%
    purrr::map(~check_supply(.x,delivery_table)) %>% unlist

  reserve_delivery <- delivery_table %>% dplyr::group_by(country,vaccine) %>%
    dplyr::summarise(doses=sum(doses)) %>%
    dplyr::filter(vaccine %in% vacnames) %>%
    dplyr::mutate(vaccine = factor(vaccine,levels = vacnames)) %>%
    reshape::cast(country~vaccine,value="doses",add.missing = TRUE)

  total_delivery <- reserve_delivery %>% dplyr::select(-c("country")) %>% colSums(na.rm = T)

  one_df <- data.frame(country=datalist_doses$country_names[types=="onetype"]) %>%
    cbind(matrix(1,ncol=length(vacnames)))
  colnames(one_df)[-1] <- vacnames

  else_df <- data.frame(country=datalist_doses$country_names[types=="else"]) %>%
    cbind(matrix(total_delivery,ncol=length(vacnames)))
  colnames(else_df)[-1] <- vacnames


  reservemat <- rbind(one_df,else_df,reserve_delivery %>%
          dplyr::filter(country %in% datalist_doses$country_names[types=="delivery"])
  )
  reservemat[,-1] %>% purrr::map(~tidyr::replace_na(.x,0)) %>%
    data.frame(country=reservemat$country,.)


}
