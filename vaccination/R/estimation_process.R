
#'
#' @export
#'
#'
multinom_estimation <- function(filename){
  res_multinom <- train_multinom()
  saveRDS(res_multinom,paste0(filename,".rds"))
}

pois_estimation <- function(filename){
  res_pois <- train_pois()
  saveRDS(res_pois,paste0(filename,".rds"))

}

vaceff_estimation <- function(){

  res_MCEM <-  argmax_marginal()

  res_EB <- hierarchical_train(hyperparam)
  saveRDS(list(res_MCEM=res_MCEM,res_EB=res_EB),"res_hier.rds")
}

