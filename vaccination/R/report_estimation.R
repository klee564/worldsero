
multinom_scatter <- function(){


}


#'
#' @export
#'
#' @examples
#' readRDS("res_multinom.rds")
#'
#'
multinom_mcmc <- function(res_multinom){
  library(ggplot2)

   sampledf <- 1:length(res_multinom$samples)%>%
    purrr::map(~data.frame(res_multinom$samples[[.x]],chain=.x,
                           iter = 1:dim(res_multinom$samples[[.x]])[1]))  %>%
     do.call("rbind",.) %>%
     dplyr::mutate(chain=as.factor(chain))

  g1 <- ggplot(sampledf,aes(y=beta,x=iter)) + geom_line(aes(col=chain))

  c(mean=mean(sampledf$beta),sd=sd(sampledf$beta),
    quantile(sampledf$beta,probs = c(0.025,0.975)))
}

pois_mcmc <- function(res_pois){

}

vaccine_visual <- function(hyperparam,res_hier){
  library(reshape2)
  res1 <- res_hier[[1]]
  hyper1 <- hyperparam[[1]]

  veind <- grep(pattern = "ve",colnames(res1$samples$chain1))

  sampledf <- (res1$samples %>% do.call("rbind",.))[,veind] %>%
    data.frame(.,general=rbeta(dim(.)[1],hyper1[["alpha"]],hyper1[["beta"]]),
               id=1:dim(.)[1])

  gres <- melt(sampledf, id.vars=c("id")) %>%
    ggplot(aes(x=as.factor(variable), y=value)) +
    geom_boxplot(fill="slateblue", alpha=0.2)

  #apply(sampledf,2,function(x){c(mean=mean(x),sd=sd(x),quantile(x,probs = c(0.025,0.975)))})

}

