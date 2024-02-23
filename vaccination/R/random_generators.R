
#'
#'
#' @export
#'
#' @examples
#' rpois_bound(1000,ub=100)
#'
#'
rpois_bound <- function(lambda,ub,lb=0){
  lb=0
  if(ub<lb){stop("err")}
  if(lb==ub) return(lb)
  tsample <- rpois(1,lambda)
  if(tsample<=ub & tsample >= lb) return(tsample)

  supp <- max(0,lb):ub
  logprob_supp <- dpois(supp,lambda = lambda,log = T)
  sample(supp,1,prob=exp(logprob_supp - matrixStats::logSumExp(logprob_supp)))
}

#size = Xtotal[i]
#prob = probmat[i,]
#lb=lbvec
#ub=ubmat[i,]

#'
#'
#' @export
#'
#' @examples
#' size = 0
#' prob = c(0.3,0.3,0.1)
#' ub = c(10,50,60)
#' rmultinom_bound(size,prob,ub=ub)
#'
rmultinom_bound <- function(size,prob,lb=rep(0,length(prob)),ub=NULL){
  return(as.vector(rmultinom(n=1,size=size,prob=prob)))
  lb=rep(0,length(prob))
  ub=NULL
  if(size==0){
    return(rep(0,length(prob)))
  }
  if(is.null(ub)){
    ub <- rep(size,length(prob))
  }
  if(sum(ub<lb)>0){stop("err")}

  tsample <- as.vector(rmultinom(n=1,size=size,prob=prob))
  if(sum(tsample > ub) == 0 & sum(tsample < lb) == 0) return(tsample)

  ind <- which.max(prob)
  tub <- as.integer(min(ub[ind], size-sum(lb[-ind])))
  tlb <- as.integer(max(lb[ind], size-sum(ub[-ind])))
  bin_sample <- rbinom_bound(size,as.double(prob[ind]/sum(prob)),
                             lb=tlb,ub=tub)

  res <- lb
  res[ind] <-  bin_sample
  res[-ind] <-  as.vector(rmultinom_bound(size=size-bin_sample,prob=prob[-ind],
                  lb=lb[-ind],ub=ub[-ind]))
  return(res)
}

#'
#' @export
#'
#' @examples
#' rbinom_bound(10000,0.8,ub=30)
#'
#'
rbinom_bound <- function(size,prob,lb=0,ub=size){
  lb =0
  up = size
  if(prob==0){return(0)}
  if(is.na(ub)) ub <- size
  if(ub<lb) stop("err_binom")
  if(lb==ub) return(lb)
  tsample <- rbinom(1,size,prob)
  if(tsample<=ub & tsample >= lb) return(tsample)

  supp <- max(0,lb):min(ub,size)
  logprob_supp <- dbinom(supp,size,prob,log = T)
  sample(supp,1,prob=exp(logprob_supp - matrixStats::logSumExp(logprob_supp)))
}
