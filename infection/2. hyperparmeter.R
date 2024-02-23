library(nimble)

### 2.1 Load data
load("Preprocessing_data.RData")

### 2.2 Nimble code
nc <- nimbleCode({
  
  for(i in 1:N){
    X[i] ~ dbin(pp[index[i]]*(theta[i])+
                  (1-pn[index[i]])*(1-(theta[i])),S[i]);
    theta[i] ~ dunif(0,1);
  }
  
  for(i in 1:M){
    clinical[i,1] ~ dbin(pp[i],tp_fn[i]);
    logit(pp[i]) ~ T(dnorm(mu_pp,sd=sig_pp),0,Inf);
    
    clinical[i,3] ~ dbin(pn[i],tn_fp[i]);
    logit(pn[i]) ~ T(dnorm(mu_pn,sd=sig_pn),0,Inf);
    
  }
  
  sig_pp ~ T(dnorm(0,sd=1),0,Inf);
  sig_pn ~ T(dnorm(0,sd=1),0,Inf);
  
  mu_pp ~ dnorm(4,sd=2);
  mu_pn ~ dnorm(4,sd=2);
  
})

### 2.3 Parameter of nimble codes
tdf_subset<-tdf%>%
  filter(!is.na(type_no2)&type_no2%in%tclinical$Type_no2)

K <- max(tdf$country_index)
M <- dim(tclinical)[1]
S_sub <- tdf_subset$Serosurvey.sample.size
X_sub <- as.integer(tdf_subset$Test.positive.size)
clinical=tclinical[,c(-1,-6)]

constants <- list(N=dim(tdf_subset)[1],M=M,S=S_sub,
                  tp_fn = clinical[,1]+clinical[,2],
                  tn_fp = clinical[,3]+clinical[,4],
                  index=tdf_subset$type_no2
)

datalist <- list(X=X_sub,clinical=clinical)

m <- nimbleModel(nc,constants = constants,data=datalist)
inits <- function() list(mu_pp = runif(1,3.5,4.5), mu_pn = rnorm(1,4,0.1))

### 2.4. Run nimble code
mcmc.out <- nimbleMCMC(nc,constants = constants,data = datalist,
                       #inits=inits,
                       inits=list(sig_pp=10,sig_pn=10,mu_pp=100,mu_pn=4,pp=0.95,pn=0.8,logit_pp = logit(0.95),logit_pn=logit(0.8)),
                       monitors = c("theta","pp","pn","logit_pp","logit_pn","mu_pp","mu_pn","sig_pp","sig_pn"),
                       nburnin=10000,niter=20000,summary = TRUE,WAIC=TRUE,
                       nchains = 4,samplesAsCodaMCMC = TRUE,thin=40)
coln<-colnames(mcmc.out$samples[[1]])

### 2.5 Save hyperparmeter
mu_pp = mcmc.out$summary$all.chains[which(coln=="mu_pp"),1]
mu_pn = mcmc.out$summary$all.chains[which(coln=="mu_pn"),1]
sd_pp = mcmc.out$summary$all.chains[which(coln=="sig_pp"),1]
sd_pn = mcmc.out$summary$all.chains[which(coln=="sig_pn"),1]

save(mu_pp,mu_pn,sd_pp,sd_pn,file="hyperpar.Rdata")





