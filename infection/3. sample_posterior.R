library(nimble)
library(purrr)
library(ggplot2)

### 3.1 Load data
load("Preprocessing_data.RData")
load("hyperpar.Rdata")

### 3.2 Nimble code
nc1 <- nimbleCode({
  
  
  for(i in 1:V){
    theta[i]<-thetaI[i]+thetaV[i]-thetaI[i]*thetaV[i];
    thetaV[i] ~ dbeta(alphaV[i],betaV[i]);
  }
  for(i in (V+1):N){
    theta[i]<-thetaI[i];
  }
  
  for(i in 1:N){
    X[i] ~ dbin(pp[index[i]]*(theta[i])+
                  (1-pn[index[i]])*(1-(theta[i])),S[i]);
    
  }
  
  for(i in 1:M){
    pp[i] ~ dbeta(clinical[i,1]+0.5,clinical[i,2]+0.5);
    pn[i] ~ dbeta(clinical[i,3]+0.5,clinical[i,4]+0.5);
    
  }
  
  for(i in (M+1):MM){
    logit(pp[i]) ~ T(dnorm(mu_pp,sd=sd_pp),0,Inf);
    logit(pn[i]) ~ T(dnorm(mu_pn,sd=sd_pn),0,Inf);
    
  }
  
  tau ~ dflat();
  beta1 ~ dflat();
  beta2 ~ dflat();
  
  for(i in 1:N){
    diff_theta[i] ~ dnorm(beta[country[i]]+beta1*log_pd[i]+beta2*log_gdp[i], sd=tau);
    thetaI[i] <- exp((diff_theta[i])+log(thetaC[i]));
    constraint_theta[i] ~ dconstraint(thetaI[i]-thetaC[i]>=0 & thetaI[i]<=1);
  }
  
  sigma ~ dflat();
  mu0 ~ dflat();
  for(i in 1:K){
    beta[i] ~ dnorm(mu0,sd=sigma);
  }
  
})


### 3.3 Parameter of nimble codes
K <- max(tdf$country_index)
M <- dim(tclinical)[1]
S <- tdf$Serosurvey.sample.size
X <- as.integer(tdf$Test.positive.size)
V <- length(which(tdf$alpha>0))
clinical<- tclinical[,-c(1,6)]

tdf[which(is.na(tdf$type_no2)),]$type_no2<-(M+1):(M+length(which(is.na(tdf$type_no2))))

constants1 <- list(N=dim(tdf)[1],S=S,
                   M=M,K=K,V=V,MM = max(tdf$type_no2),
                   index=tdf$type_no2,country=tdf$country_index,
                   country_num = table(tdf$country_index),
                   mu_pp = mu_pp,
                   mu_pn = mu_pn,
                   sd_pp = sd_pp,
                   sd_pn = sd_pn
)

datalist1 <- list(X=X,clinical=clinical,
                  log_pd = as.vector(tdf$log_popul_density_sc),
                  thetaC = tdf$thetac,
                  log_gdp = as.vector(tdf$log_gdp_sc),
                  alphaV = tdf$alpha[1:V],betaV = tdf$beta[1:V],
                  constraint_theta = rep(1,dim(tdf)[1])
)

### 3.4 Run nimble codes
m1 <- nimbleModel(nc1,constants = constants1,data=datalist1)

samplepp <- tclinical$TP/(tclinical$TP+tclinical$FN) - 0.01
for(i in (M+1):max(tdf$type_no2)){
  samplepp<-c(samplepp,mean(tclinical$TP/(tclinical$TP+tclinical$FN))-0.01)
}
samplepn <- tclinical$TN/(tclinical$TN+tclinical$FP) - 0.01
for(i in (M+1):max(tdf$type_no2)){
  samplepn<-c(samplepn,mean(tclinical$TN/(tclinical$TN+tclinical$FP))-0.01)
}

mcmc.out1 <- nimbleMCMC(nc1,constants = constants1,data = datalist1,
                        inits = list(beta=0,beta1=0,beta2=1,pp=samplepp,pn=samplepn,diff_theta=0.5,tau=1,mu0=0,sigma=1),
                        monitors = c("theta","thetaI","thetaV","beta","beta1","beta2","diff_theta","pp","pn","tau","logit_pp","logit_pn","mu0","sigma"),
                        nburnin=100000,niter=200000,summary = TRUE,WAIC=TRUE,
                        nchains = 4,samplesAsCodaMCMC = TRUE,thin=400)
coln1<-colnames(mcmc.out1$samples[[1]])

### 3.4 Check convergence of mcmc sample
df <-1:4 %>% 
  purrr::map(~data.frame(chain=.x,iter=1:dim(mcmc.out1$samples[[.x]])[1],
                         beta1=as.vector(mcmc.out1$samples[[.x]][,which(coln1=="beta1")]),
                         beta2=as.vector(mcmc.out1$samples[[.x]][,which(coln1=="beta2")]),
                         tau=as.vector(mcmc.out1$samples[[.x]][,which(coln1=="tau")]),
                         mu0=as.vector(mcmc.out1$samples[[.x]][,which(coln1=="mu0")]),
                         sigma=as.vector(mcmc.out1$samples[[.x]][,which(coln1=="sigma")])
                         )) %>% 
  do.call("rbind",.)


df$chain <- factor(df$chain)


ggplot(df,aes(x=iter,y=beta1)) + geom_line(aes(color=chain))
ggplot(df,aes(x=iter,y=beta2)) + geom_line(aes(color=chain))
ggplot(df,aes(x=iter,y=tau)) + geom_line(aes(color=chain))
ggplot(df,aes(x=iter,y=mu0)) + geom_line(aes(color=chain))
ggplot(df,aes(x=iter,y=sigma)) + geom_line(aes(color=chain))


acf(df$beta1[])
acf(df$beta2[])
acf(df$mu0[])
acf(df$sigma[])
acf(df$tau[])



### 3.5 Save sample of posterior
sample<-rbind(mcmc.out1$samples$chain1,mcmc.out1$samples$chain2,mcmc.out1$samples$chain3,mcmc.out1$samples$chain4)


save(sample,coln1,file="sample_posterior.Rdata")


