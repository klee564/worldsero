library(tidyr)
library(dplyr)
library(truncnorm)
library(coda)

### 4.1 Load data
world_data1<-read.csv("world_data_v2.csv") 
load("Preprocessing_data.RData")
load("hyperpar.Rdata")
load("sample_posterior.Rdata")

### 4.2 Merge confirmed data and nation's statistics data
world_data1<-world_data1%>%
  filter(country!='Cuba')

merge_data<-merge(conf1,world_data1,by.x='Country.Region',by.y='country')

merge_data1<-na.omit(merge_data)

merge_data2<-merge_data1%>%
  gather('date','value',starts_with('202'))

merge_data2$value<-as.numeric(merge_data2$value)

country_list<-unique(tdf$Country.or.Area)

country_df<-tdf%>%
  dplyr::select(Country.or.Area,country_index)%>%
  distinct()


### 4.3 Sampling from predictive posterior distribution
md2<-merge_data2%>%
  filter(date>="2021-01-01"&date<="2021-07-31")

md2_sero<-md2%>%
  filter(Country.Region%in%country_list)
md2_sero<-merge(md2_sero,country_df,by.x="Country.Region",by.y="Country.or.Area")
md2_sero<-md2_sero%>%
  arrange(Country.Region,date)

md2_nonsero<-md2%>%
  filter(!Country.Region%in%country_list)
md2_nonsero<-md2_nonsero%>%
  arrange(Country.Region,date)

iter<-nrow(sample)

sample<-(sample[1:iter,1:ncol(sample)])

num_coun<-nrow(merge_data1)


for(i in 1:iter){
  md2_sero<-cbind(md2_sero,thetaI=NA)
  
  for(j in 1:length(country_list)){
    idx<-which(md2_sero$country_index==j)
    sample_tn<-rtruncnorm(1,a=0,b=-log(md2_sero$value[idx]/md2_sero$pop[idx]/1000),mean=(sample[i,which(coln1==paste0("beta[",j,"]"))]+
                                                                                           sample[i,which(coln1=="beta2")]*(log(md2_sero$gdp[idx])-scale_gdp_mean)/scale_gdp_sd+
                                                                                           sample[i,which(coln1=="beta1")]*(log(md2_sero$pop[idx]/md2_sero$area[idx])-scale_pd_mean)/scale_pd_sd)
                          ,
                          sd=(sample[i,which(coln1=="tau")]))
    
    md2_sero[idx,ncol(md2_sero)]<-exp(sample_tn+log(md2_sero$value[idx]/md2_sero$pop[idx]/1000))
    
    
  }
  
  sample_tn<-rtruncnorm(1,a=0,b=-log(md2_nonsero$value/md2_nonsero$pop/1000),mean=(rnorm(1,sample[i,which(coln1=="mu0")],sample[i,which(coln1=="sigma")])+
                                                                                     sample[i,which(coln1=="beta2")]*(log(md2_nonsero$gdp)-scale_gdp_mean)/scale_gdp_sd+
                                                                                     sample[i,which(coln1=="beta1")]*(log(md2_nonsero$pop/md2_nonsero$area)-scale_pd_mean)/scale_pd_sd)
                        ,
                        sd=(sample[i,which(coln1=="tau")]))
  
  md2_nonsero<-cbind(md2_nonsero,thetaI=exp(sample_tn+log(md2_nonsero$value/md2_nonsero$pop/1000)))                
  
  
}

md2_sero1<-md2_sero[,-c(9)]
colnames(md2_nonsero)<-colnames(md2_sero1)
md2<-rbind(md2_sero1,md2_nonsero)


md2$thetaC<-md2$value/md2$pop/1000


####rm(md2_nonsero,md2_sero,mcmc.out,mcmc.out1)
thetaV<-readRDS(file = "real_thetaV.rds")
thetaV<-thetaV$ThetaV



for(i in 1:length(thetaV)){
  names(thetaV)[i]<-thetaV[[i]]$country[1]
}


thetaI_mean<-rowMeans(md2[,9:(iter+8)])
final_result<-cbind(md2[,c(1,4,7,(iter+9))],thetaI=thetaI_mean)

final_result$alpha<-NA
final_result$beta<-NA

for(i in 1:length(thetaV)){
  names(thetaV)[i]<-thetaV[[i]]$country[1]
}

idx<-match(c("country","date","shape1","shape2"),colnames(thetaV[[1]]))
idx1<-match(c("alpha","beta"),colnames(final_result[1,]))

for(i in unique(final_result$Country.Region)){
  if(i%in%names(thetaV)){
    thetaV_i = thetaV[[which(names(thetaV)==i)]][,idx]
    names(thetaV_i)<-c("Country.Region","date","alpha","beta")
    thetaV_i$date<-as.Date(thetaV_i$date)
    if(min(thetaV_i$date)<=max(final_result[which(final_result$Country.Region==i),]$date)){
      date_idx<-which(final_result[which(final_result$Country.Region==i),]$date>=min(thetaV_i$date))
      for(j in date_idx){
        tdf_date<-as.Date(final_result[which(final_result$Country.Region==i),]$date[j])
        closest_date<-which(abs(thetaV_i$date - tdf_date) == min(abs(thetaV_i$date - tdf_date)))[1]
        final_result[which(final_result$Country.Region==i),][j,idx1]<-thetaV_i[closest_date,c(3,4)]
      }
    }
  }
}
final_result$thetaV<-final_result$alpha/(final_result$alpha+final_result$beta)



### 4.4 Credibls interval of world seroprevalence at 31th July, 2021
full_result<-md2[,c(1,4,7,9:(iter+9))]  
full_result$date<-as.Date(full_result$date)

full_result$alpha<-final_result$alpha
full_result$beta<-final_result$beta

result_0731<-full_result%>%filter(date=="2021-07-31")

for(i in 4:1003){
  result_0731[,i]<-1-(1-result_0731[,i])*pmin(1,(1-rbeta(nrow(result_0731),result_0731$alpha,result_0731$beta)),na.rm=TRUE)
}

world_0731<-rep(NA,1000)
for(i in 1:1000){
  world_0731[i]<-sum(result_0731[,i+3]*result_0731$pop/sum(result_0731$pop))
}
coda::HPDinterval(as.mcmc(world_0731),prob=0.95)


### 4.5 Credibls interval of ThetaI at 31th July, 2021
result_0731<-full_result%>%filter(date=="2021-07-31")

for(i in 4:1003){
  result_0731[,i]<-1-(1-result_0731[,i])
}

world_0731<-rep(NA,1000)
for(i in 1:1000){
  world_0731[i]<-sum(result_0731[,i+3]*result_0731$pop/sum(result_0731$pop))
}
coda::HPDinterval(as.mcmc(world_0731),prob=0.95)


### 4.6 Credibls interval of ThetaV at 31th July, 2021
result_0731<-full_result%>%filter(date=="2021-07-31")

for(i in 4:1003){
  result_0731[,i]<-pmax(0,rbeta(nrow(result_0731),result_0731$alpha,result_0731$beta),na.rm=TRUE)
}

world_0731<-rep(NA,1000)
for(i in 1:1000){
  world_0731[i]<-sum(result_0731[,i+3]*result_0731$pop/sum(result_0731$pop))
}
coda::HPDinterval(as.mcmc(world_0731),prob=0.95)


### 4.7 Save data
save(final_result,full_result,thetaV,file="final_result.Rdata")

