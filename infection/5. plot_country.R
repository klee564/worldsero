library(ggplot2)
library(lubridate)
library(patchwork)


### 5.1 Load data
load("final_result.Rdata")
load("Preprocessing_data.RData")
world_data1<-read.csv("world_data_v2.csv")

### 5.2 Timeplot function for thetaI
full_result$date<-as.Date(full_result$date)
full_result<-full_result %>%
  arrange(date)

timeplot_fun<-function(country){
  cc<-full_result%>%
    filter(Country.Region==country)%>%
    rowwise()%>%
    summarize(date=date,
              mean=mean(c_across(thetaI:thetaI.999)),
              low_cred=quantile(c_across(thetaI:thetaI.999),c(0.05)),
              high_cred=quantile(c_across(thetaI:thetaI.999),c(0.95)),
              thetaC=thetaC)
  
  p<-ggplot(cc,aes(date,mean))+
    geom_line()+
    geom_ribbon(aes(ymin=low_cred,ymax=high_cred),alpha=0.3)+
    geom_line(data=cc,aes(date,thetaC),color="orangered")+
    labs(x="",y="",title=country)+
    theme(plot.title = element_text(hjust = 0.5,size=16,face="bold"),
          #bold font for legend text
          legend.text=element_text(size=16),
          #set thickness of axis ticks
          axis.ticks=element_line(size=0.4),
          #remove plot background
          plot.background=element_blank(),
          #remove plot border
          panel.border=element_blank(),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=13.5),
          axis.title = element_text(size=15),
          legend.title=element_text(size=16,face="bold"))
  
  return(p)
}


### 5.2. Trace plot of thetaI for specific country
p1<-timeplot_fun("US")
p2<-timeplot_fun("Korea, South")
p3<-timeplot_fun("Japan")
p4<-timeplot_fun("China")
p1+p2+p3+p4+plot_layout(ncol=2)


### 5.3 Timeplot function for thetaV
for(i in 1:length(thetaV)){
  names(thetaV)[i]<-thetaV[[i]]$country[1]
}

timeplot_fun_V<-function(country){
  cc<-thetaV[[which(names(thetaV)==country)]]
  
  ccc<-cc%>%
    mutate(low_cred = qbeta(0.05,shape1,shape2),
           high_cred = qbeta(0.95,shape1,shape2))%>%
    dplyr::select(country,date,thetaVhat,low_cred,high_cred)
  
  if(min(ccc$date)>"2021-01-01"){
    ccc_vec<-ccc[1,]
    ccc_vec[2]<-"2021-01-01"
    ccc_vec[3:5]<-0
    ccc<-rbind(ccc_vec,ccc)
  }else if(!"2021-01-01"%in% ccc$date){
    closest_date<-which(abs(date(ccc$date) - date("2021-01-01")) == min(abs(date(ccc$date) - date("2021-01-01"))))[1]
    ccc_vec<-ccc[closest_date,]
    ccc_vec[2]<-"2021-01-01"
    ccc<-rbind(ccc_vec,ccc)
  }
  
  if(!"2021-07-31"%in% ccc$date){
    closest_date<-which(abs(date(ccc$date) - date("2021-07-31")) == min(abs(date(ccc$date) - date("2021-07-31"))))[1]
    ccc_vec<-ccc[closest_date,]
    ccc_vec[2]<-"2021-07-31"
    ccc<-rbind(ccc,ccc_vec)
  }
  
  ccc<-ccc%>%filter(date>="2021-01-01" & date<="2021-07-31")
  ccc$date<-as.Date(ccc$date)
  
  
  p<-ggplot(ccc,aes(date,thetaVhat,group=1))+
    geom_line()+
    geom_ribbon(aes(ymin=low_cred,ymax=high_cred),alpha=0.3)+
    #geom_line(data=cc,aes(date,thetaC),color="orangered")+
    labs(x="",y="",title=country)+
    theme(plot.title = element_text(hjust = 0.5,size=16,face="bold"),
          #bold font for legend text
          legend.text=element_text(size=16),
          #set thickness of axis ticks
          axis.ticks=element_line(size=0.4),
          #remove plot background
          plot.background=element_blank(),
          #remove plot border
          panel.border=element_blank(),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=13.5),
          axis.title = element_text(size=15),
          legend.title=element_text(size=16,face="bold"))
  
  return(p)
}

### 5.4 Trace plot of thetaI for specific country
p1<-timeplot_fun_V("US")
p2<-timeplot_fun_V("Korea, South")
p3<-timeplot_fun_V("Japan")
p4<-timeplot_fun_V("China")
p1+p2+p3+p4+plot_layout(ncol=2)


### 5.5 Timeplot function for whole theta
for(i in 1:length(thetaV)){
  names(thetaV)[i]<-thetaV[[i]]$country[1]
}

emp_vac_dt<-data.frame(date = as.Date("2021-01-01"):as.Date("2021-07-31")
                       ,shape1=NA,shape2 =NA)
emp_vac_dt[,1]<-as.Date(as_date(emp_vac_dt[,1]))


timeplot_fun_VI<-function(country){
  full_result_country<-full_result%>%
    filter(Country.Region==country)
  
  cc_vaccine<-thetaV[[which(names(thetaV)==country)]][,2:4]
  cc_vaccine$date<-as.Date(cc_vaccine$date)
  min_vac_date<-min(cc_vaccine$date)
  
  vac_dt<-left_join(emp_vac_dt,cc_vaccine,by="date")[,c(1,4,5)]
  vac_dt<-vac_dt%>%
    mutate(vac = ifelse(date>=min_vac_date,1,0))
  not_na_idx<-which(!is.na(vac_dt$shape1.y))
  for(i in 1:nrow(vac_dt)){
    if(vac_dt$vac[i]==1 & is.na(vac_dt$shape1.y[i])){
      closest_date<-which(abs(not_na_idx - i) == min(abs(not_na_idx - i)))[1]
      vac_dt[i,c(2,3)]<-vac_dt[not_na_idx[closest_date],c(2,3)]
    }
  }
  
  for(i in 1:nrow(full_result_country)){
    if(vac_dt$vac[i]==1){
      ran_num<-rbeta(length(4:(ncol(full_result_country)-3)),vac_dt$shape1.y[i],vac_dt$shape2.y[i])
      aaa<-as.numeric(full_result_country[i,4:(ncol(full_result_country)-3)])
      full_result_country[i,4:(ncol(full_result_country)-3)]<-
        aaa+ran_num-ran_num*aaa
    }
  }
  
  cc<-full_result_country%>%
    rowwise()%>%
    summarize(date=date,
              mean=mean(c_across(thetaI:thetaI.999)),
              low_cred=quantile(c_across(thetaI:thetaI.999),c(0.05)),
              high_cred=quantile(c_across(thetaI:thetaI.999),c(0.95)),
              thetaC=thetaC)
  thetaV_hat<-vac_dt$shape1.y/(vac_dt$shape1.y+vac_dt$shape2.y)
  thetaV_hat[which(is.nan(thetaV_hat))]<-0
  cc<-cbind(cc,thetaV= thetaV_hat)
  
  p<-ggplot(cc,aes(date,mean))+
    geom_line()+
    geom_ribbon(aes(ymin=low_cred,ymax=high_cred),alpha=0.3)+
    geom_line(data=cc,aes(date,thetaC),color="orangered")+
    geom_line(data=cc,aes(date,thetaV),color="blue")+
    labs(x="",y="",title=country)+
    theme(plot.title = element_text(hjust = 0.5,size=16,face="bold"),
          #bold font for legend text
          legend.text=element_text(size=16),
          #set thickness of axis ticks
          axis.ticks=element_line(size=0.4),
          #remove plot background
          plot.background=element_blank(),
          #remove plot border
          panel.border=element_blank(),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=13.5),
          axis.title = element_text(size=15),
          legend.title=element_text(size=16,face="bold"))
  
  return(p)
}

### 5.6 Trace plot of whole theta for specific country
# red : confirmed ratio
# blue : seroprevalence by vaccine
# black & grey : seroprevalence and credible interval 

p1<-timeplot_fun_VI("US")
p2<-timeplot_fun_VI("Korea, South")
p3<-timeplot_fun_VI("Japan")
p4<-timeplot_fun_VI("China")
p1+p2+p3+p4+plot_layout(ncol=2)

