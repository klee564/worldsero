library(scales)
library(tidyverse)

### 6.1 Load data
load("final_result.Rdata")


### 6.2 Sample thetaV from alpha,beta
world_result<-full_result
for(i in 1:1000){
  thetaVV<-pmax(0,rbeta(nrow(world_result),world_result$alpha,world_result$beta),na.rm=TRUE)
  world_result<-cbind(world_result,thetaVV)
  world_result<-cbind(world_result,theta = 1-(1-world_result[,3+i])*(1-thetaVV))
}

ind_V<-which(colnames(world_result)=="thetaVV")
colnames(world_result)[ind_V]<-paste0("thetaV.",1:1000)

ind_total<-which(colnames(world_result)=="theta")
colnames(world_result)[ind_total]<-paste0("theta.",1:1000)

world_result<-world_result%>%dplyr::select(Country.Region,date,pop,starts_with("thetaI"),starts_with("thetaV"),starts_with("theta."))

world_result[4:ncol(world_result)]<-world_result[4:ncol(world_result)]*world_result$pop

world_result<-world_result%>%
  dplyr::select(date,starts_with("theta"))%>%
  group_by(date)%>%
  summarise(across(everything(),sum))%>%
  as.data.frame()

world_result[,-1]<-world_result[,-1]/sum(full_result%>%filter(date=="2021-01-01")%>%dplyr::select(pop))


### 6.3 Merge thetaI, thetaV, theta sample
trenddf<-matrix(NA,nrow=nrow(world_result)*3,ncol=5)
colnames(trenddf)<-c("date","low_cred","high_cred","mean","type")

trenddf<-as.data.frame(trenddf)

trenddf[,1]<-as.Date(rep(world_result$date,each=3))

for(i in 1:nrow(world_result)){
  theta_sample<-world_result%>%
    filter(date==world_result$date[i])%>%
    dplyr::select(starts_with("theta"))
  trenddf[3*i-2,2:4]<-c(HPDinterval(as.mcmc(theta_sample%>%dplyr::select(starts_with("thetaV"))%>%as.numeric()),prob=0.95),mean(theta_sample%>%dplyr::select(starts_with("thetaV"))%>%as.numeric()))
  trenddf[3*i-2,5]<-"theta^{(V)}"
  trenddf[3*i-1,2:4]<-c(HPDinterval(as.mcmc(theta_sample%>%dplyr::select(starts_with("thetaI"))%>%as.numeric()),prob=0.95),mean(theta_sample%>%dplyr::select(starts_with("thetaI"))%>%as.numeric()))
  trenddf[3*i-1,5]<-"theta^{(I)}"
  trenddf[3*i,2:4]<-c(HPDinterval(as.mcmc(theta_sample%>%dplyr::select(starts_with("theta."))%>%as.numeric()),prob=0.95),mean(theta_sample%>%dplyr::select(starts_with("theta."))%>%as.numeric()))
  trenddf[3*i,5]<-"theta"
}


### 6.4 Plotting trend of world seroprevalence
gtheta <- ggplot( trenddf%>%
                    dplyr::mutate(type = factor(type,levels = c("theta^{(V)}","theta^{(I)}","theta"))),aes((date),mean))+
  geom_line()+
  geom_ribbon(aes(ymin=low_cred,ymax=high_cred),alpha=0.3)+
  labs(x="Date",y=bquote(theta),title="Trend of world seroprevalence")+
  scale_x_date(labels= date_format("%b"))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"),
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
        legend.title=element_text(size=16,face="bold")) +
  ylab("seroprevalence")+
  facet_grid(cols = vars(type),labeller=label_parsed)

gtheta
