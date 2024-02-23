library(dplyr)
library(gtools)
library(lubridate)
library(readr)
library(tidyr)


### 1.1 Load data
getwd()
setwd("~/data/COVID/230425/")

serosurvey <- read.csv("serotracker_0731.csv")
clinicaltrial <- read.csv("clinicaltrial_final_0613.csv")
world_data1<-read.csv("world_data_v2.csv")

### 1.2 Select variable
filtered <- serosurvey %>%
  dplyr::select(Prevalence.Estimate.Name,
                URL,
                Serosurvey.sample.size,Test.positive.size,type_no,
                Last.collection.day,
                Country.or.Area,
                Denominator.Value,Serum.positive.prevalence)%>%
  dplyr::filter(Country.or.Area!="")

### 1.3 Chane name of nations
filtered$Country.or.Area[filtered$Country.or.Area=="Côte d'Ivoire"]<-"Cote d'Ivoire"
filtered$Country.or.Area[filtered$Country.or.Area=="Iran (Islamic Republic of)"]<-"Iran"
filtered$Country.or.Area[filtered$Country.or.Area=="Lao People's Democratic Republic"]<-"Laos"
filtered$Country.or.Area[filtered$Country.or.Area=="Republic of Korea"]<-"Korea, South"
filtered$Country.or.Area[filtered$Country.or.Area=="Russian Federation"]<-"Russia"
filtered$Country.or.Area[filtered$Country.or.Area=="The United Kingdom"]<-"United Kingdom"
filtered$Country.or.Area[filtered$Country.or.Area=="United States of America"]<-"US"
filtered<-filtered%>%
  filter(Country.or.Area!="French Guiana"&Country.or.Area!="Jersey")

### 1.4 Merge with nation statistics
filtered<-merge(filtered,world_data1,by.x='Country.or.Area',by.y='country') #GDP, POP, AREA

### 1.5 Change form of date
filtered$Last.collection.day<-as.Date(filtered$Last.collection.day,format = "%y.%m.%d.")
filtered$collection.day<-filtered$Last.collection.day

### 1.6 Calculate the number of positivity by using Serum positive prevalence
filtered$Test.positive.size[which(is.na(filtered$Test.positive.size))]<-
  filtered$Serosurvey.sample.size[which(is.na(filtered$Test.positive.size))]*
  as.numeric(substr(filtered$Serum.positive.prevalence[which(is.na(filtered$Test.positive.size))],1,nchar(filtered$Serum.positive.prevalence[which(is.na(filtered$Test.positive.size))])-1))/100

### 1.7 Load confirmed data
Confirmed = read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
Confirmed = Confirmed[, !colnames(Confirmed) %in% c("Lat", "Long")]
colnames(Confirmed)<-as.Date(colnames(Confirmed),format = "X%m.%d.%y")
colnames(Confirmed)[c(1,2)]<-c("Province.State","Country.Region")

conf1<-Confirmed[which(Confirmed$Province.State==""),]
country_list<-unique(Confirmed$Country.Region) # 193
setdiff(country_list,conf1$Country.Region)
conf1<-rbind(conf1,c("","Australia",as.numeric(colSums((Confirmed[9:16,3:ncol(conf1)])))))
conf1<-rbind(conf1,c("","Canada",as.numeric(colSums((Confirmed[40:55,3:ncol(conf1)])))))
conf1<-rbind(conf1,c("","China",as.numeric(colSums((Confirmed[59:92,3:ncol(conf1)])))))

### 1.8 Merge serosurvey and confirmed data
tdf<-filtered
tdf$confirmed<-0

for(i in 1:nrow(tdf)){
  coun<-tdf$Country.or.Area[i]
  date1<-tdf$collection.day[i]
  if(as.character(date1)%in%colnames(conf1)){
    tdf$confirmed[i]<-as.numeric(conf1[which(conf1$Country.Region==coun),which(colnames(conf1)==as.character(date1))])
  }
}



### 1.9 Calculate nation's population over 15 ages.
per_broad_age <- read_excel("per_broad_age_v2.xlsx", 
                            sheet = "Data")
per_broad_age<-per_broad_age[-(1:4),c(2,4)]
colnames(per_broad_age)<-c("country","percent")


per_broad_age[which(per_broad_age$country=="Bolivia (Plurinational State of)"),1]<-"Bolivia"
per_broad_age[which(per_broad_age$country=="Dem. Rep. of the Congo"),1]<-"Congo (Kinshasa)"
per_broad_age[which(per_broad_age$country=="Congo"),1]<-"Congo (Brazzaville)"
per_broad_age[which(per_broad_age$country=="Brunei Darussalam"),1]<-"Brunei"
per_broad_age[which(per_broad_age$country=="State of Palestine"),1]<-"West Bank and Gaza"
per_broad_age[which(per_broad_age$country=="Côte d'Ivoire"),1]<-"Cote d'Ivoire"
per_broad_age[which(per_broad_age$country=="Kosovo (under UNSC res. 1244)"),1]<-"Kosovo"
per_broad_age[which(per_broad_age$country=="Lao People's Dem. Republic"),1]<-"Laos"
per_broad_age[which(per_broad_age$country=="Republic of Moldova"),1]<-"Moldova"
per_broad_age[which(per_broad_age$country=="Micronesia (Fed. States of)"),1]<-"Micronesia"
per_broad_age[which(per_broad_age$country=="Syrian Arab Republic"),1]<-"Syria"
per_broad_age[which(per_broad_age$country=="Türkiye"),1]<-"Turkey"
per_broad_age[which(per_broad_age$country=="United Republic of Tanzania"),1]<-"Tanzania"
per_broad_age[which(per_broad_age$country=="Venezuela (Bolivarian Republic of)"),1]<-"Venezuela"
per_broad_age[which(per_broad_age$country=="Viet Nam"),1]<-"Vietnam"
per_broad_age[which(per_broad_age$country=="Myanmar"),1]<-"Burma"
per_broad_age[which(per_broad_age$country=="Iran (Islamic Republic of)"),1]<-"Iran"
per_broad_age[which(per_broad_age$country=="Republic of Korea"),1]<-"Korea, South"
per_broad_age[which(per_broad_age$country=="Russian Federation"),1]<-"Russia"
per_broad_age[which(per_broad_age$country=="United States of America"),1]<-"US"


### 1.10 Add new variables (confirmed_ratio,GDP per capita, adult population density) 
tdf$country_index <-
  as.integer(as.factor(tdf$Country.or.Area))

tdf$pop <- tdf$pop*as.numeric(per_broad_age$percent[match(tdf$Country.or.Area,per_broad_age$country)])/100

tdf$thetac <-
  (tdf$confirmed/
     tdf$pop/1000)


### 1.11 Exclude the serosurvey if the seroprevalence is lower than the COVID confirmed population rates
idx<-which(tdf$thetac<parse_number(tdf$Serum.positive.prevalence)/100)
tdf<-tdf[idx,]


### 1.12 Scaling the variables
tdf$popul_density<- tdf$pop/tdf$area

scale_pd_mean<-mean(log(tdf$popul_density))
scale_pd_sd<-sd(log(tdf$popul_density))
tdf$log_popul_density_sc<- scale(log(tdf$popul_density))

scale_gdp_mean<-mean(log(tdf$gdp))
scale_gdp_sd<-sd(log(tdf$gdp))
tdf$log_gdp_sc<- scale(log(tdf$gdp))

tdf$Serosurvey.sample.size<-as.numeric(tdf$Serosurvey.sample.size)


### 1.13 Add vaccination information
thetaV<-readRDS(file = "real_thetaV.rds")
thetaV<-thetaV$ThetaV

tdf$alpha<-NA
tdf$beta<-NA

for(i in 1:length(thetaV)){
  names(thetaV)[i]<-thetaV[[i]]$country[1]
}

idx<-match(c("country","date","shape1","shape2"),colnames(thetaV[[1]]))
idx1<-match(c("alpha","beta"),colnames(tdf[1,]))

for(i in unique(tdf$Country.or.Area)){
  if(i%in%names(thetaV)){
    
    thetaV_i = thetaV[[which(names(thetaV)==i)]][,idx]
    names(thetaV_i)<-c("country","date","alpha","beta")
    thetaV_i$date<-as.Date(thetaV_i$date)
    if(min(thetaV_i$date)<=max(tdf[which(tdf$Country.or.Area==i),]$collection.day)){
      date_idx<-which(tdf[which(tdf$Country.or.Area==i),]$collection.day>=min(thetaV_i$date))
      for(j in date_idx){
        tdf_date<-tdf[which(tdf$Country.or.Area==i),]$collection.day[j]
        closest_date<-which(abs(thetaV_i$date - tdf_date) == min(abs(thetaV_i$date - tdf_date)))[1]
        tdf[which(tdf$Country.or.Area==i),][j,idx1]<-thetaV_i[closest_date,c(3,4)]
      }
    }
  }
}

tdf<-tdf%>%arrange(alpha)

### 1.14 Reindexing
tdf$country_index <- as.integer(as.factor(tdf$Country.or.Area))
tdf$type_no2 <- as.integer(as.factor(tdf$type_no))

### 1.15. Serotest data
clinicaltrial<-clinicaltrial%>%
  filter(!is.na(Type_no))

tclinical <- clinicaltrial %>% filter(Type_no %in% tdf$type_no) %>%
  dplyr::select(Type_no,TP,FN,TN,FP) %>% arrange(Type_no)

tclinical$Type_no2<-as.integer(as.factor(tclinical$Type_no))

### 1.16. Save data

save(tdf,tclinical,conf1,country_list,scale_pd_mean,scale_pd_sd,scale_gdp_mean,scale_gdp_sd,file="Preprocessing_data.RData")




