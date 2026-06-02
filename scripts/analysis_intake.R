# Code to run GAMS on milk intake rates
library(readxl)
library(ggplot2)
library(tidyverse)
library(mgcv)
library(gratia)
library(here)
library(broom)


# Load data ---------------------------------------------------------------

milkintake<-readxl::read_xlsx(here("data","Milk intake data.xlsx"), sheet="Milk intake data - final") |>
  filter(is.na(Exclude) & is.na(Perinatal))|>
  mutate(Family=as.factor(Family), CommonName=as.factor(CommonName))|>
  setNames(c("CommonName","SciName","Family","Order","N","LitterSize","MaternalMass","OffspringMass","Age","LactationDuration","TimeIntoLactation",
                        "MilkIntakegday","MilkIntakeSE","MilkIntakeSD","MilkIntakeMin", "MilkIntakeMax"
                        ,"Method","Reference","Comments","Exclude","Checked","OffspringMassSD","BirthMass","FTDuration","Perinatal"))

#milkintake<-data.table::fread(here("data-sub", "Milk intake data.csv)
# Summarise data ----------------------------------------------------------

milkintakeSum<-milkintake |>
  filter((is.na(Comments) | Comments!="On a per litter basis") & TimeIntoLactation>=0.03) |>
  group_by(Order)|>
  dplyr::mutate(nOrder=length(unique(CommonName))) |>
  group_by(Order, Family)|>
  dplyr::mutate(nFamily=length(unique(CommonName)))|>
  group_by(Order,Family,CommonName)|>
  dplyr::summarise(nInd=sum(N),nTimes=length(unique(Age)), nFamily=nFamily[1], nOrder=nOrder[1])

milkintakeSpecies<-milkintake |>
  filter((is.na(Comments) | Comments!="On a per litter basis") & TimeIntoLactation>=0.03) |>
  group_by(CommonName)|>
  dplyr::summarise(across(where(is.numeric), .f=list(mean=mean, sd=sd)))

milkintakeSpeciesSum<-milkintake |>
  filter((is.na(Comments) | Comments!="On a per litter basis") & TimeIntoLactation>=0.03) |>
  group_by(CommonName)|>
  dplyr::summarise(nInd=sum(N),LactCoverage=pmin(1,pmin(1,max(TimeIntoLactation))-min(TimeIntoLactation))/0.97, 
                   LactDuration=mean(LactationDuration), mBirthMass=mean(BirthMass, na.rm=T),
                   mMatMass=mean(MaternalMass, na.rm=T),mMass=mean(OffspringMass))

milkintakeSpeciesSum2<-milkintake |>
  filter((is.na(Comments) | Comments!="On a per litter basis") & TimeIntoLactation>=0.03) |>
  group_by(Family,CommonName)|>
  dplyr::summarise(LactCoverage=pmin(1,pmin(1,max(TimeIntoLactation))-min(TimeIntoLactation))/0.97,
                   N=sum(N), LactDur=mean(LactationDuration)) |>
  ungroup()

milkintakeSpeciesSum3<-milkintakeSpeciesSum2 |>
  group_by(Family)|>
  dplyr::summarise(MLactCoverage=mean(LactCoverage,na.rm=T), MN=mean(N, na.rm=T), MinN=min(N,na.rm=T), MaxN=max(N,na.rm=T), sdN=sd(N,na.rm=T))


# Simulate data -----------------------------------------------------------

newdata<-list()
set.seed(553)
for (i in 1:50){
  newdata[[i]]<-milkintake |>
    filter(N>1)|>
    rowwise()|>
    mutate(SMilkIntakegday=list(rnorm(N, mean=MilkIntakegday, sd=ifelse(!is.na(MilkIntakeSD),MilkIntakeSD,0))),
           SOffspringMass=list(rnorm(N, mean=OffspringMass, sd=ifelse(!is.na(OffspringMassSD), OffspringMassSD,0)))) |>
    unnest(cols=c(SMilkIntakegday,SOffspringMass)) |>
    plyr::rbind.fill(filter(milkintake, N==1))|>
    mutate(Population=i)
  print(i)
}

milkintakeNew<-bind_rows(newdata) |>
  mutate(SMilkIntakegday=ifelse(N==1,MilkIntakegday,SMilkIntakegday),
         SOffspringMass=ifelse(N==1,OffspringMass, SOffspringMass),
         response=SMilkIntakegday/(SOffspringMass^0.82)) |>
  filter(TimeIntoLactation>=0.03)

# Check to see how well simulation does
simCheck<-subset(milkintakeNew, N>1) |>
  group_by(CommonName, Population,Reference,Age)|>
  summarise(n=mean(N),RealMilk=mean(MilkIntakegday), RealMass=mean(OffspringMass),
            Real=mean(MilkIntakegday/OffspringMass^0.82), SimMilk=mean(SMilkIntakegday), 
            SimMass=mean(SOffspringMass),Sim=mean(SMilkIntakegday/SOffspringMass^0.82),
            SimDiff=(Real-Sim)/Real)

# Analysis ----------------------------------------------------------------

HFamilyOtariid<-milkintakeNew |>
  filter(Family=="Otariidae")|>
  group_by(Order,Family,Population)|>
  nest()|>
  mutate(modelG = map(data, ~ mgcv::gam(response ~ s(CommonName,bs="re") + s(TimeIntoLactation,bs="tp",k=5),
                                        ,data = .x, method="REML", family=gaussian(link="log"))))

HFamilyUrsidCub<-milkintakeNew |>
  filter(Family=="Ursidae" & Comments=="On a per cub basis")|>
  group_by(Order,Family,Population)|>
  nest()|>
  mutate(modelG = map(data, ~ mgcv::gam(response ~ s(CommonName,bs="re") + s(TimeIntoLactation, bs="tp",k=5),
                                        data = .x, method="REML", family=gaussian(link="log"))))

HFamilyMustelidKit<-milkintakeNew |>
  filter(Family=="Mustelidae" & Comments=="On a per kit basis")|>
  group_by(Order,Family,Population)|>
  nest()|>
  mutate(modelG = map(data, ~ mgcv::gam(response ~ s(TimeIntoLactation, bs="tp",k=5),
                                        data = .x, method="REML", family=gaussian(link="log"))))

HFamilyMustelidLitter<-milkintakeNew |>
  filter(Family=="Mustelidae" & Comments=="On a per litter basis")|>
  group_by(Order,Family,Population)|>
  nest()|>
  mutate(modelG = map(data, ~ mgcv::gam(response ~ s(TimeIntoLactation, bs="tp",k=5),
                                        data = .x, method="REML", family=gaussian(link="log"))))

# Just see how they compare
draw(compare_smooths(HFamilyMustelidLitter$modelG[[1]],HFamilyMustelidKit$modelG[[1]]))

HFamilyPhocid<-milkintakeNew |>
  filter(Family=="Phocidae")|>
  group_by(Order,Family,Population)|>
  nest()|>
  mutate(modelG = map(data, ~ mgcv::gam(response ~ s(CommonName, bs="re") + s(TimeIntoLactation, bs="tp",k=5),
                                        ,data = .x, method="REML", family=gaussian(link="log"))))
        
HFamilyCervidae<-milkintakeNew |>
  filter(Family=="Cervidae")|>
  group_by(Order,Family,Population)|>
  nest()|>
  mutate(modelG = map(data, ~ mgcv::gam(response~  s(TimeIntoLactation, bs="tp",k=5) + s(CommonName, bs="re"),
                                        data = .x, method="REML", family=gaussian(link="log"))))

HFamilyBovidae<-milkintakeNew |>
  filter(Family=="Bovidae")|>
  group_by(Order,Family,Population)|>
  nest()|>
  mutate(modelG = map(data, ~ mgcv::gam(response~  s(TimeIntoLactation,bs="tp",k=5),
                                        data = .x, method="REML", family=gaussian(link="log"))))

HFamilyCamelidae<-milkintakeNew |>
  filter(Family=="Camelidae")|>
  group_by(Order,Family,Population)|>
  nest()|>
  mutate(modelG = map(data, ~ mgcv::gam(response~  s(TimeIntoLactation,bs="tp",k=5),
                                        data = .x, method="REML", family=gaussian(link="log"))))

HFamilyArtiodactyla<-milkintakeNew |>
  filter(Order=="Artiodactyla" & paste(Family,Population,sep="_")!="Camelidae_33")|>
  group_by(Population)|>
  nest()|>
  mutate(Order="Artiodactyla",Family="Artiodactyla",
            modelG = map(data, ~ mgcv::gam(response~  s(TimeIntoLactation, bs="tp") + s(CommonName, bs="re",k=5),
                                        data = .x, method="REML", family=gaussian(link="log"))))

draw(compare_smooths(HFamilyArtiodactyla$modelG[[1]],HFamilyBovidae$modelG[[1]],HFamilyCervidae$modelG[[1]],HFamilyCamelidae$modelG[[1]]))

# Save model outputs -------------------------------------------------

HFamily<-bind_rows(HFamilyPhocid, HFamilyOtariid, HFamilyArtiodactyla,
                   HFamilyUrsidCub,HFamilyMustelidKit)

saveRDS(HFamily, file=(here("output","Milk intake GAM output.rds")))


# Predictions -------------------------------------------------------------

linkTransform<-gratia::inv_link(HFamily$modelG[[1]])

PFamily<-HFamily |>
  mutate(predicted = map(modelG,~ tidygam::predict_gam(.x, values=list(TimeIntoLactation=seq(0,1,by=0.001)),tran_fun=linkTransform))) |>
  unnest(predicted) |>
  select(-c(modelG,data)) |>
  mutate(across(c(response, se, lower_ci, upper_ci), round,digits=3))

fst::write_fst(PFamily,here("predictions","Milk intake GAM predictions.fst"), compress = 100)

