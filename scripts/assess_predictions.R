# Code to compare independent estimates with predictions

library(tidyverse)
library(here)
library(collapse)
library(rphylopic)

pluck_multiple <- function(x, ...) {
  `[`(x, ...)
}

familycolors=pals::kelly(22)

# Read in predictions -----------------------------------------------------

PFamilyE<-fst::read_fst(here("predictions","Milk energy density GAM predictions.fst")) |>
  mutate(TimeIntoLactation=factor(TimeIntoLactation))

PFamilyI<-fst::read_fst(here("predictions","Milk intake GAM predictions.fst"))|>
  mutate(TimeIntoLactation=factor(TimeIntoLactation),
         CommonName=case_when(Family=="Mustelidae"~"American mink",.default=CommonName)) 

Family<-PFamilyI |>
  select(Order,  Family, CommonName)|>
  distinct()

PFamilyE2<-PFamilyE
colnames(PFamilyE2)[c(1:2,4)]<-paste0(colnames(PFamilyE)[c(1:2,4)],2)

# Create summaries -----------------------------------------------------

milkIntakeSum<-PFamilyI |>
  group_by(Order,Family,CommonName,TimeIntoLactation)|>
  summarise(MIntake=mean(response))|>
  ungroup(TimeIntoLactation)|> 
  summarise_at(vars(MIntake),list(mean=mean,max=max), na.rm = TRUE) |>
  data.table::setnames(c("mean","max"),c("MIntake","MaxIntake")) |>
  arrange(MIntake)

milkEnergySum<-PFamilyE|>
  group_by(Order,Family,CommonName,TimeIntoLactation)|> # Average across populations
  summarise(MEnergy=mean(response)) |>
  ungroup(TimeIntoLactation)|> # Average across time
  summarise_at(vars(MEnergy),list(mean=mean,max=max), na.rm = TRUE) |>
  data.table::setnames(c("mean","max"),c("MEnergy","MaxEnergy")) |>
  arrange(MEnergy)


# Spotted seal --------------------------------------------------

# Data from Zhang et al 2014
phlaWean<-18

# Join with individual species data
pl<-data.frame(Age=1:10, OffspringMass=c(8.5,8.5,9.1,9.8,10.6,11.8,13.1,14.3,15.7,17)*1000) |>
  mutate(TimeIntoLactation=factor(round(Age/phlaWean,digits=3))) |>
  left_join(PFamilyI,relationship = "many-to-many")|>
  mutate(MilkIntake=response*(OffspringMass^0.82), MilkIntakeL=lower_ci*(OffspringMass^0.82), MilkIntakeU=upper_ci*(OffspringMass^0.82))|>
  arrange(Population, Age) |>
  select(-c(response, se, lower_ci, upper_ci))|>
  collapse::join(PFamilyE2,how="full",multiple=T)|>
  mutate(TimeIntoLactation=as.numeric(paste(TimeIntoLactation)), MilkEnergyD=response, 
         MilkEnergy=MilkIntake*MilkEnergyD, MilkEnergyL=MilkIntakeL*lower_ci, MilkEnergyU=MilkIntakeU*upper_ci,Reference="Prediction", Type="Model",
         CommonNameInt=paste(CommonName, CommonName2)) |>
  select(-c(response, se, lower_ci, upper_ci)) |>
  filter(!is.na(MilkEnergy) & (as.character(CommonName)==as.character(CommonName2) | CommonName2=="Spotted seal")) 

plSum<- pl |>
  filter(Family=="Phocidae" & Family2=="Phocidae")|>
  group_by(Order,Family,CommonName,Family2,CommonName2,Population,Reference)|>
  dplyr::summarise(MilkEnergyT=sum(MilkEnergy))|>
  ungroup()|>
  group_by(Order, Family, CommonName, Family2,CommonName2,Reference)|>
  dplyr::summarise(MilkEnergy=mean(MilkEnergyT))|>
  ungroup()|>
  mutate(CompareZhangMJ=(MilkEnergy-395200)/1000, CompareAbs=abs(CompareZhangMJ),
         CompareAvePerc=(CompareZhangMJ/395.200)) |>
  arrange(CompareAbs) 

# Comparison of daily estimates
pldata<-data.frame(Age=1:10, ActualMilkIntake=c(22.2,28.4,36.1,38.8,41.4,42,43.6,45.2,47.6,49.9)*1000)

pldaily<-pl |>
  filter(Family=="Phocidae" & Family2=="Phocidae")|>
  left_join(pldata) |>
  mutate(DailyDiff=MilkEnergy-ActualMilkIntake, DailyDiffProp=DailyDiff/ActualMilkIntake) |>
  group_by(CommonName, CommonName2, Age)|>
  dplyr::summarise(MDiff=mean(DailyDiff)/1000, MPropDiff=mean(DailyDiffProp))  |>
  ungroup()

# Calculate overlap 
ploverlap<-pl |>
  filter(Family=="Phocidae" & Family2=="Phocidae")|>
  group_by(CommonName, CommonName2,CommonNameInt,Age,TimeIntoLactation)|>
  summarise_at(c("MilkEnergy"),.funs="mean") |>
  ungroup()|>
  nest_by(CommonName,CommonName2,CommonNameInt,.key="prediction") |>
  mutate(temp=list(list(pluck(prediction$MilkEnergy) |>as.vector(),pldata$ActualMilkIntake)),
         overlap=flatten(list(overlapping::overlap(temp,type="1")))) |>
  unnest(overlap) |>
  select(-c(prediction, temp))|>
  arrange(overlap)

# These two species bracket the independent estimate when combined with spotted seal milk data
subset(plSum, CommonName=="Bearded seal" | CommonName=="Harp seal")

# Southern elephant seal --------------------------------------------------
sesM<-data.frame(Age=0:22, growthRate=3.8, LactationDuration=23,
                 BirthMass=45.5)
sesM$Mass<-(sesM$BirthMass+sesM$Age*sesM$growthRate)*1000
sesM$Sex<-"Male"

sesF<-data.frame(Age=0:22, growthRate=3.6, LactationDuration=23, BirthMass=40.5)
sesF$Mass<-(sesF$BirthMass+sesF$Age*sesF$growthRate)*1000
sesF$Sex<-"Female"

# From Fedak et al 1996 - Tables 2 and 3
sesdata1<-data.frame(EnergyExpend=1965, EnergyLoss=4414) |>
  mutate(MilkEnergy=(EnergyLoss-EnergyExpend)*0.95*1000, Reference="Fedak et al. 1996", CommonName="Southern elephant seal")

# From Arnbom et al 1997. Non-lactating female mass loss was truncated using relative changes for lactating females
#23.1 MJ/kg is the mean of values from Fedak et al 1996
sesdata2<-data.frame(NLFDMassLoss=3.4, LFDMassLoss=7.9)|>
  mutate(MilkEnergy=(LFDMassLoss-NLFDMassLoss)*23*23.1*1000*0.95,Reference="Arnbom et al. 1997", CommonName="Southern elephant seal")

# Combined with milk intake rates
ses<-rbind(sesM, sesF) |>
  group_by(Age) |>
  dplyr::summarise(OffspringMass=mean(Mass), LactationDuration=mean(LactationDuration)) |>
  ungroup()|>
  mutate(Species="Southern elephant seal",TimeIntoLactation=factor(round(Age/LactationDuration, digits=3)))|>
  left_join(PFamilyI,relationship = "many-to-many")|>
  mutate(MilkIntake=response*(OffspringMass^0.82))|>
  arrange(Population, Age) |>
  select(-c(response, se, lower_ci, upper_ci))|>
  collapse::join(PFamilyE2,how="full",multiple=T)|>
  mutate(TimeIntoLactation=as.numeric(paste(TimeIntoLactation)), MilkEnergyD=response, 
         MilkEnergy=MilkIntake*MilkEnergyD, Reference="Prediction", Type="Model") |>
  select(-c(response, se, lower_ci, upper_ci)) |>
  filter(!is.na(MilkEnergy) & (as.character(CommonName)==as.character(CommonName2) | CommonName2=="Southern elephant seal")) 

sesSum<-ses |>
  group_by(Species,Order,Family,CommonName,Family2,CommonName2,Population,Reference)|>
  dplyr::summarise(MilkEnergyT=sum(MilkEnergy))|>
  ungroup()|>
  group_by(Species,Order, Family, CommonName, Family2,CommonName2,Reference)|>
  dplyr::summarise(MilkEnergy=mean(MilkEnergyT))|>
  ungroup() |>
  mutate(CompareFedakMJ=(MilkEnergy-sesdata1$MilkEnergy)/1000, CompareArnbomMJ=(MilkEnergy-sesdata2$MilkEnergy)/1000,
         CompareAvg=(MilkEnergy-((sesdata1$MilkEnergy+sesdata2$MilkEnergy)/2))/1000,
         CompareAvgAbs=abs(CompareAvg), CompareAvePerc=CompareAvg/(((sesdata1$MilkEnergy+sesdata2$MilkEnergy)/2)/1000)) |>
  arrange(CompareAvgAbs) |>
  mutate(Substrate="Species", MilkEnergyType=case_when(CommonName2=="Southern elephant seal"~"Southern elephant seal",.default="Species"))  |>
  filter(Family=="Phocidae" & Family2=="Phocidae")

milkIntakeSum |> filter(CommonName=="Crabeater seal" | CommonName=="Gray seal") 

# Sea otter ---------------------------------------------------------------
# Nicole provided data
eldata<-readxl::read_xlsx(here("data","Female & Pup Figure 4 Data 5.1.14.xlsx"), sheet="Sheet2") |>
  janitor::clean_names()|>
  mutate(ActualMilkEnergy=mj_day_2*1.19*1000, Age=days_post_partum)|>
  filter(days_post_partum>=0 & days_post_partum<28)

elWean<-6/12*365
elAge<-0:elWean
elMilkSpecies<-milkEnergySum |>
  filter(MaxEnergy>=11 & MaxEnergy<13)

el<-data.frame(Mass=eldata$pup_mass, Age=eldata$Age)|>
  mutate(OffspringMass=Mass*1000, TimeIntoLactation=as.factor(round(Age/elWean, digits=3))) |>
  left_join(PFamilyI)|>
  mutate(MilkIntake=response*(OffspringMass^0.82), MilkIntakeL=lower_ci*(OffspringMass^0.82), MilkIntakeU=upper_ci*(OffspringMass^0.82))|>
  arrange(Population, Age) |>
  select(-c(response, se, lower_ci, upper_ci))|>
  collapse::join(PFamilyE2,how="full",multiple=T)|>
  collapse::fmutate(TimeIntoLactation=as.numeric(paste(TimeIntoLactation)), MilkEnergyD=response, 
                    MilkEnergy=MilkIntake*MilkEnergyD,MilkEnergyL=MilkIntakeL*lower_ci, MilkEnergyU=MilkIntakeU*upper_ci,
                    Reference="Prediction",Type="Model",
                    CommonNameInt=paste(CommonName, CommonName2)) |>
  collapse::fselect(-c(response, se, lower_ci, upper_ci))|>
  filter(!is.na(MilkEnergy))|>
  droplevels() |>
  filter(CommonName==CommonName2 | CommonName2 %in% elMilkSpecies$CommonName)

# Daily differences
elDayDiff<-el |>
  group_by(CommonName, CommonName2,Age,TimeIntoLactation)|>
  summarise_at(c("MilkEnergy","MilkEnergyL","MilkEnergyU"),.funs="mean") |>
  ungroup()|>
  left_join(eldata) |>
  mutate(DailyDiff=MilkEnergy-ActualMilkEnergy, DailyDiffProp=DailyDiff/ActualMilkEnergy) |>
  filter(Age<=27) |>
  group_by(CommonName, CommonName2)|>
  dplyr::summarise(MMilkEnergy=mean(MilkEnergy),MDiff=mean(abs(DailyDiff))/1000, MPropDiff=mean(abs(DailyDiffProp)), MinPropDiff=min(DailyDiffProp), MaxPropDiff=max(DailyDiffProp))  |>
  arrange(abs(MPropDiff))

elCompare<-eldata |>
  filter(Age<=27)|>
  summarise_at(c("ActualMilkEnergy"),"sum")

# Total differences
elTotalDiff<-  el |>
  filter(Age<=27)|>
  group_by(CommonName, CommonName2,Population)|>
  summarise_at(c("MilkEnergy","MilkEnergyL","MilkEnergyU"),.funs="sum") |>
  ungroup() |>
  group_by(CommonName, CommonName2)|>
  summarise_at(c("MilkEnergy","MilkEnergyL","MilkEnergyU"),.funs="mean") |>
  group_by(CommonName, CommonName2)|>
  summarize(MilkEnergy=MilkEnergy,TDiff=(MilkEnergy-elCompare$ActualMilkEnergy)/1000,PropDiff=TDiff/(elCompare$ActualMilkEnergy/1000))|>
  arrange(abs(PropDiff), descending=F) |>
  ungroup()

# Calculate overlap
eloverlap<-el |>
  group_by(CommonName, CommonName2,Age,TimeIntoLactation)|>
  filter(Age<=27)|>
  summarise_at(c("MilkEnergy"),.funs="mean") |>
  ungroup()|>
  nest_by(CommonName,CommonName2,.key="prediction") |>
  mutate(temp=list(list(pluck(prediction$MilkEnergy) |>as.vector(),eldata$ActualMilkEnergy[eldata$Age<=27])),
         overlap=flatten(list(overlapping::overlap(temp,type="1")))) |>
  unnest(overlap) |>
  select(-c(prediction, temp)) |>
  mutate(CommonNameInt=paste(CommonName, CommonName2))|>
  arrange(overlap)

# Atlantic bottlenose dolphins --------------------------------------------
ttWean=803

# Uses a mass curve that is more consistent with where I think they weighed at 2 years
tt<-data.frame(Age=seq(1,803, by=1)/365, G=0.93, g=1.25,k=0.8, Lo=102.6)|> 
  mutate(TimeIntoLactation=as.factor(round((1:ttWean/ttWean), digits=3)),Length=Lo*exp(G/g*(1-exp(-g*Age))), OffspringMass=(0.000011114*Length^3.021)*1000)|>
  left_join(PFamilyI)|>
  mutate(MilkIntake=response*(OffspringMass^0.82),MilkIntakeL=lower_ci*(OffspringMass^0.82), MilkIntakeU=upper_ci*(OffspringMass^0.82)) |>
  arrange(Population, Age) |>
  select(-c(response, se, lower_ci, upper_ci))|>
  collapse::join(PFamilyE2,how="full",multiple=T)|>
  collapse::fmutate(TimeIntoLactation=as.numeric(paste(TimeIntoLactation)), MilkEnergyD=response, 
                    MilkEnergy=MilkIntake*MilkEnergyD,MilkEnergyL=MilkIntakeL*lower_ci, MilkEnergyU=MilkIntakeU*upper_ci,Reference="Prediction",Type="Model",
                    CommonNameInt=paste(CommonName, CommonName2)) |>
  collapse::fselect(-c(response, se, lower_ci, upper_ci))|>
  filter(!is.na(MilkEnergy))|>
  droplevels()|>
  filter(CommonName==CommonName2 | CommonName2=="Bottlenose dolphin")

# Empirical data
ttdata1<-read_csv(file.path(here("data","Reddy et al 1991_Female energetics.csv")),show_col_types = FALSE)|>
  filter(x2>=0)|>
  mutate(MilkEnergy=ydiffMJ*1000*0.95, Family="Delphinidae/Phocoenidae", CommonName="Bottlenose dolphin", Reference="Reddy et al 1991") |>
  select(MilkEnergy,CommonName, Family, Reference,dolphin,x2) |>
  mutate(LactationDur=c(rep(23, 24), rep(22, 24), rep(23, 24)), TimeIntoLactation=round(x2/LactationDur,digits=3),
         Age=x2*30/365) 

ttdata2<-read_csv(file.path(here("data","Kastelein et al 2003_Female energetics.csv")),show_col_types = FALSE) |>
  mutate(MilkEnergy=ydiff*0.95, Family="Delphinidae/Phocoenidae", CommonName="Bottlenose dolphin", Reference="Kastelein et al 2003") |>
  filter(ydiff>0 & x2>=0)|>
  select(MilkEnergy,CommonName, Family, Reference,x2) |>
  mutate(dolphin=c(rep("Tt017",29),rep("Tt023",26)),
         LactationDur=c(rep(33, 29), rep(26, 26)), TimeIntoLactation=round(x2/LactationDur, digits=3),
         Age=x2*30/365)

ttdata3<-read_csv(file.path(here(),"data","Kastelein et al 2002_Female energetics.csv"),show_col_types = FALSE) |>
  mutate(MilkEnergy=ydiff*0.95, Family="Delphinidae/Phocoenidae", CommonName="Bottlenose dolphin", Reference="Kastelein et al 2002") |>
  filter(ydiff>0 & x2>=0)|>
  select(MilkEnergy,CommonName, Family, Reference,x2) |>
  mutate(dolphin=c(rep("002a",35),rep("002b",15),rep('002c',26),rep("003a",35),rep("003b",19),rep("003c",16)),
         LactationDur=c(rep(35,35), rep(14,15),rep(26,26),rep(35,35),rep(19,19),rep(16,16)),
         TimeIntoLactation=round(x2/LactationDur, digits=3), Age=x2*30/365)

ttdata<-rbind(ttdata1,ttdata2,ttdata3) |>
  select(-c(x2,Reference,Age,TimeIntoLactation,LactationDur)) |>
  mutate(Type="Empirical") |>
  pivot_wider(values_from=MilkEnergy, names_from=dolphin)

ttdataSum<-rbind(ttdata1,ttdata2,ttdata3) |>
  group_by(dolphin)|>
  dplyr::summarise(MMilkEnergy=mean(MilkEnergy),MaxMilkEnergy=max(MilkEnergy))|>
  ungroup()|>
  dplyr::summarise(MilkEnergy=mean(MMilkEnergy),MMaxMilkEnergy=mean(MaxMilkEnergy))

# Calculate overlap
ttoverlap<-tt|>
  group_by(Order, Order2, Family, Family2, CommonName, CommonName2, TimeIntoLactation)|>
  dplyr::summarise(MilkEnergy=mean(MilkEnergy)) |>
  ungroup()|>
  nest_by(CommonName,CommonName2,.key="prediction") |>
  mutate(temp=list(list(P=pluck(prediction$MilkEnergy) |>as.vector(),E1=unlist(ttdata$`tt001`),
                        E2=unlist(ttdata$`tt497`),E3=unlist(ttdata$`tt453`),E4=unlist(ttdata$`Tt017`),
                        E5=unlist(ttdata$`Tt023`),E6=unlist(ttdata$`002a`),E7=unlist(ttdata$`002b`),
                        E8=unlist(ttdata$`002c`),E9=unlist(ttdata$`003a`),E10=unlist(ttdata$`003b`),
                        E11=unlist(ttdata$`003c`))),
         overlap=list(overlapping::overlap(temp)$OVPairs)) 

ttoverlap2<-ttoverlap |>
  mutate(names=list(names(overlap)))|>
  unnest(c(overlap,names))|>
  filter(names=="P-E1" | names=="P-E2" | names=="P-E3" | names=="P-E4" | names=="P-E5"
         | names=="P-E6"| names=="P-E7"| names=="P-E8"| names=="P-E9"| names=="P-E10"| names=="P-E11")|>
  select(-c(temp, prediction)) |>
  left_join(milkEnergySum,by=c("CommonName2"="CommonName")) |>
  left_join(milkIntakeSum, by="CommonName")|>
  mutate(Alpha=case_when(overlap<0.6~0,overlap>=0.6~1,.default=NA),
         overlap2=case_when(overlap<0.6~NA,.default=overlap))|>
  group_by(Order.x, Order.y, Family.x, Family.y,CommonName, CommonName2,MIntake,MEnergy) |>
  dplyr::summarise(CountOver=sum(Alpha), MOverlap=mean(overlap2, na.rm=T)) |>
  ungroup() |>
  mutate(CountBin=cut(CountOver, breaks=c(1,3,6,10), include.lowest=T), CommonNameInt=paste(CommonName, CommonName2)) |>
  #filter(!is.na(CountBin)) |> 
  arrange(CountBin) |>
  setNames(c("OrderEnergy","OrderIntake","FamilyEnergy","FamilyIntake","CommonNameIntake","CommonNameEnergy","MIntake",
             "MEnergy","Count","MOverlap","CountBin", "CommonNameInt"))

ttoverlap3<-ttoverlap2 |>
  filter(Count>=7)

# Humpback whales --------------------------------------------

mndata<-readxl::read_xlsx(here("data","Fig.5C.data.for.Liz.xlsx")) |>
  filter(Day<91) |>
  mutate(ActualMilkEnergy=(Calf.cost.MJ.day+Lost.milk.energy.MJ.day)*1000, Age=Day)

mnMilkSpecies<-milkEnergySum |>
  filter(MaxEnergy>=19 & MaxEnergy<21)

mnWean<-10/12*365
mnAge<-(1:365)/365

mn<-data.frame(Species="Humpback whale",OffspringMass=mndata$BM.kg*1000, Age=mndata$Day) |>
  mutate(TimeIntoLactation=as.factor(round(Age/mnWean, digits=3))) |>
  left_join(PFamilyI)|>
  mutate(MilkIntake=response*(OffspringMass^0.82), MilkIntakeL=lower_ci*(OffspringMass^0.82), MilkIntakeU=upper_ci*(OffspringMass^0.82)) |>
  arrange(Population, Age) |>
  select(-c(response, se, lower_ci, upper_ci))|>
  collapse::join(PFamilyE2,how="full",multiple=T)|>
  collapse::fmutate(TimeIntoLactation=as.numeric(paste(TimeIntoLactation)), MilkEnergyD=response, 
                    MilkEnergy=MilkIntake*MilkEnergyD,MilkEnergyL=MilkIntakeL*lower_ci, MilkEnergyU=MilkIntakeU*upper_ci,Reference="Prediction",Type="Model",
                    CommonNameInt=paste(CommonName,CommonName2)) |>
  collapse::fselect(-c(response, se, lower_ci, upper_ci))|>
  filter(!is.na(MilkEnergy))|>
  droplevels() |>
  filter(CommonName==CommonName2 | CommonName2 %in% mnMilkSpecies$CommonName)

# Daily differences
mnDayDiff<-mn |>
  group_by(CommonName, CommonName2,Age,TimeIntoLactation)|>
  summarise_at(c("MilkEnergy","MilkEnergyL","MilkEnergyU"),.funs="mean") |>
  ungroup()|>
  left_join(mndata) |>
  mutate(DailyDiff=MilkEnergy-ActualMilkEnergy, DailyDiffProp=DailyDiff/ActualMilkEnergy) |>
  group_by(CommonName, CommonName2)|>
  dplyr::summarise(mMilkEnergy=mean(MilkEnergy),MDiff=mean(abs(DailyDiff))/1000, MPropDiff=mean(abs(DailyDiffProp)), MinPropDiff=min(DailyDiffProp), MaxPropDiff=max(DailyDiffProp))  |>
  arrange(abs(MPropDiff))

mnCompare<-mndata |>
  summarise_at(c("ActualMilkEnergy"),"sum")

# Total differences
mnTotalDiff<-  mn|>
  group_by(CommonName, CommonName2,Population)|>
  summarise_at(c("MilkEnergy","MilkEnergyL","MilkEnergyU"),.funs="sum") |>
  ungroup(Population) |>
  summarise_at(c("MilkEnergy","MilkEnergyL","MilkEnergyU"),.funs="mean") |>
  ungroup()|>
  group_by(CommonName, CommonName2) |>
  summarize(TDiff=(MilkEnergy-mnCompare$ActualMilkEnergy)/1000,PropDiff=TDiff/(mnCompare$ActualMilkEnergy/1000))|>
  arrange(abs(PropDiff), descending=F)

# Calculate overlap
mnoverlap<-mn|>
  group_by(Order, Order2, Family, Family2, CommonName, CommonName2, TimeIntoLactation)|>
  dplyr::summarise(MilkEnergy=mean(MilkEnergy)) |>
  ungroup()|>
  nest_by(CommonName,CommonName2,.key="prediction") |>
  mutate(temp=list(list(pluck(prediction$MilkEnergy) |>as.vector(),mndata$ActualMilkEnergy)),
         overlap=flatten(list(overlapping::overlap(temp,type="1")))) |>
  unnest(overlap)|>
  select(-c(temp, prediction)) |>
  left_join(milkEnergySum,by=c("CommonName2"="CommonName")) |>
  left_join(milkIntakeSum, by="CommonName")|>
  ungroup()|>
  mutate(CommonNameInt=paste(CommonName,CommonName2)) |>
  setNames(c("CommonName","CommonName2","overlap","Order2","Family2","MEnergy","MaxEnergy",
             "Order","Family","MIntake","MaxIntake","CommonNameInt")) |>
  arrange(overlap)

mnoverlap2<-mn |>
  group_by(Order, Order2, Family, Family2, CommonName, CommonName2,CommonNameInt,TimeIntoLactation) |>
  dplyr::summarise_at(c("MilkEnergy"),.funs="mean")|>
  ungroup(TimeIntoLactation)|>
  dplyr::summarise_at(c("MilkEnergy"),.funs="sum")|>
  ungroup() |>
  dplyr::left_join(mnoverlap)

mnoverlap3<-mnoverlap2 |>
  filter(overlap>0.6)

# Estimate what curve could be based on milk energy and mass estimates
mnCurve<-mn |>
  select(OffspringMass, Age,Population, CommonName2, Family2,MilkEnergyD)|>
  filter(CommonName2 %in% mnMilkSpecies$CommonName) |>
  dplyr::right_join(y=mndata, by="Age")  |>
  mutate(MilkIntakegday=(ActualMilkEnergy/MilkEnergyD)/(OffspringMass^0.82)) |>
  group_by(Family2,CommonName2,Age) |>
  dplyr::summarise(MMilkIntakegday=mean(MilkIntakegday), MMilkEnergyD=mean(MilkEnergyD))|>
  ungroup()|>
  mutate(TimeIntoLactation=Age/mnWean)

mnCurveSum<-mnCurve |>
  group_by(Family2,CommonName2)|>
  summarise(meanIntake=mean(MMilkIntakegday), maxIntake=max(MMilkIntakegday), minIntake=min(MMilkIntakegday))|>
  ungroup()


# Plots -------------------------------------------------------------------

# Spotted seals
Gpl<-pl |>
  group_by(Order, Order2, Family, Family2, CommonName, CommonName2,CommonNameInt,TimeIntoLactation,Age) |>
  dplyr::summarise_at(c("MilkEnergy","MilkEnergyL","MilkEnergyU"),.funs="mean") |>
  filter(CommonNameInt %in% ploverlap$CommonNameInt[ploverlap$overlap>=0.6])|>
  ggplot(aes(x=Age, y=MilkEnergy/1000))+
  geom_line(aes(group=CommonNameInt),color="gray30",linewidth=0.5, lty=2, alpha=0.5) +
  geom_ribbon(aes(ymin=MilkEnergyL/1000, ymax=MilkEnergyU/1000, group=CommonNameInt), alpha=0.02)+
  geom_smooth(data=pldata, aes(y=ActualMilkIntake/1000), se=F, linewidth=0.5, color=familycolors[7])+
  xlab("Age (days)")+
  ylab(bquote("Milk energy intake (MJ"~day^-1*")"))+
  ggthemes::theme_few()+
  theme(legend.position="none",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(color="black",linewidth=0.25)) +
  add_phylopic(uuid="157b50ce-b330-4315-b2cd-53a0fa681d10/", x=2.5,y=80,height=6.5, fill=familycolors[7],verbose=T)+
  scale_x_continuous(expand=c(0.0,0.0), lim=c(0.95,10.05))

# Sea otter
Gel<-el |>
  group_by(Order, Order2, Family, Family2, CommonName, CommonName2,CommonNameInt,TimeIntoLactation,Age) |>
  dplyr::summarise_at(c("MilkEnergy","MilkEnergyL","MilkEnergyU"),.funs="mean") |>
  filter(Age<28 & ((CommonNameInt %in% eloverlap$CommonNameInt[eloverlap$overlap>0.6]| (CommonName=="American mink" & CommonName2=="American mink"))))|>
  ggplot(aes(x=Age, y=MilkEnergy/1000))+
  geom_line(aes(group=CommonNameInt, color=Family),linewidth=0.5, lty=2, alpha=0.5) +
  geom_ribbon(aes(ymin=MilkEnergyL/1000, ymax=MilkEnergyU/1000, group=CommonNameInt), alpha=0.02)+
  geom_smooth(data=eldata |> filter(Age<=28), aes(y=ActualMilkEnergy/1000,x=Age), se=F, linewidth=0.5, color=familycolors[16])+
  xlab("Age (days)")+
  ylab(bquote("Milk energy intake (MJ"~day^-1*")"))+
  ggthemes::theme_few()+
  theme(legend.position="none",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(color="black",linewidth=0.25)) +
  add_phylopic(uuid="0070ddbf-fdcd-4a7b-97c3-4670504dc06f", x=3.5,y=5,height=0.3, fill=familycolors[16], verbose=T)+
  scale_x_continuous(expand=c(0,0), lim=c(-0.05,27.05))+
  scale_color_manual(values=c(familycolors[16],"gray30","gray30"))

# Bottlenose dolphins
ttorig<-plyr::rbind.fill(ttdata1,ttdata2,ttdata3)
ttplot<-tt |>
  filter(CommonNameInt %in% ttoverlap3$CommonNameInt) |>
  group_by(CommonNameInt, Family,Age, TimeIntoLactation)|>
  dplyr::summarise_at(c("MilkEnergy","MilkEnergyL","MilkEnergyU"),.funs="mean") |>
  ungroup()

Gtt<-ggplot(ttplot, aes(x=Age*365,y=MilkEnergy/1000))+
  geom_line(aes(group=CommonNameInt), color="gray30",linewidth=0.5, lty=2, alpha=0.5)+
  geom_point(data=ttorig, size=0.75, aes(x=Age*365), color=familycolors[10])+
  geom_ribbon(aes(ymin=MilkEnergyL/1000, ymax=MilkEnergyU/1000, group=CommonNameInt), alpha=0.05)+
  ggthemes::theme_few()+
  scale_x_continuous(expand=c(0.01,0.01), lim=c(0,806))+
  scale_y_continuous(expand=c(0.01,0.01), lim=c(0,130))+
  xlab("Age (days)")+
  ylab(bquote("Milk energy intake (MJ"~day^-1*")"))+
  theme(legend.position="none",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(color="black",linewidth=0.25)) +
  add_phylopic(uuid="388e792c-f8fd-4bd8-ad38-e251b5244ac0", x=100,y=120,height=20, fill=familycolors[10], verbose=T)

# Humpback whales  
mnplot<-mn |>
  group_by(Order, Order2, Family, Family2, CommonName, CommonName2,CommonNameInt,TimeIntoLactation,Age) |>
  dplyr::summarise_at(c("MilkEnergy","MilkEnergyL","MilkEnergyU"),.funs="mean") |>
  filter(CommonNameInt %in% mnoverlap$CommonNameInt[mnoverlap$overlap>=0.5]) 

Gmn<- ggplot(mnplot, aes(x=Age, y=MilkEnergy/1000))+
  geom_line(aes(group=CommonNameInt), linewidth=0.5, color="gray30", lty=2, alpha=0.5)+
  geom_smooth(data=mndata, aes(y=ActualMilkEnergy/1000), color=familycolors[12], linewidth=0.5, se=F)+
  geom_ribbon(aes(ymin=MilkEnergyL/1000, ymax=MilkEnergyU/1000, group=CommonNameInt), alpha=0.02)+
  xlab("Age (days)")+
  ggthemes::theme_few()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(color="black",linewidth=0.25))+
  scale_x_continuous(expand=c(0.0,0.01), lim=c(-1,93), breaks=c(0,30,60,90))+
  ylab(bquote("Milk energy intake (MJ"~day^-1*")"))+
  add_phylopic(uuid="ce70490a-79a5-47fc-afb9-834e45803ab4", x=20,y=2000,height=300, fill=familycolors[12], verbose=T)


# Combined figure ---------------------------------------------------------

GAll<-Gpl + theme(plot.tag.position=c(0.1,1.02)) +(Gel +theme(plot.tag.position=c(0,1.02)))+(Gtt+ theme(plot.tag.position=c(0.1,1.02)))+(Gmn+ theme(plot.tag.position=c(0,1.02)))+
  plot_annotation(tag_levels="a")+
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(plot.tag = element_text(face = 'bold', size=10), axis.title=element_text(size=10), axis.text=element_text(size=10)) 

ggsave(here("figures","Fig. 3.tiff"), GAll,compression="lzw",units="in", dpi=300,width=6, height=4.8)

# Derived curve plot ------------------------------------------------------
otherSpecies<-PFamilyI |>
  group_by(Family,CommonName,TimeIntoLactation)|>
  dplyr::mutate(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci),
                TimeIntoLactation=as.numeric(paste(TimeIntoLactation))) |>
  ungroup()|>
  filter((Family=="Mustelidae" |CommonName=="Gray seal" | CommonName=="Antarctic fur seal"|
            CommonName=="Reindeer"|
            CommonName=="American black bear") & TimeIntoLactation<=max(mnCurve$TimeIntoLactation)) 

GMnCurve<-ggplot()+
  geom_line(data=subset(otherSpecies), aes(y=mResponse, x=TimeIntoLactation,group=CommonName, color=Family),linetype=3, alpha=0.65,linewidth=0.75)+
  geom_line(data=mnCurve, aes(x=TimeIntoLactation,y=MMilkIntakegday, group=CommonName2),color=pals::kelly(22)[12])+
  scale_x_continuous(label=scales::percent_format(), lim=c(0,0.3), expand=c(0.01,0.005))+
  ggthemes::theme_few()+
  xlab("Time into lactation")+
  ylab(bquote("Milk intake rate (g"~day^-1~"g Offspring"^-0.82*")"))+
  scale_color_manual(values=pals::kelly(22)[c(12,16,3,7,21)],name="")+
  theme(legend.position="none", axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(color="black",size=0.25))+
  add_phylopic(uuid="0070ddbf-fdcd-4a7b-97c3-4670504dc06f",x=0.25,y=0.865, height=0.075,fill=familycolors[16], verbose=T) +
  add_phylopic(uuid="5a5dafa2-6388-43b8-a15a-4fd21cd17594",x=0.25,y=0.25, height=0.10, fill=familycolors[21], verbose=T)+
  add_phylopic(uuid="4e95fdec-8a1c-45be-b63f-dab265d4623a",x=0.25,y=0.36, height=0.10, fill=familycolors[3], verbose=T)+
  add_phylopic(uuid="157b50ce-b330-4315-b2cd-53a0fa681d10/",x=0.25,y=0.7, height=0.095, fill=familycolors[7], verbose=T)+
  add_phylopic(uuid="72f2f997-e474-4caf-bbd5-72fc8dbcc40d",x=0.25,y=0.55,height=0.195,fill=familycolors[12], verbose=T)+
  add_phylopic(uuid="ce70490a-79a5-47fc-afb9-834e45803ab4", x=0.04,y=0.9,height=0.15, fill=familycolors[12], verbose=T)

ggsave(here("figures","Fig. 4.tiff"),GMnCurve, compression="lzw",width=5.25, height=4.5, dpi=300)


