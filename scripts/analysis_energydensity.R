
library(readxl)
library(tidyverse)
library(mgcv)
library(gratia)
library(here)


# Load data ---------------------------------------------------------------

milkenergy<-readxl::read_xlsx(file.path(here("data","Milk intake data.xlsx")), sheet="Milk energy density data",
                              col_types = c(rep("text",4), rep("numeric",9), "text","text","text"))|>
  mutate(CommonName=factor(CommonName),
         Family=case_when(Family=="Delphinidae" | Family=="Phocoenidae"~"Delphinidae/Phocoenidae",.default=Family))|>
  filter(is.na(Exclude)) |>
  setNames(c("CommonName","SciName","Family","Order","N","Age","LactationDuration","TimeIntoLactation","Milkenergykjg","MilkenergySD","MilkenergySE",
                        "MilkenergyMin","MilkenergyMax","Reference","Comments","Exclude"))

milkenergy$MilkenergySD[is.na(milkenergy$MilkenergySD) | milkenergy$MilkenergySD==0]<-0.001
milkenergy$MilkenergyMin[is.na(milkenergy$MilkenergyMin)]<--Inf
milkenergy$MilkenergyMax[is.na(milkenergy$MilkenergyMax)]<-Inf


# Summarise data ----------------------------------------------------------

# Species summary
milkenergySpeciesSum<-milkenergy |>
  group_by(CommonName)|>
  dplyr::summarise(LactCoverage=pmin(1,max(TimeIntoLactation)-min(TimeIntoLactation)),
                   N=sum(N), minTime=min(TimeIntoLactation), maxTime=max(TimeIntoLactation))

# Species summary excluding very early periods
milkenergySpeciesSum2<-milkenergy |>
  filter(TimeIntoLactation>=0.03)|>
  group_by(Order,Family,CommonName)|>
  dplyr::summarise(LactCoverage=pmin(1,pmin(max(TimeIntoLactation),1)-min(TimeIntoLactation))/0.97,
                   N=sum(N),LactDur=mean(LactationDuration),
                   minTime=min(TimeIntoLactation), maxTime=max(TimeIntoLactation))

# Family summary
milkenergySpeciesSum3<-milkenergy |>
  filter(TimeIntoLactation>=0.03)|>
  group_by(Family,CommonName)|>
  dplyr::summarise(LactCoverage=(pmin(1,pmin(max(TimeIntoLactation),1)-min(TimeIntoLactation))/0.97),
                   N=sum(N)) |>
  ungroup()|>
  group_by(Family)|>
  dplyr::summarise(MLactCoverage=mean(LactCoverage), MN=mean(N), MinN=min(N), MaxN=max(N), sdN=sd(N))


# Simulate data -----------------------------------------------------------

newdata<-list()
set.seed(543)
for (i in 1:50){
  newdata[[i]]<-milkenergy |>
    filter(N>1)|>
    rowwise()|>
    mutate(SMilkenergykjg=list(MCMCglmm::rtnorm(N, mean=Milkenergykjg, 
                                                sd=MilkenergySD,
                                                lower=MilkenergyMin,
                                                upper=MilkenergyMax))) |>
    unnest(cols=c(SMilkenergykjg)) |>
    plyr::rbind.fill(filter(milkenergy, N==1))|>
    mutate(Population=i)
  print(i)
}

milkenergyNew<-bind_rows(newdata) |>
  mutate(response=ifelse(N==1,Milkenergykjg,SMilkenergykjg)) |>
  filter(TimeIntoLactation>=0.03)


# Analysis ----------------------------------------------------------------

HFamilyOtariid<-milkenergyNew |>
  filter(Family=="Otariidae")|>
  group_by(Order,Family,Population)|>
  nest()|>
  mutate(modelG=map(data,~ mgcv::gam(response ~ s(CommonName, bs="re") + s(TimeIntoLactation, bs="tp",k=5),
                                     data = .x, method="REML", family=gaussian())))

HFamilyPhocid<-milkenergyNew |>
  filter(Family=="Phocidae")|>
  group_by(Order,Family,Population)|>
  nest()|>
  mutate(modelG=map(data,~ mgcv::gam(response ~ s(CommonName, bs="re") + s(TimeIntoLactation, bs="tp",k=5),
                                     data = .x, method="REML", family=gaussian())))
         
HFamilyUrsid<-milkenergyNew |>
  filter(Family=="Ursidae")|>
  group_by(Order,Family,Population)|>
  nest()|>
  mutate(modelG=map(data,~ mgcv::gam(response ~ s(CommonName, bs="re") + s(TimeIntoLactation, bs="tp",k=5),
                                     data = .x, method="REML", family=gaussian())))

HFamilyMustelid<-milkenergyNew |>
  filter(Family=="Mustelidae")|>
  group_by(Order,Family,Population)|>
  nest()|>
  mutate(modelG=map(data,~ mgcv::gam(response ~ s(CommonName, bs="re") + s(TimeIntoLactation, bs="tp",k=5),
                                     data = .x, method="REML", family=gaussian())))

HFamilyBovidae<-milkenergyNew |>
  filter(Family=="Bovidae" & Population==1)|>
  group_by(Order,Family,Population)|>
  nest()|>
  mutate(modelG=map(data,~ mgcv::gam(response ~ s(CommonName, bs="re") + s(TimeIntoLactation, bs="tp",k=5),
                                     data = .x, method="REML", family=gaussian())))

HFamilyCervidae<-milkenergyNew |>
  filter(Family=="Cervidae")|>
  group_by(Order,Family,Population)|>
  nest()|>
  mutate(modelG=map(data,~ mgcv::gam(response ~ s(CommonName, bs="re") + s(TimeIntoLactation, bs="tp",k=5),
                                     data = .x, method="REML", family=gaussian())))
        
HFamilyDelphinidae<-milkenergyNew |>
  filter(Family=="Delphinidae/Phocoenidae")|>
  group_by(Order,Population)|>
  nest()|>
  mutate(modelG=map(data,~ mgcv::gam(response ~ s(CommonName, bs="re") + s(TimeIntoLactation, bs="tp",k=5),
                                     data = .x, method="REML", family=gaussian())),
         Family="Delphinidae/Phocoenidae")

HFamilyCamelidae<-milkenergyNew |>
  filter(Family=="Camelidae")|>
  group_by(Order,Family,Population)|>
  nest()|>
  mutate(modelG=map(data,~ mgcv::gam(response ~ s(CommonName, bs="re") + s(TimeIntoLactation, bs="tp",k=5),
                                     data = .x, method="REML", family=gaussian())))
        
HFamilyGiraffidae<-milkenergyNew |>
  filter(Family=="Giraffidae" & Population==1)|>
  group_by(Order,Family,Population)|>
  nest()|>
  mutate(modelG=map(data,~ mgcv::gam(response ~s(CommonName, bs="re") + s(TimeIntoLactation, bs="tp",k=5),
                                     data = .x, method="REML", family=gaussian())))

HFamilyElephantidae<-milkenergyNew |>
  filter(Family=="Elephantidae"  & Population==1)|>
  group_by(Order,Family, Population)|>
  nest()|>
  mutate(modelG=map(data,~ mgcv::gam(response ~ s(CommonName, bs="re") + s(TimeIntoLactation, bs="tp",k=5),
                                     data = .x, method="REML", family=gaussian())))
         
# Save model outputs -------------------------------------------------

HFamily<-bind_rows(HFamilyPhocid, HFamilyOtariid, HFamilyBovidae,
                   HFamilyCamelidae, HFamilyCervidae, HFamilyDelphinidae,
                   HFamilyElephantidae, HFamilyGiraffidae, HFamilyUrsid, HFamilyMustelid)

saveRDS(HFamily, file=here("output","Milk energy density GAM output.rds"))


# Predictions -------------------------------------------------------------
linkTransform<-gratia::inv_link(HFamily$modelG[[1]])

PFamily<-HFamily |>
  mutate(predicted = map(modelG,~ tidygam::predict_gam(.x, values=list(TimeIntoLactation=seq(0,1,by=0.001)),tran_fun=linkTransform))) |>
  unnest(predicted) |>
  select(-c(modelG,data))|>
  mutate(across(c(response, se, lower_ci, upper_ci), round,digits=3))


fst::write_fst(PFamily,here("predictions","Milk energy density GAM predictions.fst"), compress = 100)

