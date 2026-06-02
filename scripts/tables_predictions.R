# Code to create supplemental table csv files
library(tidyverse)
library(here)


# Table S3. Milk intake predictions ---------------------------------------

PFamilyI<-fst::read_fst(here("predictions","Milk intake GAM predictions.fst"))|>
  mutate(CommonName=case_when(Family=="Mustelidae"~"American mink",.default=CommonName)) 

PFamilyISum<-PFamilyI |>
  group_by(Family,Order,TimeIntoLactation)|>
  summarise(n=length(unique(response)), N=length(unique(CommonName))) |>
  ungroup()|>
  select(-TimeIntoLactation)|>
  filter(n==N) |>
  distinct()

milkIPredict<-PFamilyI |>
  mutate(AnalysisGroup=Family)|>
  filter(!(Family %in% PFamilyISum$Family) | (Population==1 & Family %in% PFamilyISum$Family)) |>
  relocate(Order,AnalysisGroup, CommonName,Population,TimeIntoLactation,response, lower_ci, upper_ci) |>
  select(-c(se,Family)) |>
  mutate_if(is.numeric, round, digits=3)

data.table::fwrite(milkIPredict, here("predictions","Supplemental Table S3_Milk intake predictions.csv"))

# Table S4. Milk energy predictions ---------------------------------------

PFamilyE<-fst::read_fst(here("predictions","Milk energy density GAM predictions.fst")) |>
  mutate(TimeIntoLactation=factor(TimeIntoLactation))

PFamilyESum<-PFamilyE|>
  group_by(Family,Order,TimeIntoLactation)|>
  summarise(n=length(unique(response)), N=length(unique(CommonName))) |>
  ungroup()|>
  select(-TimeIntoLactation)|>
  filter(n==N) |>
  distinct()

milkEPredict<-PFamilyE |>
  mutate(AnalysisGroup=Family)|>
  filter(!(Family %in% PFamilyESum$Family) | (Population==1 & Family %in% PFamilyESum$Family)) |>
  relocate(Order,AnalysisGroup, CommonName,Population,TimeIntoLactation,response, lower_ci, upper_ci)|>
  select(-c(Family,se)) |>
  mutate_if(is.numeric, round, digits=3)

milkEPredictA<-milkEPredict |>
  filter(Order!="Carnivora")

milkEPredictC<-milkEPredict |>
  filter(Order=="Carnivora") 


data.table::fwrite(milkEPredictA, here("predictions","Supplemental Table S4_Milk energy density predictions.csv"))
data.table::fwrite(milkEPredictC, here("predictions","Supplemental Table S5_Milk energy density predictions.csv"))
