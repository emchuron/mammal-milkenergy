# Code to create data subset, excluding data that were provided by study authors

library(readxl)
library(tidyverse)

# Load in data ------------------------------------------------------------


milkenergy<-readxl::read_xlsx(here("data","Milk intake data.xlsx"), sheet="Milk energy density data",
                              col_types = c(rep("text",4), rep("numeric",9), "text","text","text"))|>
  mutate(CommonName=factor(CommonName),
         Family=case_when(Family=="Delphinidae" | Family=="Phocoenidae"~"Delphinidae/Phocoenidae",.default=Family))|>
  filter(is.na(Exclude)) |>
  setNames(c("CommonName","SciName","Family","Order","N","Age","LactationDuration","TimeIntoLactation","Milkenergykjg","MilkenergySD","MilkenergySE",
             "MilkenergyMin","MilkenergyMax","Reference","Comments","Exclude"))

milkintake<-readxl::read_xlsx(here("data","Milk intake data.xlsx"), sheet="Milk intake data - final") |>
  filter(is.na(Exclude) & is.na(Perinatal))|>
  mutate(Family=as.factor(Family), CommonName=as.factor(CommonName))|>
  setNames(c("CommonName","SciName","Family","Order","N","LitterSize","MaternalMass","OffspringMass","Age","LactationDuration","TimeIntoLactation",
             "MilkIntakegday","MilkIntakeSE","MilkIntakeSD","MilkIntakeMin", "MilkIntakeMax"
             ,"Method","Reference","Comments","Exclude","Checked","OffspringMassSD","BirthMass","FTDuration","Perinatal"))


# Create subsets ----------------------------------------------------------

milkenergySub<-milkenergy |>
  filter(CommonName!="Harbor porpoise")

milkintakeSub<-milkintake |>
  filter(Reference!="Crocker 1995, Crocker et al. 2001; Hooper et al 2019" &
           Reference!="Donohue et al 2002" & Reference!="Arnould et al 2003" &
           Reference!="Arnould et al. 1996" & Reference!="McDonald et al 2012")


# Save output  ------------------------------------------------------------


data.table::fwrite(milkintakeSub, here("data-sub","Milk intake data.csv"))
data.table::fwrite(milkenergySub, here("data-sub","Milk energy density data.csv"))

