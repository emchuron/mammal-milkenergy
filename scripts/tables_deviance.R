library(tidyverse)
library(mgcv)

# Load in data ------------------------------------------------------------


HFamilyE<-readRDS(file=here("output","Milk energy density GAM output.rds"))
HFamilyI<-readRDS(file=here("output","Milk intake GAM output.rds"))


# wrapper function --------------------------------------------------------

get_tidy_deviance_explained <- function(model) {
  tibble(
    metric = "Deviance Explained (%)",
    value = utils.add::dev_expl(model)
  )
}

# Table 2 - Milk energy density -----------------------------------------------------

# Deviance
energySum<-HFamilyE|>
  select(-c(data))|>
  mutate(gamhp=purrr::map(modelG, gam.hp::gam.hp),
         ExpDev=purrr::map(gamhp, pluck,"Explained.deviance"),
         IndPerc=purrr::map(gamhp, pluck,"hierarchical.partitioning") |> reshape2::melt() |> 
           filter(Var2=="I.perc(%)") |> list()) |>
  select(-c(modelG, gamhp))|>
  unnest(c(IndPerc)) |>
  select(-c(Var2,L1)) |>
  pivot_wider(names_from=Var1) |>
  unnest(ExpDev)|>
  mutate(ExpDev=ExpDev*100)|>
  pivot_longer(-c(Order,Family,Population))|>
  mutate(name=case_when(str_detect(name,"CommonName")~"CommonName", str_detect(name,"Lactation")~"TimeIntoLactation",.default=name))|>
  group_by(Order, Family,name)|>
  summarise_at("value",list(mean=mean, sd=sd, min=min,max=max)) |>
  mutate_if(is.numeric, round, 1)

# Sample sizes
energyObs<-HFamilyE |>
  mutate(Glance=purrr::map(modelG, broom::glance)) |>
  select(-c(modelG))|>
  unnest(c(Glance)) |>
  ungroup()|>
  dplyr::select(Order, Family,nobs)|>
  distinct()

# Table 2 - Milk intake -------------------------------------------------------------
intakeSum<-HFamilyI|>
  select(-c(data))|>
  filter(Family!="Mustelidae") |># No random effect
  mutate(gamhp=purrr::map(modelG, gam.hp::gam.hp),
         ExpDev=purrr::map(gamhp, pluck,"Explained.deviance"),
         IndPerc=purrr::map(gamhp, pluck,"hierarchical.partitioning") |> reshape2::melt() |> 
           filter(Var2=="I.perc(%)") |> list()) |>
  select(-c(modelG, gamhp))|>
  unnest(c(IndPerc)) |>
  select(-c(Var2,L1)) |>
  pivot_wider(names_from=Var1)|>
  unnest(ExpDev)|>
  mutate(ExpDev=ExpDev*100)|>
  pivot_longer(-c(Order,Family,Population))|>
  mutate(name=case_when(str_detect(name,"CommonName")~"CommonName", str_detect(name,"Lactation")~"TimeIntoLactation",.default=name))|>
  filter(!is.na(value))|>
  group_by(Order, Family,name)|>
  summarise_at("value",list(mean=mean, sd=sd, min=min,max=max)) |>
  mutate_if(is.numeric, round, 1)

intakeObs<-HFamilyI |>
  mutate(Glance=purrr::map(modelG, broom::glance)) |>
  select(-c(modelG))|>
  unnest(c(Glance)) |>
  ungroup()|>
  select(Family, Order,nobs)|>
  distinct()

HFamilyI |>
  filter(Family=="Mustelidae")|>
  mutate(DevExp=map(modelG,get_tidy_deviance_explained)) |>
  select(-c(data))|>
  unnest(DevExp) |>
  ungroup()|>
  select(-Population)|>
  distinct()