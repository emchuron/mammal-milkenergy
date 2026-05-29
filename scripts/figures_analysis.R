library(tidyverse)
library(ggplot2)
library(here)
library(patchwork)
library(marginaleffects)
library(DHARMa)
library(rphylopic)

familycolors=pals::kelly(22)

# Load data ---------------------------------------------------------------

# Model and predicted output
HFamilyE<-readRDS(file=here("output","Milk energy density GAM output.rds"))
HFamilyI<-readRDS(file=here("output","Milk intake GAM output.rds"))

PFamilyI<-fst::read_fst(here("predictions","Milk intake GAM predictions.fst"))
PFamilyE<-fst::read_fst(here("predictions","Milk energy density GAM predictions.fst"))

# Original data
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


# Figure 1 - Species effect -----------------------------------------------

GEnergy<-PFamilyE |>
  filter(TimeIntoLactation>=0.03)|>
  group_by(Family,CommonName)|>
  dplyr::summarise(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci)) |>
  ungroup() |>
  mutate(Family2=factor(Family, levels=c("Bovidae","Camelidae","Cervidae","Elephantidae",
                                         "Giraffidae","Mustelidae","Ursidae","Otariidae","Phocidae","Delphinidae/Phocoenidae")),
         Type=case_when(Family=="Phocidae" | Family=="Otariidae" | Family=="Delphinidae/Phocoenidae"| CommonName=="Polar bear"~2, .default=1),
         CommonName=tidytext::reorder_within(CommonName,by=Type, within=Family))|>
  arrange(Family2,mResponse)|>
  ggplot(aes(y=mResponse, x=Family2))+
  geom_point(aes(color=Family2, group=interaction(CommonName,Family2)),size=2, position=position_dodge(width=0.5))+
  geom_linerange(aes(ymin=mlwrCI, ymax=muprCI, color=Family2,group=interaction(CommonName,Family2)),position=position_dodge(width=0.5))+
  scale_color_manual(values=familycolors[c(17,6,12,14,19,16,21,3,7,10)],name="")+
  ylab(bquote("Milk energy density (kJ"~g^-1*")"))+
  ggthemes::theme_few()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position="bottom"
        ,legend.text=element_text(size=10), axis.text.y=element_text(size=10),
        axis.title=element_text(size=10),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(color="black",size=0.25))+
  xlab(NULL)+
  coord_transform(y="log2", clip="off")+
  add_phylopic(uuid="10423604-8361-4cf1-8808-91f6d0c20164", x=1, y=2.06, height=0.35, fill=familycolors[17],verbose=T)+
  add_phylopic(name="Lama glama",x=2,y=2.06,height=0.4,fill=familycolors[6],verbose=T)+
  add_phylopic(uuid="72f2f997-e474-4caf-bbd5-72fc8dbcc40d",x=3,y=2.08,height=0.425,fill=familycolors[12],verbose=T)+
  add_phylopic(uuid="62398ac0-f0c3-48f8-8455-53512a05fbc4", x=4, y=2.06, height=0.35, fill=familycolors[14],verbose=T)+
  add_phylopic(uuid="bbce74cf-4df3-4b7d-8b1d-f5b24dd3264a",x=5,y=2.2, height=0.5,fill=familycolors[19],verbose=T) +
  add_phylopic(uuid="0070ddbf-fdcd-4a7b-97c3-4670504dc06f",x=6,y=2, height=0.175,fill=familycolors[16],verbose=T) +
  add_phylopic(uuid="5a5dafa2-6388-43b8-a15a-4fd21cd17594",x=7,y=2, height=0.25, fill=familycolors[21],verbose=T)+
  add_phylopic(uuid="4e95fdec-8a1c-45be-b63f-dab265d4623a",x=8,y=2, height=0.25, fill=familycolors[3],verbose=T)+
  add_phylopic(uuid="157b50ce-b330-4315-b2cd-53a0fa681d10",x=9,y=2, height=0.2, fill=familycolors[7],verbose=T)+
  add_phylopic(uuid="388e792c-f8fd-4bd8-ad38-e251b5244ac0",x=10,y=2, height=0.3, fill=familycolors[10],verbose=T)+
  ggpubr::geom_bracket(xmin =c(0.5,7.125), xmax = c(7.025,10.5), y.position = 20,bracket.nudge.y=5,
                       label = c("Terrestrial","Marine"),label.size = 3)

GIntake<-PFamilyI |>
  filter(TimeIntoLactation>=0.03)|>
  group_by(Family,CommonName)|>
  dplyr::summarise(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci)) |>
  ungroup() |>
  mutate(Family2=factor(Family, levels=c("Artiodactyla","Mustelidae","Ursidae","Otariidae","Phocidae")),
         Type=case_when(Family=="Phocidae" | Family=="Otariidae" | Family=="Delphinidae/Phocoenidae"| CommonName=="Polar bear"~2, .default=1),
         CommonName=tidytext::reorder_within(CommonName,by=Type, within=Family))|>
  arrange(Family2,mResponse)|>
  ggplot(aes(y=mResponse, x=Family2))+
  geom_point(aes(color=Family2, group=interaction(CommonName,Family2)),size=2, position=position_dodge(width=0.5))+
  geom_linerange(aes(ymin=mlwrCI, ymax=muprCI, color=Family2,group=interaction(CommonName,Family2)),position=position_dodge(width=0.5))+
  coord_transform(y="log2")+
  scale_y_continuous(breaks=c(0.5,1.0,1.5))+
  scale_color_manual(values=pals::kelly(22)[c(12,16,21,3,7)],name="")+
  ylab(bquote("Milk intake rate (g"~day^-1~"g offspring"^-0.82*")"))+
  ggthemes::theme_few()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position="bottom",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(color="black",size=0.25),
        axis.title=element_text(size=10), axis.text.y=element_text(size=10))+
  xlab(NULL)+
  add_phylopic(uuid="0070ddbf-fdcd-4a7b-97c3-4670504dc06f",x=2,y=0.008, height=0.35,fill=familycolors[16]) +
  add_phylopic(uuid="5a5dafa2-6388-43b8-a15a-4fd21cd17594",x=3,y=0.008, height=0.5, fill=familycolors[21])+
  add_phylopic(uuid="4e95fdec-8a1c-45be-b63f-dab265d4623a",x=4,y=0.008, height=0.5, fill=familycolors[3])+
  add_phylopic(uuid="157b50ce-b330-4315-b2cd-53a0fa681d10",x=5,y=0.008, height=0.435, fill=familycolors[7])+
  add_phylopic(uuid="72f2f997-e474-4caf-bbd5-72fc8dbcc40d",x=1,y=0.009,height=0.9,fill=familycolors[12])+
  ggpubr::geom_bracket(xmin =c(0.5,3.1), xmax = c(3,5.25), y.position = 1.7,bracket.nudge.y=0.2,
                       label = c("Terrestrial","Marine"),label.size = 3)


GCombo<-GIntake+ GEnergy+ plot_layout(widths=c(0.8,1.4))+ plot_annotation(tag_levels = 'a') &
  theme(legend.position="none", plot.tag = element_text(face="bold",size=10))

ggsave(here("figures","Fig. 1.tiff"), GCombo, width=7.5, height=3.1,units="in", dpi=300, compression="lzw")


# Figure 2 - Temporal changes ---------------------------------------------

MIFamlimits<-milkintake |>
  filter(TimeIntoLactation>=0.03)|>
  mutate(Family2=case_when(Order=="Artiodactyla"~"Artiodactyla", .default=Family))|>
  group_by(Family2)|>
  summarise(minTime=round(min(TimeIntoLactation,na.rm=T), digits=3), maxTime=pmin(1,round(max(TimeIntoLactation,na.rm=T),digits=3))) |>
  pivot_longer(-c(Family2),values_to="TimeIntoLactation") |>
  data.table::setnames("Family2","Family")|>
  ungroup()|>
  mutate(TimeIntoLactation=factor(TimeIntoLactation)) |>
  left_join(PFamilyI|> mutate(TimeIntoLactation=factor(TimeIntoLactation)))|>
  mutate(TimeIntoLactation=as.numeric(paste(TimeIntoLactation)))|>
  group_by(Family,CommonName,TimeIntoLactation)|>
  dplyr::mutate(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci)) |>
  ungroup()

GIntake2<-PFamilyI |>
  filter(Family=="Mustelidae" | (CommonName=="Gray seal" & Population!=8) | CommonName=="Antarctic fur seal"|
           CommonName=="Reindeer"|
           CommonName=="American black bear")|>
  group_by(Family,CommonName,TimeIntoLactation)|>
  dplyr::mutate(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci)) |>
  ungroup()|>
  ggplot(aes(x=TimeIntoLactation, y=mResponse))+
  geom_ribbon(aes(ymin=mlwrCI, ymax=muprCI, fill=Family, x=TimeIntoLactation), alpha=0.2, inherit.aes=F)+
  geom_line(aes(color=Family, group=interaction(CommonName)))+
  ggthemes::theme_few() +
  scale_x_continuous(label=scales::percent_format(), lim=c(0,1), expand=c(0.01,0))+
  xlab("Time into lactation")+
  ylab(bquote("Milk intake rate (g"~day^-1~"g offspring"^-0.82*")"))+
  scale_color_manual(values=pals::kelly(22)[c(12,16,3,7,21)],name="")+
  scale_fill_manual(values=pals::kelly(22)[c(12,16,3,7,21)],name="")+
  theme(legend.position="none", axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(color="black",size=0.25))+
  geom_point(data=subset(MIFamlimits, Family=="Mustelidae" | CommonName=="Gray seal" | CommonName=="Antarctic fur seal" |
                           CommonName=="American black bear" | CommonName=="Reindeer"),aes(color=Family),size=1)+
  add_phylopic(uuid="0070ddbf-fdcd-4a7b-97c3-4670504dc06f",x=0.23,y=0.895, height=0.05,fill=familycolors[16]) +
  add_phylopic(uuid="5a5dafa2-6388-43b8-a15a-4fd21cd17594",x=0.25,y=0.25, height=0.075, fill=familycolors[21])+
  add_phylopic(uuid="4e95fdec-8a1c-45be-b63f-dab265d4623a",x=0.25,y=0.375, height=0.075, fill=familycolors[3])+
  add_phylopic(uuid="157b50ce-b330-4315-b2cd-53a0fa681d10/",x=0.25,y=0.705, height=0.065, fill=familycolors[7])+
  add_phylopic(uuid="72f2f997-e474-4caf-bbd5-72fc8dbcc40d",x=0.25,y=0.505,height=0.15,fill=familycolors[12])


MEFamlimits<-milkenergy |>
  filter(TimeIntoLactation>=0.03)|>
  group_by(Family)|>
  summarise(minTime=round(min(TimeIntoLactation), digits=3), maxTime=pmin(1,round(max(TimeIntoLactation, na.rm=T),digits=3))) |>
  pivot_longer(-c(Family),values_to="TimeIntoLactation") |>
  ungroup()|>
  left_join(PFamilyE |> mutate(TimeIntoLactation=as.numeric(paste(TimeIntoLactation))))|>
  group_by(Family,CommonName,TimeIntoLactation)|>
  dplyr::mutate(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci)) |>
  ungroup()

GEnergy2<-PFamilyE |>
  filter(CommonName=="Domestic ferret" | CommonName=="Gray seal" | CommonName=="Antarctic fur seal" | CommonName=='Bottlenose dolphin' |
           CommonName=="American black bear" | CommonName=="African elephant" | CommonName=="Reindeer"  )|>
  group_by(Family,CommonName,TimeIntoLactation)|>
  dplyr::mutate(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci)) |>
  ungroup()|>
  ggplot(aes(x=TimeIntoLactation, y=mResponse))+
  #geom_line(aes(color=Family))+
  geom_ribbon(aes(ymin=mlwrCI, ymax=muprCI, fill=Family, x=TimeIntoLactation), alpha=0.2, inherit.aes=F)+
  geom_line(aes(color=Family, group=CommonName))+
  ggthemes::theme_few() +
  scale_x_continuous(label=scales::percent_format(), lim=c(0.0,1), expand=c(0.01,0))+
  ylab(bquote("Milk energy density (kJ"~g^-1*")"))+
  xlab("Time into lactation")+
  scale_color_manual(values=familycolors[c(12,10,14,16,3,7,21)],name="")+
  scale_fill_manual(values=familycolors[c(12,10,14,16,3,7,21)],name="")+
  theme(legend.position="none",
        axis.title=element_text(size=10),
        axis.text=element_text(size=10),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(color="black",size=0.25))+
  guides(color = guide_legend(override.aes = list(linewidth = 3)))+
  geom_point(data=subset(MEFamlimits, CommonName=="Domestic ferret" |CommonName=="Gray seal" | CommonName=="Antarctic fur seal" | CommonName=='Bottlenose dolphin' |
                           CommonName=="American black bear" | CommonName=="African elephant" | CommonName=="Reindeer"),aes(color=Family),size=1)+
  add_phylopic(uuid="0070ddbf-fdcd-4a7b-97c3-4670504dc06f",x=0.115,y=5.5, height=1,fill=familycolors[16], verbose=T) + #Mustelid
  add_phylopic(uuid="4e95fdec-8a1c-45be-b63f-dab265d4623a",x=0.1,y=15.5, height=1.5, fill=familycolors[3], verbose=T)+ # Otariid
  add_phylopic(uuid="157b50ce-b330-4315-b2cd-53a0fa681d10/",x=0.1,y=18, height=1.2, fill=familycolors[7], verbose=T)+ # Phocid
  add_phylopic(uuid="72f2f997-e474-4caf-bbd5-72fc8dbcc40d",x=0.1,y=9,height=2.5,fill=familycolors[12], verbose=T)+ #Cervid
  add_phylopic(uuid="5a5dafa2-6388-43b8-a15a-4fd21cd17594",fill=familycolors[21], x=0.1, y=11.2, height=1.5, verbose=T)+ # Ursid
  add_phylopic(uuid="62398ac0-f0c3-48f8-8455-53512a05fbc4",fill=familycolors[14],x=0.1, y=3.75, height=1.75, verbose=T)+ # Elephant
  add_phylopic(uuid="388e792c-f8fd-4bd8-ad38-e251b5244ac0",x=0.1,y=6.85, height=1.5, fill=familycolors[10], verbose=T) # dolphin

GCombo2<-GIntake2 + GEnergy2 + plot_layout(axis_titles = "collect") +
  plot_annotation(tag_levels="a") &
  theme(plot.tag = element_text(face="bold", size=10),
        plot.tag.position=c(-0.03,1.0))

ggsave(here("figures","Fig. 2.tiff"), GCombo2, width=7.3, height=3,units="in",dpi=300, compression="lzw")

# Supplemental - Intake --------------------------------------------------------

gPhocidI<-PFamilyI |>
  filter(Family=="Phocidae" & TimeIntoLactation>=0.03)|>
  group_by(Family,CommonName,TimeIntoLactation)|>
  dplyr::mutate(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci)) |>
  ungroup()|>
  ggplot(aes(x=TimeIntoLactation, y=mResponse))+
  geom_point(data=milkintake|>filter(Family=="Phocidae" & TimeIntoLactation>=0.03),aes(y=MilkIntakegday/OffspringMass^0.82,size=N), color=familycolors[7],alpha=0.3)+
  geom_line(color=familycolors[7], linewidth=1, alpha=1)+
  facet_wrap(~CommonName)+
  guides(color = "none")+
  ggthemes::theme_few()+
  ylab(bquote("Milk intake rate (g"~day^-1~"g offspring"^-0.82*")"))+
  xlab("Time into lactation")+
  scale_x_continuous(label=scales::percent_format(), expand=c(0.01,0),
                     breaks=c(0,0.25,0.5,0.75,1.0))

ggsave(here("figures","S1 Fig.tiff"), gPhocidI, dpi=300, compression="lzw",width=10, height=7)

# Supplemental - Energy --------------------------------------------------------

gPhocidE<-PFamilyE |>
  filter(Family=="Phocidae" & TimeIntoLactation>=0.03)|>
  group_by(Family,CommonName,TimeIntoLactation)|>
  dplyr::mutate(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci)) |>
  ungroup()|>
  ggplot(aes(x=TimeIntoLactation, y=mResponse))+
  geom_point(data=milkenergy|>filter(Family=="Phocidae" & TimeIntoLactation>=0.03),aes(y=Milkenergykjg,size=N), color=familycolors[7],alpha=0.3)+
  geom_line(color=familycolors[7], linewidth=1, alpha=1)+
  facet_wrap(~CommonName)+
  guides(color = "none")+
  ggthemes::theme_few()+
  ylab(bquote("Milk energy density (kJ"~g^-1*")"))+
  xlab("Time into lactation")+
  scale_x_continuous(label=scales::percent_format(), expand=c(0.01,0),
                     breaks=c(0,0.25,0.5,0.75,1.0))

ggsave(here("figures","S6 Fig.tiff"), gPhocidE, dpi=300, compression="lzw",width=10, height=6)

gOtariidE<-PFamilyE |>
  filter(Family=="Otariidae" & TimeIntoLactation>=0.03)|>
  group_by(Family,CommonName,TimeIntoLactation)|>
  dplyr::mutate(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci)) |>
  ungroup()|>
  ggplot(aes(x=TimeIntoLactation, y=mResponse))+
  geom_point(data=milkenergy|>filter(Family=="Otariidae" & TimeIntoLactation>=0.03),aes(y=Milkenergykjg,size=N), color=familycolors[3],alpha=0.3)+
  geom_line(color=familycolors[3], linewidth=1, alpha=1)+
  facet_wrap(~CommonName)+
  guides(color = "none")+
  ggthemes::theme_few()+
  ylab(bquote("Milk energy density (kJ"~g^-1*")"))+
  xlab("Time into lactation")+
  scale_x_continuous(label=scales::percent_format(), expand=c(0.01,0),
                     breaks=c(0,0.25,0.5,0.75,1.0))

ggsave(here("figures","S7 Fig.tiff"), gOtariidE, dpi=300, compression="lzw",width=10, height=6)

gMustelidE<-PFamilyE |>
  filter(Family=="Mustelidae" & TimeIntoLactation>=0.03)|>
  group_by(Family,CommonName,TimeIntoLactation)|>
  dplyr::mutate(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci)) |>
  ungroup()|>
  ggplot(aes(x=TimeIntoLactation, y=mResponse))+
  geom_point(data=milkenergy|>filter(Family=="Mustelidae" & TimeIntoLactation>=0.03),aes(y=Milkenergykjg,size=N), color=familycolors[16],alpha=0.3)+
  geom_line(color=familycolors[16], linewidth=1, alpha=1)+
  facet_wrap(~CommonName)+
  guides(color = "none")+
  ggthemes::theme_few()+
  ylab(bquote("Milk energy density (kJ"~g^-1*")"))+
  xlab("Time into lactation")+
  scale_x_continuous(label=scales::percent_format(), expand=c(0.01,0),
                     breaks=c(0,0.25,0.5,0.75,1.0))

ggsave(here("figures","S8 Fig.tiff"), gMustelidE, dpi=300, compression="lzw",width=8, height=4)

gUrsidE<-PFamilyE |>
  filter(Family=="Ursidae" & TimeIntoLactation>=0.03)|>
  group_by(Family,CommonName,TimeIntoLactation)|>
  dplyr::mutate(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci)) |>
  ungroup()|>
  ggplot(aes(x=TimeIntoLactation, y=mResponse))+
  geom_point(data=milkenergy|>filter(Family=="Ursidae" & TimeIntoLactation>=0.03),aes(y=Milkenergykjg,size=N), color=familycolors[21],alpha=0.3)+
  geom_line(color=familycolors[21], linewidth=1, alpha=1)+
  facet_wrap(~CommonName)+
  guides(color = "none")+
  ggthemes::theme_few()+
  ylab(bquote("Milk energy density (kJ"~g^-1*")"))+
  xlab("Time into lactation")+
  scale_x_continuous(label=scales::percent_format(), expand=c(0.01,0),
                     breaks=c(0,0.25,0.5,0.75,1.0))

ggsave(here("figures","S9 Fig.tiff"), gUrsidE, dpi=300, compression="lzw",width=8, height=5)

gBovidE<-PFamilyE |>
  filter(Family=="Bovidae" & TimeIntoLactation>=0.03)|>
  group_by(Family,CommonName,TimeIntoLactation)|>
  dplyr::mutate(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci)) |>
  ungroup()|>
  ggplot(aes(x=TimeIntoLactation, y=mResponse))+
  geom_point(data=milkenergy|>filter(Family=="Bovidae" & TimeIntoLactation>=0.03),aes(y=Milkenergykjg,size=N), color=familycolors[17],alpha=0.3)+
  geom_line(color=familycolors[17], linewidth=1, alpha=1)+
  facet_wrap(~CommonName)+
  guides(size = "none")+
  ggthemes::theme_few()+
  ylab(bquote("Milk energy density (kJ"~g^-1*")"))+
  xlab("Time into lactation")+
  scale_x_continuous(label=scales::percent_format(), expand=c(0.01,0),
                     breaks=c(0,0.25,0.5,0.75,1.0))

ggsave(here("figures","S10 Fig.tiff"), gBovidE, dpi=300, compression="lzw",width=8, height=4)

gCamelidaeE<-PFamilyE |>
  filter(Family=="Camelidae" & TimeIntoLactation>=0.03)|>
  group_by(Family,CommonName,TimeIntoLactation)|>
  dplyr::mutate(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci)) |>
  ungroup()|>
  ggplot(aes(x=TimeIntoLactation, y=mResponse))+
  geom_point(data=milkenergy|>filter(Family=="Camelidae" & TimeIntoLactation>=0.03),aes(y=Milkenergykjg,size=N), color=familycolors[6],alpha=0.3)+
  geom_line(color=familycolors[6], linewidth=1, alpha=1)+
  facet_wrap(~CommonName)+
  guides(color = "none")+
  ggthemes::theme_few()+
  ylab(bquote("Milk energy density (kJ"~g^-1*")"))+
  xlab("Time into lactation")+
  scale_x_continuous(label=scales::percent_format(), expand=c(0.01,0),
                     breaks=c(0,0.25,0.5,0.75,1.0))

ggsave(here("figures","S11 Fig.tiff"), gCamelidaeE, dpi=300, compression="lzw",width=8, height=4)

gElephantidaeE<-PFamilyE |>
  filter(Family=="Elephantidae" & TimeIntoLactation>=0.03)|>
  group_by(Family,CommonName,TimeIntoLactation)|>
  dplyr::mutate(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci)) |>
  ungroup()|>
  ggplot(aes(x=TimeIntoLactation, y=mResponse))+
  geom_point(data=milkenergy|>filter(Family=="Elephantidae" & TimeIntoLactation>=0.03),aes(y=Milkenergykjg, size=N), color=familycolors[14],alpha=0.3)+
  geom_line(color=familycolors[14], linewidth=1, alpha=1)+
  facet_wrap(~CommonName)+
  guides(size = "none")+
  ggthemes::theme_few()+
  ylab(bquote("Milk energy density (kJ"~g^-1*")"))+
  xlab("Time into lactation")+
  scale_x_continuous(label=scales::percent_format(), expand=c(0.01,0),
                     breaks=c(0,0.25,0.5,0.75,1.0))

ggsave(here("figures","S12 Fig.tiff"), gElephantidaeE, dpi=300, compression="lzw",width=8, height=4)


gCervidaeE<-PFamilyE |>
  filter(Family=="Cervidae" & TimeIntoLactation>=0.03)|>
  group_by(Family,CommonName,TimeIntoLactation)|>
  dplyr::mutate(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci)) |>
  ungroup()|>
  ggplot(aes(x=TimeIntoLactation, y=mResponse))+
  geom_point(data=milkenergy|>filter(Family=="Cervidae" & TimeIntoLactation>=0.03),aes(y=Milkenergykjg,size=N), color=familycolors[12],alpha=0.3)+
  geom_line(color=familycolors[12], linewidth=1, alpha=1)+
  facet_wrap(~CommonName)+
  guides(color = "none")+
  ggthemes::theme_few()+
  ylab(bquote("Milk energy density (kJ"~g^-1*")"))+
  xlab("Time into lactation")+
  scale_x_continuous(label=scales::percent_format(), expand=c(0.01,0),
                     breaks=c(0,0.25,0.5,0.75,1.0))
  
ggsave(here("figures","S13 Fig.tiff"), gCervidaeE, dpi=300, compression="lzw",width=10, height=6)

gDolphinE<-PFamilyE |>
  filter(Family=="Delphinidae/Phocoenidae" & TimeIntoLactation>=0.03)|>
  group_by(Family,CommonName,TimeIntoLactation)|>
  dplyr::mutate(mResponse=mean(response),mlwrCI=mean(lower_ci), muprCI=mean(upper_ci)) |>
  ungroup()|>
  ggplot(aes(x=TimeIntoLactation, y=mResponse))+
  geom_point(data=milkenergy|>filter(Family=="Delphinidae/Phocoenidae" & TimeIntoLactation>=0.03),aes(y=Milkenergykjg,size=N), color=familycolors[10],alpha=0.3)+
  geom_line(color=familycolors[10], linewidth=1, alpha=1)+
  facet_wrap(~CommonName)+
  guides(size = "none")+
  ggthemes::theme_few()+
  ylab(bquote("Milk energy density (kJ"~g^-1*")"))+
  xlab("Time into lactation")+
  scale_x_continuous(label=scales::percent_format(), expand=c(0.01,0),
                     breaks=c(0,0.25,0.5,0.75,1.0))

ggsave(here("figures","S14 Fig.tiff"), gDolphinE, dpi=300, compression="lzw",width=9, height=4)

