######## This script is for the analysis of data and creation of figures and outputs for the article "Consumer mobility predicts impacts of herbivory across an environmental stress gradient" ##########


## Creation: March 12, 2019

## Authors: Robert W. Lamb,
          # Franz Smith
          # Jon D. Witman

# Affiliation: Brown University

## 1. Set up the core functionality
##
# clean up
rm(list=ls())

# call to core packages for data manipulation
library(data.table)
library(fuzzyjoin)
library(bayesbio)
library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(assertr)
library(lubridate)
library(forcats)
library(broom) 
library(reshape2) 
library(stringr)
library(kSamples)

# for importing different formats
library(readr)
library(readxl)
library(openxlsx)  

# call to visualisation & output generation
library(ggplot2)
library(GGally)
#library(Cairo)
library(extrafont) 
library(RColorBrewer)
library(viridis)
library(gridExtra)

# call to community and modelling functionality
library(vegan)
library(effsize)
library(Hmisc)
library(marelac)

# point to working directory        ## -- will need to adjust for local copy. download the feeding_and_flow folder from the Dryad Repository link [LINK WILL BE ADDED UPON ARTICLE PUBLICATION]. replace the file path in the following command with the path to the feeding_and_flow folder on your computer -- ##
#setwd("~/Desktop/feeding_and_flow_test")

# Set sites of interest

sites_of_interest <- 
  c("camano_protected", 
    "camano_exposed", 
    "palmas_protected", 
    "palmas_exposed")



# Read in data sources

load("data_intermediate/flow/all_flow.rda")
load("data_intermediate/community/fish_bio.rda") # call to fish biomass data
load("data_intermediate/community/quad_summary.rda") # call to quadrat data

# call to taxonomy
taxonomy <- 
  read_excel("data_raw/concordances/taxonomy/algal_species.xlsx") %>%
  dplyr::rename(taxa = Morphospecies)

# urchin abundance from quadrats

urchins <- read_excel("data_raw/surveys/quadrat_urchins.xlsx") %>%
  mutate(taxa = "Eucidaris") %>%
  rename(abundance = Urchins) %>%
  mutate(z_score = scale(abundance, center= TRUE, scale=TRUE)) %>%
  dplyr::select(Site, taxa, abundance, z_score) %>%
  mutate(abundance = abundance*4) # get density per square meter

# urchin exposure data

urch_exp <- read_excel("data_raw/surveys/Urchin Exposure.xlsx") %>%
  filter(Site %in% sites_of_interest) %>%
  mutate(refuge = if_else(`Distance to Refuge` > 0, 1, 0)) %>%
  group_by(Site, Date, Plot) %>%
  summarise(Refuge = sum(refuge/length(`Distance to Refuge`))) %>%
  ungroup() %>%
  mutate(z_score = scale(Refuge, center= TRUE, scale=TRUE)) %>%
  mutate(taxa = "Refuge") %>%
  rename(abundance = Refuge) %>%
  dplyr::select(Site, taxa, abundance, z_score) %>%
  group_by(Site)%>%
  summarise(mean_refuge = mean(abundance, na.rm = T),
            sd_refuge = sd(abundance, na.rm = T))



## 2. Create Figure S1: Distribution of flow speeds measured at each site and wave exposure

wave_counts <- all_flow %>%
  mutate(Season = if_else(Month < 6, "warm", "cool")) %>%
  mutate(flow3D = round(flow3D, digits = 1)) %>%
  dplyr::select(flow3D, Site) %>%
  group_by(Site) %>%
  mutate(total = length(flow3D)) %>%
  ungroup() %>%
  group_by(Site, flow3D) %>%
  summarise(`Observations/minute` = length(flow3D)/20,
            total = mean(total)) %>%
  filter(Site %in% sites_of_interest) %>%
  ungroup() %>%
  mutate(Site = fct_recode(Site, `Caamaño exposed` = "camano_exposed", `Caamaño sheltered` = "camano_protected", `Palmas exposed` = "palmas_exposed", `Palmas sheltered` = "palmas_protected")) %>%
  ggplot(aes(flow3D, `Observations/minute`)) +
  geom_point(aes(color = Site), alpha = 0.4) +
  #geom_line(aes(color = Site)) +
  geom_smooth(aes(color = Site), se = F, size = 2) +
  geom_smooth(aes(alpha = Site), color = "black", se = F, size = .5) +
  scale_alpha_manual(values = c(1, 1, 1,1)) +
  scale_color_viridis(discrete = T, option = "plasma") +
  theme_bw() +
  lims(x = c(0, 400)) +
  theme(legend.position = c(.8,.8)) +
  labs(x = "Absolute flow speed (cm/s)")


## 3. Create Table S1: Kolmogorov-Smirnov tests for differences in the distribution of wave-induced flow speeds between sites and wave exposures 

# collate, summarize, and export flow data

# Flow conditions for each run of the experiment

site_flows <- all_flow %>% # collated flow speeds from all measurements at all sites
  mutate(date = paste(Year, Month, Day, sep = "-")) %>%
  mutate(date = ymd(date))

warm <- site_flows[site_flows$date >= "2016-11-01" & site_flows$date <= "2017-02-01",] %>%
  mutate(Season = "warm") # select warm season flow speeds

cool <- site_flows[site_flows$date >= "2017-06-01" & site_flows$date <= "2017-08-31",] %>%
  mutate(Season = "cool") # select cool season flow speeds

all_flows <- bind_rows(warm, cool) # bring data back together

# Tests of signfiicantly different distributions by Location, period

summ_flow4 <- all_flows %>%
  filter(Site != "camano_exposed_shallow",
         Site != "baltra_exposed") %>%
  mutate(site_season = paste(Site, Season, sep = "_"))


# Pairwise Anderson-Darling tests with bonferroni correction

siteseasons <- unique(summ_flow4$site_season)

dat <- as_tibble()

for(h in 1:length(unique(siteseasons))) {
  for(i in 1:length(unique(siteseasons))) {
    mod = ad.test(summ_flow4$flow3D[which(summ_flow4$site_season==siteseasons[h])], summ_flow4$flow3D[which(summ_flow4$site_season==siteseasons[i])])
    mod.output = as_tibble(mod$ad)
    dat1 <- data.frame(site1 = siteseasons[h], site2 = siteseasons[i], AD = mod.output$T.AD[1], P = mod.output$` asympt. P-value`[1])
    dat %<>% bind_rows(dat1)
  }
}

# add Bonferroni correction

dat %<>% filter(site1!=site2) %>% distinct() 
adjusted.p = p.adjust(dat$P, method = "bonferroni")
dat$adj.P = adjusted.p

Table_S1 = dat %>%
  filter(adj.P<0.05)

## 4. Create Table S2: Flow and temperature parameters at each site, wave exposure, and season
# summarise flow

# summarise flow
    flow_summary <- 
      all_flow %>%
      mutate(Season = if_else(Month < 6, "warm", "cool")) %>%
      group_by(Site, Season) %>% 
      summarise(#sig_wave_height  = 4*sd(rel_depth), # Calculate significant wave height
                #sig_wave_height2 = swh(rel_depth = rel_depth),
                mean_flow        = flow3D %>% mean(na.rm = TRUE),
                sd_flow          = flow3D %>%   sd(na.rm = TRUE),
                #max_flow         = flow3D %>%  max(na.rm = TRUE),
                #var_flow         = flow3D %>%  var(na.rm = TRUE),
                flow_eucidaris_inhibition = 1/(length(flow3D[which(flow3D>=40)])/(length(flow3D)/20)/60),
                flow_eucidaris_dislodgement = 1/(length(flow3D[which(flow3D>=532.6)])/(length(flow3D)/20)/60),
           	    flow_damselfish = 1/(length(flow3D[which(flow3D>=69)])/(length(flow3D)/20)/60),
            	flow_surgeonfish = 1/(length(flow3D[which(flow3D>=76)])/(length(flow3D)/20)/60),
            	flow_parrotfish = 1/(length(flow3D[which(flow3D>=83)])/(length(flow3D)/20)/60),
            	percentile_50 = quantile(flow3D, .50),
            	percentile_75 = quantile(flow3D, .75),
            	percentile_99 = quantile(flow3D, .99),
                mean_temp        = mean(Temperature),
                sd_temp        = sd(Temperature),
                time_measured    = length(flow3D)/20)

##
    ## generate table for supplement 

    table <- flow_summary %>%
      filter(Site %in% sites_of_interest,
             Site != "camano_exposed_shallow")
    
    table2 <- flow_summary %>%
      ungroup() %>%
      filter(Site %in% sites_of_interest,
             Site != "camano_exposed_shallow") %>%
      dplyr::select(-time_measured) %>%
      separate(Site, into = c("Site", "Exposure"), sep = "_")
 
      
    correlation = cor(table2[,-c(1:3, 7)], method = "pearson")
      

## 5. Create Figure 1: Mean and standard error of abundance/biomass of benthic algae and the herbivores that eat it across a wave stress gradient

# create benthic data
benthic <- 
  quad_summary %>%
  left_join(taxonomy, by = "taxa") %>%
  group_by(Site, 
           File, 
           `Functional Group`) %>%
  summarise(abundance = abundance %>% sum(na.rm = TRUE)) %>%
  ungroup() %>%
  spread(`Functional Group`, abundance) %>%
  mutate(Macroalgae = Filamentous + Foliose + `Corticated Branching/Bushy`) %>%
  rename(CCA = `Crustose Calcareous`) %>%
  dplyr::select(Site, 
                Macroalgae, 
                CCA) %>%
  gather(Macroalgae, CCA, key = "taxa", value = "abundance") %>%
  mutate(abundance = abundance * 100) %>% # get percent cover
  group_by(taxa) %>%
  mutate(z_score = abundance %>% scale(center = TRUE, 
                                       scale  = TRUE)) %>%
  ungroup()

# Specify fish families of interest 
families_of_interest <- 
  c("acanthuridae", 
    "labridae_scarini",
    "pomacentridae")

# create fish biomass
fish <- 
  fish_bio %>%
  filter(Site  %in% sites_of_interest,
         taxa %in% families_of_interest) %>%
  group_by(taxa) %>%
  mutate(z_score = abundance %>% scale(center = TRUE, 
                                     scale  = TRUE)) %>%
  ungroup() %>%
  mutate(abundance = abundance / 1000) %>% # get biomass in kg/50 m2
  dplyr::select(Site, 
                taxa, 
                abundance, 
                z_score)


## Create figure

# create taxa order
taxa_order <-
  c("CCA", 
    "Macroalgae", 
    "Surgeonfish", 
    "Parrotfish", 
    "Damselfish",
    "Urchins")

# create site order
site_order <-
  c("camano_exposed", 
    "palmas_exposed", 
    "camano_protected", 
    "palmas_protected")

#create data object

surveys <- bind_rows(benthic, fish, urchins) %>%
  filter(Site!="NA") %>%
  mutate(taxa = fct_recode(taxa, Parrotfish = "labridae_scarini", Surgeonfish = "acanthuridae", Damselfish = "Pomacentridae", Urchins = "Eucidaris")) %>%
  mutate(taxa = factor(taxa, levels = c("CCA", "Macroalgae", "Surgeonfish", "Parrotfish", "Damselfish", "Urchins"))) %>%
  mutate(Site = factor(Site, levels = c("camano_exposed", "palmas_exposed", "camano_protected", "palmas_protected"))) %>%
  filter(Site != "NA") 

# set sites of interest
sites_of_interest <-
  c("camano_protected", 
    "camano_exposed", 
    "palmas_protected", 
    "palmas_exposed")

# correct site namnes
site <- gsub("_protected", "", sites_of_interest)
site <- gsub("_exposed", "", site)

# correct exposure names
exposure <- gsub("camano_", "", sites_of_interest)
exposure <- gsub("palmas_", "", exposure)

# create site x exposure object
site_exposure <- 
  tibble(sites_of_interest, 
         site, 
         exposure) %>%
  rename(Site = sites_of_interest)

#create plot

barplot <- surveys %>%
  group_by(taxa, Site) %>%
  summarise(mean_ab = mean(abundance, na.rm = T),
            serr_ab = sd(abundance, na.rm = T)/sqrt(length(abundance)),
            n = length(abundance)) %>%
  mutate(upper = mean_ab + serr_ab,
         lower = mean_ab - serr_ab) %>%
  right_join(site_exposure, by = "Site") %>%
  filter(!is.na(taxa)) %>%
  ggplot(aes(site, mean_ab, fill = exposure)) +
  geom_bar(aes(fill = exposure),  position = position_dodge(width = 1), stat = "identity") +
  geom_errorbar(aes(ymax = upper,
                    ymin = lower),
                position = position_dodge(width = 1),
                stat = "identity",
                width = .25) +
  scale_fill_manual(values = c("grey30", "grey80")) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, size = 10, vjust = 0.15),
        axis.text = element_text(size = 12),
        axis.title = element_blank()) +
  facet_grid(taxa~., scales = "free_y") +
  labs(y = "Mean abundance +/- SE",
       x = "")

# statistical models

taxa <- unique(surveys$taxa)
for(i in 1:length(taxa)){
  dat = surveys[which(surveys$taxa==taxa[i]),]
  mod = lm(abundance ~ Site, data = dat)
  print(taxa[i])
  print(anova(mod))
  print(TukeyHSD(aov(mod)))
}

## 6 Multivariate analysis of benthic community structure (Appendix Fig. S2, Table S3

robbie_quads <- read_excel("data_raw/surveys/benthic_quadrat_algae_robbie.xlsx") %>%
dplyr::select(-TOTAL) %>%
  tidyr::gather(-Site, -Date_taken, -Observer, -Date_processed, -File, -Quadrat, -Cell, key = "taxa", value = "abundance")

maya_quads <- read_excel("data_raw/surveys/benthic_quadrat_algae_maya.xlsx") %>%
  dplyr::rename(Calcareous_polychaete = `Calcareous polychaete`) %>%
  dplyr::filter(is.na(notes)==T) %>%
  dplyr::select(-notes) %>%
  tidyr::gather(-Site, -Date_taken, -Observer, -Date_processed, -File, -Quadrat, -Cell, key = "taxa", value = "abundance") 

quads <- bind_rows(maya_quads, robbie_quads) %>%
  group_by(taxa) %>%
  mutate(Total = sum(abundance, na.rm = T)) %>%
  filter(Total > 0) %>%
  dplyr::select(-Total) %>%
  ungroup()

taxonomy <- read_excel("data_raw/concordances/taxonomy/algal_species.xlsx") %>%
  dplyr::rename(taxa = "Morphospecies")

# set NA to 0
quads$abundance[which(is.na(quads$abundance)==T)] = 0

##### multivariate analysis

# Create dissimilarity matrix using a bray-curtis approach
ord_data <- quad_summary %>%
  dplyr::select(-Total) %>%
  mutate(abundance = asin(sqrt(abundance))) %>%
  spread(taxa, abundance) %>%
  dplyr::select(-Sand_shell_rubble, -Bare_rock)

ord_data[is.na(ord_data)] = 0

mod <- metaMDS(ord_data[,-c(1:2)], distance = "bray", autotransform = F)

# get scores for species contributing to the MDS

species_scores=scores(mod, display="species")


plot = ordiplot(mod, type = "none")
points(mod, "sites", pch=21, col=as.factor(ord_data$Site))
text(mod, "species", col="blue", cex=0.9)

# Merge points with labelling data for MDS plot

mod_points=cbind(as.data.frame(ord_data[,1:2]), as.data.frame(mod$points)) 

# Create convex hull


# Use dplyr to create outlines of convex hulls around points for each site

library(plyr)
findhulls = function(mod_points) mod_points[chull(mod_points$MDS1, mod_points$MDS2), ]
df = mod_points
hull_plot = ddply(df, "Site", findhulls)

detach("package:plyr", unload=TRUE)

# Make labelling object

scores = species_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Species") %>%
  dplyr::rename(MDS1 = "NMDS1", MDS2 = "NMDS2")

# Generate MDS visual

Figure_S2 <- ggplot(mod_points, aes(x=MDS1, y=MDS2)) + 
  geom_polygon(data=hull_plot, aes(fill=Site), alpha = .7) +
  geom_point(color = "black") +
  scale_fill_viridis(discrete = T) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(data = scores, aes(MDS1, MDS2, label = Species))

# PERMANOVA style statistical analysis for significant differences between sites
perm = adonis(ord_data[,-c(1:2)] ~ ord_data$Site)

# SIMPER analysis to look at contribution of each algal species to difference between sites

site_simper = simper(ord_data[,-c(1:2)], ord_data$Site)
Table_S3 = summary(site_simper)

# BetaDisper (dispersion) analysis to look at breadth in community composition

dispersion.test <- betadisper(vegdist(ord_data[,-c(1:2)], method = "bray"), group = ord_data$Site)
range(ord_data[,-c(1:2)])

anova(dispersion.test)

permutest(dispersion.test, permutations = 99, pairwise = T)


## 7 Regression analyses of consumer abundance and biomass against mean flow speed (Fig. 2)

load("data_intermediate/plate_biomass.rda")

site_flows <- all_flow %>%
  mutate(date = paste(Year, Month, Day, sep = "-")) %>%
  mutate(date = ymd(date))

warm <- site_flows[site_flows$date >= "2016-11-01" & site_flows$date <= "2017-02-01",] %>%
  mutate(Season = "warm")

cool <- site_flows[site_flows$date >= "2017-06-01" & site_flows$date <= "2017-08-31",] %>%
  mutate(Season = "cool")

all_flows <- bind_rows(warm, cool) 

summ_flows <- all_flows %>%
  # filter(flow3D<=100) %>%
  group_by(Site, Season) %>% 
  summarise(mean_flow = mean(flow3D, na.rm = T),
            sd_flow = sd(flow3D, na.rm = T),
            mean_temp = mean(Temperature),
            sd_temp = sd(Temperature)) %>%
  gather(-Site, -Season, key = "variable", value = "value")


# algal biomass for each run of the experiment

biomass_mod <- biomass %>%
  dplyr::filter(`Biomass category` %in% c("Filamentous", "Foliose", "Thick_foliose", "Branch")) %>%
  mutate(Site = tolower(Locality),
         Season = tolower(Season)) %>%
  group_by(Site, Season, Treatment, Plot, Plate) %>%
  summarise(macroalgae=sum(Biomass, na.rm =T)) %>%
  ungroup() %>%
  mutate(Treatment = fct_recode(Treatment, `u+f+` = "Open", `u+f-` = "Roof", `u-f+` = "Fence", `u-f-` = "Cage")) %>%
  dplyr::filter(Treatment == "u-f-") %>%
  group_by(Site, Season) %>%
  summarise(mean_macro = mean(macroalgae, na.rm =T)) %>%
  mutate(log_macro = log10(mean_macro+1)) %>%
  right_join(summ_flows)

# Plot

mean_plot <- biomass_mod %>%
  dplyr::filter(variable == "mean_flow") %>%
  ggplot(aes(value, log_macro)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F, color = "black") +
  theme_bw() +
  labs(x = "Mean flow speed (cm/s)",
       y = "Mean macroalgal biomass (g)") +
  theme(axis.text.x  = element_blank(),
        axis.text.y  = element_text(colour = "black", size = 11),
        axis.title.y = element_text(colour = "black", size = 13),
        axis.title.x = element_blank()) +
  # geom_text(aes(x = 15, y = 0, label = "y = 0.022x - 0.143; p = 0.002; R2 = 0.8")) +
  xlim(8, 42)



# fish biomass for each run of the experiment

# in-house fish traits
lab_ids <- 
  paste0("data_raw/concordances/taxonomy/",
         "fish_ids.xlsx") %>%
  read_excel()

# gloabl fish trait database from David Mouillot
gaspar <- 
  paste0("data_raw/concordances/traits/",
         "global_traits_species.xlsx") %>%
  read_excel() %>%
  mutate(genus_species = capitalize(as.character(genus_species)))

fish <- read_excel("data_raw/surveys/herbivory_flow_fish_transects.xlsx") %>%
  gather(-Site, -Date, -Observer, -Depth, -Temperature, -Water.movement, -Substrate, -Transect, -Block, -Transect.Length, -Block.Length, -Species, -Sick, -Site.Exposure, key = Size, value = Abundance) %>%
  mutate(Size      = gsub("X", "", Size),
         Size      = as.numeric(Size),
         Abundance = ifelse(is.na(Abundance), 0, Abundance)) %>%
  mutate(Abundance = Abundance*10/Block.Length) %>%
  right_join(lab_ids, fish, by="Species", all=TRUE) %>%
  right_join(gaspar, by="genus_species", all=TRUE)


# Extract biomass
genera_of_interest <- c("Prionurus", "Scarus", "Stegastes", "Ophioblennius", "Thalassoma", "Halichoeres", "Holacanthus", "Microspathodon")
families_of_interest <- c("acanthuridae", "labridae_scarini", "pomacentridae")
methods <- c("lm", "lm", "gam", "gam", "gam", "gam", "gam", "gam")

fish_bio <- fish %>%
  filter(Species != "NA") %>%
  dplyr::select(-Species, -`Max weight (kg)`, -`Common TL (cm)`, -`Max TL`, -`Min Depth (m)`, -`Max depth (m)`, -`Usual Min depth`, -`Usual max depth`, -`Max lat`, - `Min lat`, -`Max long`, -`Min long`, -Sick) %>%
  mutate(biomass=Abundance*(constant_a*Size^allometric_coefficient_b)*10/`Block.Length`) %>%
  separate(genus_species,
           into = c("genus", "species"),
           sep = "_") %>%
  filter(family %in% families_of_interest) %>%
  group_by(Site, Date, Transect, Block, family) %>%
  summarise(Biomass = sum(biomass, na.rm = T)) %>%
  #Abundance = sum(Abundance, na.rm = T)) %>%
  #complete(family, nesting(Site, Date, Transect, Block), fill = list(Abundance = 0, Biomass = 0)) %>%
  spread(family, Biomass, fill =0) %>%
  gather(-Site, -Date, -Transect, - Block, key = family, value = Biomass) %>%
  group_by(Site, Date, family) %>%
  summarise(mean_bio = mean(Biomass, na.rm = T))


fish_bio <- fish_bio %>%
  filter(Site %in% sites_of_interest) %>%
  rename(taxa = family,
         abundance = mean_bio) %>%
  mutate(abundance = abundance/1000) %>% # get biomass in kg/50 m2 
  rename(date = Date)

fwarm <- fish_bio[fish_bio$date >= "2016-11-01" & fish_bio$date <= "2017-02-01",] %>%
  mutate(Season = "warm")

fcool <- fish_bio[fish_bio$date >= "2017-06-01" & fish_bio$date <= "2017-08-31",] %>%
  mutate(Season = "cool")

fish_plot <- bind_rows(fwarm, fcool) %>%
  filter(taxa %in% families_of_interest) %>%
  group_by(Site, Season, taxa) %>%
  summarise(mean_fish = mean(abundance, na.rm = T)) %>%
  mutate(mean_fish = log10(mean_fish+1)) %>%
  right_join(summ_flows) %>%
  filter(variable == "mean_flow") %>%
  #mutate(sig = as.numeric(fct_recode(taxa, "1" = "acanthuridae", "1" = "labridae_scarini", "0" = "pomacentridae")),
  #       sig2 = as.character(fct_recode(taxa, black = "acanthuridae", black = "labridae_scarini", grey10 = "pomacentridae"))) %>%
  filter(Site != "camano_exposed_shallow",
         !is.na(taxa)) %>%
  ggplot(aes(x = value, y = mean_fish)) +
  geom_point(size = 2) +
  geom_smooth(se = F, method = "lm", color = "black", aes(linetype = taxa)) +
  #scale_color_viridis(discrete = T, end = .7) +
  facet_grid(taxa~., scales = "free_y") +
  theme_bw() +
  scale_linetype_manual(values = c(1, 1, 0)) +
  labs(x="Mean flow (cm/s)",
       y = "Mean fish biomass (kg/50m2)") +
  #geom_text(aes(x = 15, y = 0, label = c(rep("y = 0.03x -0.25; p = 0.033; R2 = 0.48", 8), rep("y = 0.02x -0.09; p = 0.016; R2 = 0.59", 8), rep("", 8)))) +
  theme(strip.text = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))+
  xlim(5, 41)


# models

fishes <- bind_rows(fwarm, fcool) %>%
  filter(taxa %in% families_of_interest) %>%
  group_by(Site, Season, taxa) %>%
  summarise(mean_fish = mean(abundance, na.rm = T)) %>%
  mutate(mean_fish = log10(mean_fish+1)) %>%
  right_join(summ_flows) %>%
  filter(variable == "mean_flow") 

pomacentridae_mod <- lm(mean_fish~value, data = fishes[which(fishes$taxa=="pomacentridae"),])
summary(pomacentridae_mod)

parrotfish_mod <- lm(mean_fish~value, data = fishes[which(fishes$taxa=="labridae_scarini"),])
summary(parrotfish_mod)

acanthuridae_mod <- lm(mean_fish~value, data = fishes[which(fishes$taxa=="acanthuridae"),])
summary(acanthuridae_mod)

## urchins from field census data

load("data_intermediate/community/urchin_abundance.rda")

urchin_plot <- urchins %>%
  filter(Site != "camano_exposed_shallow") %>% # remove extra survey data points
  ggplot(aes(x = value, y = log10(mean_ab+1.1))) +
  theme_bw()+
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F, color = "black") +
  theme_bw() +
  labs(x="Mean flow (cm/s)",
       y = "Mean urchin abundance (individuals/5m2)") +
  #  geom_text(aes(x = 15, y = 0, label = "y = -0.01x + 1.62; p < 0.001; R2 = 0.31")) +
  theme(axis.title.x = element_text(size = 13),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13)) +
  xlim(8, 42)


# model

urchmod <- lm(mean_ab~value, data = urchins %>%filter(Site != "camano_exposed_shallow"))
summary(urchmod)

## 9. Rapid algal feeding on Ulva spp. trials under different wave exposure and herbivore access treamtents (Figure 3)

#Bring in the data 

ulva <- read.csv("data_raw/experiments/ulva_sanduches/rapid_algal_assays_08_23_17_complete.xlsx") %>%
  #filter(Site == "Las Palmas") %>%
  mutate(mass_lost = Before_weight_total-After_weight_ulva-After_weight_mesh, perc_mass_lost = mass_lost/(Before_weight_total-After_weight_mesh)*100, perc_mass_rem = 100-mass_lost/(Before_weight_total-After_weight_mesh)*100)

# barplot

Figure_3 <- ulva %>%
  mutate(Treatment = fct_recode(Treatment, `+all` = "Control", `+urchins` = "Roof", `+fish` = "Fence", `-all` = "Exclusion")) %>%
  mutate(Treatment = factor(Treatment, levels = c("+all", "+fish", "+urchins", "-all"))) %>%
  rename(Exposure = Location) %>%
  group_by(Treatment, Site, Exposure) %>%
  summarise(mean_macro = mean(perc_mass_rem, na.rm =T),
            serr_macro = sd(perc_mass_rem, na.rm=T)/sqrt(length(perc_mass_rem))) %>%
  mutate(upper = mean_macro+serr_macro,
         lower = mean_macro-serr_macro) %>%
  ggplot(aes(x = Treatment, y = mean_macro, fill = Treatment, ymax = upper, ymin = lower)) +
  geom_bar(aes(fill = Treatment), stat = "identity", color = "black") +
  facet_grid(Site~Exposure) +
  geom_linerange() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("white", "grey40", "grey80", "black")) +
  theme_bw() +
  labs(x = "", 
       y= "Ulva % biomass remaining +/- SE") +
  theme(text = element_text(size=16),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")


# Procure summary stats
# Statistical analyses

# Tests for ANOVA assumptions

# Homogeneity of variance

# Log transform to fit model assumptions

mod=lm(perc_mass_lost~Site*Location*Treatment, dat=ulva)
summary(mod)
anova(mod)            

TukeyHSD(aov(mod))

##10. Response of algal community assembly to different herbivore exclusion treatments after 2 months (Figure 4, Table S5)

load("data_intermediate/plate_biomass.rda")

# summary date exploration


sites = c("Palmas", "Palmas","Camano", "Camano")
seasons = c("warm", "cool", "warm", "cool")
start_date = c("2016-11-16", "2017-06-15","2016-11-17","2017-06-14")
end_date = c ("2017-01-17","2017-08-23","2017-01-18","2017-08-25")
dates <- tibble(sites, seasons, start_date, end_date) %>%
  mutate(elapsed_time = int_length(interval(start_date, end_date))/60/60/24)

# Exploratory plotting and getting to know the data

# models

biomass_mod <- biomass %>%
  filter(`Biomass category` %in% c("Filamentous", "Foliose", "Thick_foliose", "Branch")) %>%
  group_by(Site, Season, Exposure, Treatment, Plot, Plate) %>%
  summarise(macroalgae=sum(Biomass, na.rm =T)) %>%
  ungroup() %>%
  mutate(Treatment = fct_recode(Treatment, `+all` = "Open", `+urchins` = "Roof", `+fish` = "Fence", `-all` = "Cage")) %>%
  mutate(Treatment = factor(Treatment, levels = c("+all", "+urchins", "+fish", "-all")))

#set color by season
season_color <- c("blue", "blue", "blue", "blue", "red", "red", "red", "red")

Figure_4 <- biomass_mod %>%
  group_by(Treatment, Site, Season, Exposure) %>%
  summarise(mean_macro = mean(macroalgae, na.rm =T),
            serr_macro = sd(macroalgae, na.rm=T)/sqrt(length(macroalgae))) %>%
  mutate(upper = mean_macro+serr_macro,
         lower = mean_macro-serr_macro) %>%
  ggplot(aes(x = Treatment, y = mean_macro, fill = Treatment)) +
  geom_bar(aes(fill = Treatment), stat = "identity", color = "black") +
  #facet_grid(Exposure~Season*Site, scales = "free_y") +
  facet_wrap(~Season*Site*Exposure, ncol = 2, scales = "free_y") +
  geom_linerange(aes(ymax = upper,
                     ymin = lower)) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("white", "grey40", "grey80", "black")) +
  theme_bw() +
 # theme(strip.background = element_blank(),
  #      strip.text = element_blank()) +
  labs(x = "", 
       y= expression(paste('Macroalgal biomass (g*169cm',""^-2,') +/-SE'))) +
  theme(text = element_text(size=16),
        legend.position = "none")



#models

warm <- biomass_mod %>%
  filter(Season == "Warm")

warm_mod <- lm(log10(macroalgae+1.1) ~ Site*Exposure*Treatment, data = warm)
anova(warm_mod)
TukeyHSD(aov(warm_mod))

cool <- biomass_mod %>%
  filter(Season == "Cool")

cool_mod <- lm(log10(macroalgae+1.1) ~ Site*Exposure*Treatment, data = cool)
anova(cool_mod)
TukeyHSD(aov(cool_mod))


# or all together

sites = unique(biomass_mod$Site)
seasons = unique(biomass_mod$Season)
for (i in 1:length(seasons)) {
  dat = biomass_mod[which(biomass_mod$Season == seasons[i]),]
  for(j in 1:length(sites)){
    dat2 = dat[which(dat$Site == sites[j]),]
    mod = lm(log10(macroalgae+1.1) ~ Exposure*Treatment, data = dat2)
    print(seasons[i])
    print(sites[j])
    print(anova(mod))
    print(TukeyHSD(aov(mod)))
  }
}


## 11. Analysis of contribution of each fish species observed biting the benthos (Appendix Table S6)
## This script is to collate all feeding observations from video analysis for the purpose of examining the total herbivory pressure by each fish species for the 2-month herbivore exclusion experiments ##

## read in concordance files

taxonomy <- read_excel("data_raw/concordances/taxonomy/fish_ids.xlsx") %>%
  mutate(`Common Name` = as.character(`Common Name`))

# point to data locale
data_locale <- "data_intermediate/boris/"

# read in intermediate data object of collated feeding observations from BORIS analyses

load(paste(data_locale, "feeding_data_from_boris.rda", sep = ""))

## summarize to get average values per species per site, exposure and sampling period


conv.data <- feeding_data_from_boris %>%
  group_by(Site, Exposure, Date, Total_seconds, Species) %>%
  summarise(total_bites = length(Bites)/unique(Total_seconds)*60) %>%
  ungroup() %>%
  #rename(Species = Species2) %>%
  mutate(Species = as.factor(str_trim(Species))) %>%
  mutate(Species = fct_recode(Species, BODI = "Mexican hogfish",
                              PRLA = "Surgeonfish",
                              STBE = "Ringtail damselfish",
                              HADI = "Chameleon wrasse",
                              OPST = "Fanged blenny",
                              MIDO = "Giant damselfish",
                              THLU = "Rainbow wrasse",
                              HANI = "Spinster wrasse",
                              ANSC = "Peruvian grunt",
                              HOPA = "King angelfish",
                              SCRU = "Bicolored parrotfish",
                              SCGH = "Bluechin parrotfish",
                              SPAN = "Bullseye puffer",
                              PLAZ = "Sabertoothed blenny",
                              SPAN = "Sphoeroides annulatus",
                              STBE = "Stegastes beebei",
                              HANI = "Halichoeres nicholsi",
                              BODI = "Bodianus diplotaenia")) %>%
  left_join(taxonomy, by = "Species") %>%
  select(Site, Exposure, Species, `Scientific Name`, total_bites) 


herbivory <- read_csv("data_raw/herbivory/compiled_vid_herbivory_11_09.csv", col_names = T) %>%
  as_tibble() %>%
  group_by(Site, Exposure, Segment, `Video number`, Species) %>%
  summarise(total_bites = sum(Bites)) %>%
  ungroup() %>%
  spread(Species, total_bites, fill = 0) %>%
  gather(-Site, -Exposure, -Segment, -`Video number`, key = Species, value = total_bites) %>%
  mutate(total_bites = total_bites/10) %>%  # to get bites per minute
  left_join(taxonomy, by = "Species") %>%
  select(Site, Exposure, Species, `Scientific Name`, total_bites) %>%
  mutate(Site = substr(Site, start = 1, stop = 6))



## merge data and analyze

sites_to_analyze <- c("camano", "palmas")

# Extract biomass

fish_bio <- fish %>%
  filter(Species != "NA") %>%
  dplyr::select(-Species, -`Max weight (kg)`, -`Common TL (cm)`, -`Max TL`, -`Min Depth (m)`, -`Max depth (m)`, -`Usual Min depth`, -`Usual max depth`, -`Max lat`, - `Min lat`, -`Max long`, -`Min long`, -Sick) %>%
  mutate(biomass=Abundance*(constant_a*Size^allometric_coefficient_b)*10/`Block.Length`) %>%
  filter(diet_gaspar == "hd") %>%
  group_by(Site, Date, Transect, Block, `Scientific Name`) %>%
  summarise(Biomass = sum(biomass, na.rm = T)) %>%
  #Abundance = sum(Abundance, na.rm = T)) %>%
  #complete(family, nesting(Site, Date, Transect, Block), fill = list(Abundance = 0, Biomass = 0)) %>%
  spread(`Scientific Name`, Biomass, fill =0) %>%
  gather(-Site, -Date, -Transect, - Block, key = `Scientific Name`, value = Biomass) %>%
  group_by(Site, `Scientific Name`) %>%
  summarise(mean_bio = mean(Biomass, na.rm = T)) %>%
  filter(Site != "NA") %>%
  ungroup()

# get site and exposure into separate columns

locality = fish_bio$Site
site = gsub("_protected", "", locality)
site = gsub("_exposed", "", site)

exposure = gsub("camano_", "", locality)
exposure = gsub("palmas_", "", exposure)

# set sites of interest

sites_of_interest <- c("camano", "palmas")

fish_bio_1 <- fish_bio %>%
  ungroup() %>%
  dplyr::select(-Site) %>%
  mutate(Site = site,
         Exposure = exposure) %>%
  filter(Site %in% sites_of_interest) %>%
  mutate(Site = capitalize(Site),
         Exposure = capitalize(Exposure))

## To create a table of mean and standard deviation of bites for each species and site (Table 2 in manuscript)

bites_table <- bind_rows(herbivory, conv.data) %>%
  filter(Site %in% (sites_to_analyze)) %>%
  mutate(Site = capitalize(Site)) %>%
  complete(`Scientific Name`, nesting(Site, Exposure), fill = list(total_bites = 0)) %>%
  #total_bites = log10(total_bites+1)) %>%
  group_by(Site, Exposure, `Scientific Name`) %>%
  summarise (mean_bites = mean(total_bites, na.rm = T),
             sd_bites = sd(total_bites, na.rm = T)) %>%
  left_join(fish_bio_1) %>%
  replace_na(list(mean_bio = 0,
                  sd_bites = 0)) %>%
  #filter(Site == "Palmas", Exposure == "Protected") %>% ## currently this code requires changing the site and exposure parameters for each column of the table, can be coded in to produce the whole table in one run
  ungroup() %>%
  group_by(Site, Exposure) %>%
  mutate(mean_sd = paste(round(mean_bites, 2), " ", "(", round(sd_bites, 2), ")", sep = "")) %>%
  ungroup() %>%
  mutate(location = paste(Site, Exposure)) %>%
  dplyr::select(location, `Scientific Name`, mean_sd) %>%
  spread(key = location, value = mean_sd)

## 12. Analysis of flow speed distributions when fish were feeding on the benthos (Figure 5, Appendix Table S7)

load("data_intermediate/feeding_and_flow/combined_feeding_and_flow.no_summary.rda")

#visualize each species across all flow observations

# set cut off for flow data
flow_cut <- 100

# set palette options
n_palette <- combined_feeding_and_flow.no_summary$Species %>% unique() %>% length()

## Look at histogram with global flow values as comparison


Figure_5 <- combined_feeding_and_flow.no_summary %>%
  filter(!is.na(vertical)) %>%
  #rename(Site = Site.x, Exposure = Exposure.x, Date = Date.x, true_time = true_time.x, true_sec = true_sec.x) %>%
  dplyr::select(Site, Exposure, Date, true_time, true_sec, `onshore-offshore`, vertical, `along-shore`, flow3D, Species, Bites, Subject_number) %>%
  mutate(`onshore-offshore` = abs(`onshore-offshore`), `along-shore` = abs(`along-shore`), vertical = abs(vertical)) %>%
  group_by(true_sec) %>%
  mutate(max_flow = max(`onshore-offshore`, vertical, `along-shore`)) %>%
  ungroup() %>%
  filter(max_flow<=100) %>%
  dplyr::select(-max_flow) %>%
  mutate(Species = fct_recode(Species, Parrotfish = "Bluechin parrotfish", Parrotfish = "Bicolored parrotfish")) %>%
  filter(Species != "Sabertoothed blenny",
         Species != "Bullseye puffer",
         Species != "Mexican hogfish") %>%
  group_by(Species, Site, Exposure, Date) %>%
  mutate(number = max(Subject_number, na.rm = T)) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(total_biters = sum(unique(number))) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(med3D = median(flow3D)) %>%
  arrange(med3D) %>%
  ungroup() %>%
  mutate(Species = factor(Species, levels = unique(Species))) %>%
  gather(`onshore-offshore`, vertical, `along-shore`, flow3D, key = "Dimension", value = "Velocity") %>%
  mutate(Dimension = factor(as.character(Dimension), levels = c("flow3D", "onshore-offshore", "along-shore", "vertical"))) %>%
  uncount(Bites)  %>%
  group_by(Dimension, Species) %>%
  mutate(medians = median(Velocity, na.rm = T),
         average = mean(Velocity, na.rm = T)) %>%
  ungroup()

# Reshape data

histogram_data <- combined_feeding_and_flow.no_summary %>%
  filter(!is.na(vertical)) %>%
  #rename(Site = Site.x, Exposure = Exposure.x, Date = Date.x, true_time = true_time.x, true_sec = true_sec.x) %>%
  dplyr::select(Site, Exposure, Date, true_time, true_sec, `onshore-offshore`, vertical, `along-shore`, flow3D, Species, Bites, Subject_number) %>%
  mutate(`onshore-offshore` = abs(`onshore-offshore`), `along-shore` = abs(`along-shore`), vertical = abs(vertical)) %>%
  group_by(true_sec) %>%
  mutate(max_flow = max(`onshore-offshore`, vertical, `along-shore`)) %>%
  ungroup() %>%
  filter(max_flow<=100) %>%
  dplyr::select(-max_flow) %>%
  mutate(Species = fct_recode(Species, Parrotfish = "Bluechin parrotfish", Parrotfish = "Bicolored parrotfish")) %>%
  filter(Species != "Sabertoothed blenny",
         Species != "Bullseye puffer",
         Species != "Mexican hogfish") %>%
  group_by(Species, Site, Exposure, Date) %>%
  mutate(number = max(Subject_number, na.rm = T)) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(total_biters = sum(unique(number))) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(med3D = median(flow3D)) %>%
  arrange(med3D) %>%
  ungroup() %>%
  mutate(Species = factor(Species, levels = unique(Species))) %>%
  gather(`onshore-offshore`, vertical, `along-shore`, flow3D, key = "Dimension", value = "Velocity") %>%
  mutate(Dimension = factor(as.character(Dimension), levels = c("flow3D", "onshore-offshore", "along-shore", "vertical"))) %>%
  uncount(Bites)  %>%
  group_by(Dimension, Species) %>%
  mutate(medians = median(Velocity, na.rm = T),
         average = mean(Velocity, na.rm = T)) %>%
  ungroup()

# Set up plotting settings

Species <- unique(histogram_data$Species)
latin_names <- c("Scarus spp.", "Stegastes beebei", "Halichoeres dispilus", "Stegastes arcifrons", "Halichoeres nicholsi", "Holacanthus passer", "Thalassoma lutescens", "Microspathodon dorsalis", "Ophioblennius steindachneri", "Prionurus laticlavius")
group <- c("Low", "Low", "Low", "Low", "Mid", "Mid", "High", "High", "High", "High")
family <- c("Labridae-Scarini", "Pomacentridae", "Labridae", "Pomacentridae", "Labridae", "Pomacanthidae", "Labridae", "Pomacentridae", "Blenniidae", "Acanthuridae")
specs <- data_frame(family, latin_names, Species, group)



all_dat <- histogram_data %>%
  left_join(specs, by = "Species") %>%
  filter(Dimension!="flow3D") %>%
  mutate(group = factor(group, levels = unique(group)),
         Species = factor(Species, levels = unique(Species)),
         Dimension = factor(Dimension, levels = unique(Dimension)))



# PLot histogram

# set up breaks for y-axes

myBreaks <- function(x){
  breaks <- c(round(min(x)), round(median(x)), round(max(x)))
  names(breaks) <- attr(breaks,"labels")
  breaks
}

histogram_plot <- all_dat %>%
  filter(Dimension == "onshore-offshore",
         Velocity <= 150) %>%
  arrange(Species, Dimension, Velocity) %>%
  group_by(Species, Dimension) %>%
  mutate(Bites = n()) %>%
  mutate(flow_total_bites = cumsum(Bites)) %>%
  mutate(rel_flow_total_bites = flow_total_bites/max(flow_total_bites)) %>%
  mutate(med3D = median(Velocity),
         n95_3D = Velocity[which.min(abs(rel_flow_total_bites-.95))]) %>%
  ungroup() %>%
  ggplot(aes(Velocity)) +
  geom_histogram(aes(fill = group),binwidth = 1) +
  #stat_smooth() +
  #facet_grid(Species~Dimension, scales = "free_y") +
  facet_wrap(Species~Dimension, scales = "free_y", nrow = 10) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 150)) +
  scale_fill_manual(values = viridis_pal()(4)) +
  theme(strip.background = element_blank(),
        #strip.text.y     = element_text(angle = 0),
        strip.text = element_blank(),
        legend.position  = "none",
        plot.title       = element_text(size  = 16,
                                        hjust = 0.5),
        axis.title       = element_text(size = 14)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
  labs(x = "Flow speed (cm/s)",
       y = "Bites")  



# Produce range of times that flow conditions were measured during bite observations

flow_times <- histogram_data %>%
  group_by(Site,Exposure,Date) %>%
  summarise(max_true_time = max(true_time),
            min_true_time = min(true_time),
            max_true_sec = max(true_sec),
            min_true_sec = min(true_sec))

# Bring in flow data

# point to flow locale
flow_locale <- "data_intermediate/flow/"

# call to data
load(paste0(flow_locale, "all_flow.rda"))

# create time field for flow data
all_flow %<>%
  mutate(Date = paste0(Year, "-", 
                       Month, "-", 
                       Day) %>% as.Date()) %>%
  mutate(dateTime = paste0(Year, "-", 
                           Month, "-", 
                           Day, " ", 
                           Hour,   ":", 
                           Minute, ":",
                           Second) %>% as_datetime()) %>%
  dplyr::select(Site,
                Date,
                dateTime,
                `vertical`         = VelocityX,
                `onshore-offshore` = VelocityY,
                `along-shore`      = VelocityZ,
                flow3D,
                depth)

# create site and exposure
all_flow %<>%
  mutate(Site = gsub("exposed_shallow", "exposed-shallow", Site)) %>%
  separate(Site,
           into = c("Site", "Exposure"),
           sep  = "_") %>%
  mutate(Exposure = Exposure %>% str_to_title())



flow <- all_flow %>%
  rename(true_time = dateTime) %>%
  mutate(true_sec = as.numeric(seconds(ymd_hms(true_time)))) %>%
  dplyr::select(-vertical, -`along-shore`, -flow3D, -depth)


# Anderson-Darling tests for differences in empirical distribution function with bonferroni correction

spec <- unique(histogram_data$Species)


dat <- as_tibble()

  data = histogram_data %>% filter (Dimension == "flow3D")
  for(i in 1:length(unique(spec))) {
    for(j in 1:length(unique(spec))) {
      mod = ad.test(data$Velocity[which(data$Species==spec[i])], data$Velocity[which(data$Species==spec[j])])
    mod.output = as_tibble(mod$ad)
    dat1 <- data.frame(species1 = spec[i], species2 = spec[j], AD = mod.output$T.AD[1], P = mod.output$` asympt. P-value`[1])
    dat %<>% bind_rows(dat1) 
  }
}

# add Bonferroni correction

dat %<>% filter(species1!=species2) %>% distinct() 
adjusted.p = p.adjust(dat$P, method = "bonferroni")
dat$adj.P = adjusted.p

Table_S7 = dat %>%
  filter(adj.P<0.05)
  
## 14. Modelling of fish diversity effects on herbivory rates across flow speed gradient (Figure)

# Step 1: read in data

## -- call to feeding and flow data -- ##
# point to data locale
data_locale <- "data_intermediate/feeding_and_flow/"

# load data
load(paste0(data_locale, 
            "combined_feeding_and_flow.no_summary.rda"))
            
# Step 2: randomize species identity and obtain the aaverage bite value for each flow interval. First have to intervalize flow speeds and add 0-values.

# set cut off for flow data
flow_cut <- 400

intervalized <- combined_feeding_and_flow.no_summary %>%
  			mutate(Species = fct_recode(Species, Parrotfish = "Bluechin parrotfish", 					Parrotfish = "Bicolored parrotfish")) %>%
			mutate(int_onshore_offshore = round(`onshore-offshore`, digits = 0),
			int_vertical = round(vertical, digits = 0),
			int_along_shore = round(`along-shore`, digits = 0),
			int_3D = round(flow3D, digits = 0)) %>%
			gather(int_onshore_offshore, int_vertical, int_along_shore, int_3D, key = 					"Dimension", value = "Velocity") %>%
			filter(abs(Velocity) <= flow_cut) %>%
			mutate(Velocity = abs(Velocity)) %>%
			group_by(Dimension, Velocity, Species) %>%
			summarise(bites = length(Species)) %>%
			ungroup() %>%
			complete(Species, Velocity, Dimension, fill = list(bites = 0)) %>%
			dplyr::filter(!is.na(Velocity)) %>%
			filter(Species != "Sabertoothed blenny",
         	Species != "Bullseye puffer",
         	Species != "Mexican hogfish",
         	!is.na(Species))			
						

randomized <- intervalized %>%
         mutate(spec_codes = fct_recode(Species, "1" = "Chameleon wrasse", "2" = "Fanged blenny", "3" = "Giant damselfish", "4" = "Rainbow wrasse", "5" = "Ringtail damselfish", "6" = "Spinster wrasse", "7" = "Surgeonfish", "8" = "Yellowtail damselfish", "9" = "Parrotfish", "10" = "King angelfish")) %>%
         filter(Dimension == "int_onshore_offshore")


iterations = 100
velocities <- (unique(randomized$Velocity))
species = 5
store <- vector("list", length = length(velocities))
store2 <- vector("list", length = iterations)


## Single species


for(i in 1:iterations) {
	code <- round(runif(1, min = 1, max = 10), digits = 0)
	data <- randomized %>%
	filter(spec_codes == code)
for(j in 1:length(velocities)) {
		dat <- data %>%
		filter(Velocity == velocities[j]) %>%
		dplyr::select(Velocity, bites)
		store[[j]] = dat
		output <- bind_rows(store)
		store2[[i]] = bind_rows(output)
		}
	
}

one_results <- bind_rows(store2)

## Two species

for(i in 1:iterations) {
	code1 <- round(runif(1, min = 1, max = 10), digits = 0)
	code2 <- round(runif(1, min = 1, max = 10), digits = 0)
	codes <- c(code1, code2)
	data <- randomized %>%
	filter(spec_codes %in% codes) %>%
	group_by(Velocity) %>%
	summarise(bites = sum(bites))
for(j in 1:length(velocities)) {
		dat <- data %>%
		filter(Velocity == velocities[j]) %>%
		dplyr::select(Velocity, bites)
		store[[j]] = dat
		output <- bind_rows(store)
		store2[[i]] = bind_rows(output)
		}
	
}
two_results <- bind_rows(store2)

  ## three species 
  
for(i in 1:iterations) {
	code1 <- round(runif(1, min = 1, max = 10), digits = 0)
	code2 <- round(runif(1, min = 1, max = 10), digits = 0)
	code3 <- round(runif(1, min = 1, max = 10), digits = 0)
	codes <- c(code1, code2, code3)
	data <- randomized %>%
	filter(spec_codes %in% codes) %>%
	group_by(Velocity) %>%
	summarise(bites = sum(bites))
for(j in 1:length(velocities)) {
		dat <- data %>%
		filter(Velocity == velocities[j]) %>%
		dplyr::select(Velocity, bites)
		store[[j]] = dat
		output <- bind_rows(store)
		store2[[i]] = bind_rows(output)
		}
	
}

three_results <- bind_rows(store2)


  ## four species 
  
for(i in 1:iterations) {
	code1 <- round(runif(1, min = 1, max = 10), digits = 0)
	code2 <- round(runif(1, min = 1, max = 10), digits = 0)
	code3 <- round(runif(1, min = 1, max = 10), digits = 0)
	code4 <- round(runif(1, min = 1, max = 10), digits = 0)
	codes <- c(code1, code2, code3, code4)
	data <- randomized %>%
	filter(spec_codes %in% codes) %>%
	group_by(Velocity) %>%
	summarise(bites = sum(bites))
for(j in 1:length(velocities)) {
		dat <- data %>%
		filter(Velocity == velocities[j]) %>%
		dplyr::select(Velocity, bites)
		store[[j]] = dat
		output <- bind_rows(store)
		store2[[i]] = bind_rows(output)
		}
	
}

four_results <- bind_rows(store2)



  ## five species 
  
for(i in 1:iterations) {
	code1 <- round(runif(1, min = 1, max = 10), digits = 0)
	code2 <- round(runif(1, min = 1, max = 10), digits = 0)
	code3 <- round(runif(1, min = 1, max = 10), digits = 0)
	code4 <- round(runif(1, min = 1, max = 10), digits = 0)
	code5 <- round(runif(1, min = 1, max = 10), digits = 0)
	codes <- c(code1, code2, code3, code4, code5)
	data <- randomized %>%
	filter(spec_codes %in% codes) %>%
	group_by(Velocity) %>%
	summarise(bites = sum(bites))
for(j in 1:length(velocities)) {
		dat <- data %>%
		filter(Velocity == velocities[j]) %>%
		dplyr::select(Velocity, bites)
		store[[j]] = dat
		output <- bind_rows(store)
		store2[[i]] = bind_rows(output)
		}
	
}


five_results <- bind_rows(store2)



## Step 3: analyze bootstrapped data

one_results <- one_results %>%
mutate(species = "1")

two_results <- two_results %>%
mutate(species = "2") 

three_results <- three_results %>%
mutate(species = "3")

four_results <- four_results %>%
mutate(species = "4") 

five_results <- five_results %>%
mutate(species = "5") 

## load all observed flow times to get relative number of observations at each flow speed

load("data_intermediate/flow/flow_times_with_feeding_obs.rda")

flow_weights <- flow_times_with_feeding_obs %>%
		  mutate(int_onshore_offshore = round(`onshore-offshore`, digits = 0)) %>%
		  rename(Velocity = int_onshore_offshore) %>%
		  filter(abs(Velocity) <= flow_cut) %>%
		  mutate(Velocity = abs(Velocity),
		  total_obs = length(Site)) %>%
		  group_by(Velocity) %>%
		  summarise(vel.weight = length(Velocity)) %>%
		  mutate(rel_vel_weight = vel.weight/5964)

Figure_5B <- bind_rows(one_results, two_results, three_results, four_results, five_results) %>%
left_join(flow_weights, by = "Velocity") %>%
filter(Velocity <= 200) %>%
mutate(weighted_bites = if_else(Velocity <= 130, bites/vel.weight, 0)) %>%
group_by(species, Velocity) %>%
summarise(mean_bites = mean(weighted_bites),
		  sd_bites = sd(weighted_bites)) %>%
mutate(upper = mean_bites+sd_bites,
	   lower = mean_bites - sd_bites) %>%
ggplot(aes(Velocity, mean_bites)) +
geom_point(aes(color = species)) +
geom_smooth(aes(color = species), se = F) +
scale_color_viridis(discrete = T, option = "inferno", begin = 0, end = .8) +
theme_classic() +
scale_x_continuous(limits = c(0, 160)) +
labs(x = "Flow speed (cm/s)",
	 y = "Bites/observation")


## 13. 

## Regression of flow speeds against urchin and fish effects in 2-month exclusion and 5-day ulva sandwich experiments Figure 6A ##

load("data_intermediate/plate_biomass.rda")

## Use ulva sandwich data for Palmas exposed due to positive effects of consumers there

ulva <- read_excel("data_raw/experiments/ulva_sanduches/rapid_algal_assays_08_23_17_complete.xlsx") %>%
  ungroup() %>%
  #filter(Site == "Las Palmas") %>%
  mutate(mass_lost = Before_weight_total-After_weight_ulva-After_weight_mesh, perc_mass_lost = mass_lost/(Before_weight_total-After_weight_mesh)*100, perc_mass_rem = 100-mass_lost/(Before_weight_total-After_weight_mesh)*100) %>%
  group_by(Treatment, Site, Location) %>%
  summarise(mean_macro = mean(perc_mass_rem, na.rm =T)) %>%
  ungroup() %>%
  mutate(Treatment = fct_recode(Treatment, Open = "Control", Cage = "Exclusion")) %>%
  mutate(Site = fct_recode(Site, camano = "Camaño", palmas = "Las Palmas")) %>%
  mutate(Exposure = tolower(Location)) %>%
  mutate(Site2 = paste(Site, Exposure, sep = "_")) %>%
  dplyr::select(mean_macro, Site2, Treatment) %>%
  rename(Site = Site2) %>%
  spread(Treatment, mean_macro) %>%
  gather(Fence, Open, Roof, key = "Treatment", value = "biomass") %>%
  mutate(rel.effect = (Cage - biomass)/Cage)


#bind flow data from date of experiment


ulva_flows <- all_flows %>%
  filter(date > "2017-08-21",
         date < "2017-09-01") %>%
  group_by(Site) %>%
  summarise(mean_flow = mean(flow3D, na.rm = T))

# combine experiment results and flow

rel.ulva <- ulva %>%
  left_join(ulva_flows) %>%
  dplyr::select(rel.effect, mean_flow, Treatment)


# select variables for plotting

sites = c("Palmas", "Palmas","Camano", "Camano")
seasons = c("warm", "cool", "warm", "cool")
start_date = c("2016-11-16", "2017-06-15","2016-11-17","2017-06-14")
end_date = c ("2017-01-17","2017-08-23","2017-01-18","2017-08-25")
dates <- tibble(sites, seasons, start_date, end_date) %>%
  mutate(elapsed_time = int_length(interval(start_date, end_date))/60/60/24)

## 2. Exploratory plotting and getting to know the data

# plot

biomass_mod <- biomass %>%
  filter(`Biomass category` %in% c("Filamentous", "Foliose", "Thick_foliose", "Branch")) %>%
  group_by(Site, Season, Exposure, Treatment, Plot, Plate) %>%
  summarise(macroalgae=sum(Biomass, na.rm =T)) %>%
  ungroup() %>%
  group_by(Treatment, Site, Season, Exposure) %>%
  summarise(mean_macro = mean(macroalgae, na.rm =T)) %>%
  ungroup() %>%
  mutate(Site2 = paste(tolower(Site), tolower(Exposure), sep = "_"),
         Season = tolower(Season)) %>%
  dplyr::select(-Site) %>%
  rename(Site = Site2) %>%
  spread(Treatment, mean_macro) %>%
  gather(Fence, Open, Roof, key = "Treatment", value = "biomass") %>%
  mutate(rel.effect = (Cage - biomass)/Cage) %>%
  filter(rel.effect > 0) %>%
  #mutate(rel.effect = if_else(rel.effect <= 0, 0, rel.effect)) %>%
  left_join(summ_flows) %>%
  filter(variable == "mean_flow") %>%
  dplyr::select(rel.effect, value, Treatment) %>%
  rename(mean_flow = value) %>%
  bind_rows(rel.ulva) %>%
  mutate(Treatment = fct_recode(Treatment, Fish = "Fence", Urchins = "Roof", `All consumers` = "Open")) %>%
  filter(Treatment != "All consumers") 

Figure_6A <- biomass_mod %>%
  ggplot(aes(mean_flow, rel.effect)) +
  geom_smooth(se = F, color = "black", size = 2, aes(linetype = Treatment)) +
  geom_point(aes(pch = Treatment), size = 2.5) +
  theme_bw() +
  labs(x = "Mean flow speed (cm/s)",
       y = "Proportion of algal biomass removed") +
  scale_linetype_manual(values = c("dashed", "solid")) +
  theme(legend.title = element_blank(),
        legend.position = c(.875, .85),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.size = unit(4,"line"),
        legend.background = element_rect(fill = alpha("white", alpha = 0)))
