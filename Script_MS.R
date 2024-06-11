#--------------------------------------------------------------------------------------------------------------
#                                 META DATA
#------------------------------------------------------------------------------------------------------------
# Script by: Ellinor Jakobsson, 2023-2024
# Script for plotting all figures and doing all analyses in manuscript:
#"Effects of changing snow- and ice cover conditions on phytoplankton biomass and community composition in Lake Erken"

#--------------------------------------------------------------------------------------------------------------
#                                 IMPORT DATA
#------------------------------------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(forestmangr)
library(purrr)
library(broom)
library(polynom)
library(tibble)
library(data.table)
library(vegan)
library(janitor)
library(rstudioapi)

#Speify if any conflicting functions
conflicted::conflicts_prefer(dplyr::count)
conflicted::conflicts_prefer(lubridate::yday)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::mutate)
conflicted::conflicts_prefer(dplyr::rename)
conflicted::conflicts_prefer(dplyr::summarise)
conflicted::conflicts_prefer(vegan::scores)
conflicted::conflicts_prefer(dplyr::arrange)

#Set working directory where all the files are found. Preferably a copy of the folder "Edited_files_for_use" found on github.
dir <- paste0(dirname(getActiveDocumentContext()$path), "/Edited_files_for_use/Files used/")
setwd(dir)
getwd()

#Dir for writing figures and tables
fig_dir <- paste0(dirname(getActiveDocumentContext()$path), "/Results/Figures/")
table_dir <- paste0(dirname(getActiveDocumentContext()$path), "/Results/Tables/")

#List all files in the directory
list.files()

#Load all data files
#Snow on ice data between 1997-2019
Snow_on_ice <- read_excel(paste(dir, "Snow_data_on_ice.xlsx", sep = ""))

#Ice thickness data 1999-2019
Ice_thickness <- read_excel(paste(dir, "Ice_thickness.xlsx", sep = ""))

#Ice cover period data 1997-2019
Ice_cover_period <- read_excel(paste(dir, "/Ice_cover_period.xlsx", sep = ""))

#Phytoplankton data below ice 1997-2019
Phytoplankton_below_ice <- read_excel(paste(dir, "/Phytoplankton_below_ice.xlsx", sep = ""))

#Phytoplankton biomass after ice-out 1997-2019
Phytoplankton_after_ice_out_mean <- read_excel(paste0(dir, "/Phytoplankton_after_ice_out.xlsx"))

#Chla data after ice-out 1997-2019
Chla_after_ice_out <- read_excel(paste0(dir, "/Chla_after_ice_out.xlsx"))

#Write chla within ice-period
Chla_below_ice <- read_excel(paste0(dir, "/Chla_below_ice.xlsx"))

#Write summary data
Big_data_summary <- read_excel(paste0(dir, "/Big_data_summary.xlsx")) 

#Precipitation data from all stations
#Data from Svanberga precip and Vallnora precip
precip_data <- read_excel(paste0(dir, "/precip_data.xlsx"))

#Data from Svanberga precip and Norrveda precip
precip_data2 <- read_excel(paste0(dir, "/precip_data2.xlsx")) 

#Useful functions and defined parameters for the script
#Read function to replace NA with 0. Used later on.
hybrd.rplc_if <- function(x) { mutate_if(x, is.numeric, ~replace(., is.na(.), 0)) }

#Set within which dates to filter i.e., ice on and ice off dates from ice cover period data
Lower_range <- as.Date(Ice_cover_period$`Beginning of ice cover`)
Upper_range <- as.Date(Ice_cover_period$`Ice break up`)

#Format to date format
Snow_on_ice$Date <- as.Date(Snow_on_ice$Date)
#Calculate monthly mean of snow and chla per year
Mean_snow <- Snow_on_ice %>% select(-Datum, -Day) %>% group_by(Winter_year, Month) %>% drop_na(Snow_depth) %>%
  reframe(Mean_snow = mean(Snow_depth))
Mean_chla <- Chla_below_ice %>% select(-Date, -Day) %>% group_by(Winter_year, Month) %>%
  drop_na(`Chl a_µg/l`) %>% reframe(Mean_chla = mean(`Chl a_µg/l`))
Mean_ice <- Ice_thickness %>% group_by(Winter_year, Month) %>% drop_na(Ice_thickness) %>%
  reframe(Mean_ice = mean(Ice_thickness))

#--------------------------------------------------------------------------------------------------------------
#                                 FIGURE 1: PHYTOPLANKTON GROWTH (SNOW AND ICE)
#------------------------------------------------------------------------------------------------------------
#Join mean snow and mean chla
Snow_chla_data <- left_join(Mean_chla, Mean_snow, by = c("Winter_year", "Month"))
Snow_chla_data$Month <- as.numeric(Snow_chla_data$Month)
Snow_chla_data$Month_abbrev <- month.abb[Snow_chla_data$Month]
Snow_chla_data$Month_abbrev <- factor(Snow_chla_data$Month_abbrev, levels = c("Dec", "Jan", "Feb", "Mar", "Apr"))
Snow_chla_data$Winter_year <- as.numeric(Snow_chla_data$Winter_year)
#------------------------------------------------------Plot snow and chl-a
#Jan
data = filter(Snow_chla_data, Month_abbrev == "Jan")
x <- data$Mean_snow
y <- data$Mean_chla
data <- data.frame(x, y) %>% drop_na()
m1 <- lm(y~x, data=data)
#m2 <- nls(y ~ SSasymp(x, Asym, R0, lrc))
m3 <- lm(log(y)~x, data=data)
AIC(m1,m3) 
BIC(m1,m3)
shapiro.test(m1$residuals) #Keep model 1, aka. linear model (non-log)

data = filter(Snow_chla_data, Month_abbrev == "Jan")
m_1 <- m1
p_val <- round(summary(m_1)$coefficients[,4][2], 3)
rsqrd <- round(summary(m_1)$r.squared, 2)
coeff <- as.data.frame(cbind(p_val, rsqrd))
p20 <- ggplot(data, aes(Mean_snow, y = Mean_chla)) +
  geom_smooth(mapping = aes(x = Mean_snow, y = Mean_chla), method = "lm", linetype = "dashed", se = F, col = "black", linewidth = 0.3) + 
  geom_point(mapping = aes(x = Mean_snow, y = Mean_chla), size = 2) + 
  lims(x = c(0,55), y=c(0,NA)) + theme_bw() + ylab("January")  

#Feb
data = filter(Snow_chla_data, Month_abbrev == "Feb")
x <- data$Mean_snow
y <- data$Mean_chla
data <- data.frame(x, y) %>% drop_na()
m2 <- lm(y~x, data=data)
#m3 <- nls(y ~ SSasymp(x, Asym, R0, lrc))
m1 <- lm(log(y)~x, data=data)
AIC(m1,m2)
BIC(m1,m2)
shapiro.test(m1$residuals) #Keep model 1, aka. linear model (non-log)
m_2 <- m1

#Test if model is significant
summary(m_2)
data = filter(Snow_chla_data, Month_abbrev == "Feb")
#Plot the data
p <- ggplot(data, aes(Mean_snow, Mean_chla))
p <- p + geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), se = F, color = "black", linewidth = 0.3, linetype = "dashed") + 
  geom_point(size = 2)
p21 <- p + lims(x=c(0,55)) + scale_y_continuous(trans = "log10", name = "February") +
  annotation_logticks(sides = "l") + theme_bw()

#Mar
data = filter(Snow_chla_data, Month_abbrev == "Mar")
#Find function coefficients for exponential decay
x <- data$Mean_snow
y <- data$Mean_chla
data <- data.frame(x, y) %>% drop_na()
#Test goodness of fit
m1 <- lm(y~x, data=data)
m2 <- lm(log(y)~x, data=data)
#m2 <- nls(y ~ SSasymp(x, Asym, R0, lrc))
AIC(m1,m2)
BIC(m1,m2)
shapiro.test(m1$residuals) #Keep model 2, linear model (log)
m_3 <- m2
summary(m_3)

data = filter(Snow_chla_data, Month_abbrev == "Mar")
#Plot everything
p <- ggplot(data, aes(Mean_snow, Mean_chla))
p <- p + geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), se = F, color = "black", linewidth = 0.3, linetype = "dashed") +
  geom_point(size = 2)
p22 <- p + #geom_function(fun = function(x) (24.905*exp(-0.07*x)), colour = "black", linetype=2) +  
  scale_y_continuous(trans = "log10", name = "March") + #,
                     #sec.axis = sec_axis(trans=~.*1)) + 
                     theme_bw() +
  annotation_logticks(sides = "l") + lims(x=c(0,55))

#---------------------------------------------------------------Join plots and annotate axes
# Functions for facet wrap theme
element_textbox <- function(...) {
  el <- element_text(...)
  class(el) <- c("element_textbox", class(el))
  el
}
element_grob.element_textbox <- function(element, ...) {
  text_grob <- NextMethod()
  rect_grob <- element_grob(calc_element("strip.background", theme_bw()))
  
  ggplot2:::absoluteGrob(
    grid::gList(
      element_grob(calc_element("strip.background", theme_bw())),
      text_grob
    ),
    height = grid::grobHeight(text_grob), 
    width = grid::unit(1, "npc"),
  )
}

p20 <- p20 + theme(axis.title = element_blank()) + ggtitle("Jan") + 
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5), size = 8, face = "bold"),
        plot.title.position = "panel") + lims(x = c(0,55))
p21 <- p21 + theme(axis.title = element_blank()) + ggtitle("Feb") + 
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5), size = 8, face = "bold"),
        plot.title.position = "panel") + lims(x = c(0,55))
p22 <- p22 + theme(axis.title = element_blank()) + ggtitle("Mar") + 
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5), size = 8, face = "bold"),
        plot.title.position = "panel") + lims(x = c(0,55))
# Remove space between plots
p21 <- p21 + theme(plot.margin = unit(c(0,0.1,0.2,0.1), 'lines'), axis.title = element_blank()) + lims(x = c(0,55))
p22 <- p22 + theme(plot.margin = unit(c(0,0.1,0.2,0.1), 'lines'), axis.title = element_blank()) + lims(x = c(0,55))
p20 <- p20 + theme(plot.margin = unit(c(0,0.1,0.2,0.1), 'lines'), axis.title = element_blank()) + lims(x = c(0,55))

# Arrange and annotate big figure
p_chla_snow <- ggarrange(p20, p21, p22, ncol = 3, nrow = 1, align = "hv")

p_final_production_snow_depth <- annotate_figure(p_chla_snow, left = text_grob(expression(paste("Chl-"*alpha*" ("*mu*"g "*L^-1*")")), rot = 90, hjust = 0.5),
                                                                          bottom = text_grob("Snow depth (cm)"))
# Write plot into folder
setwd(fig_dir)
pdf("Figure_1.pdf", width = 6, height = 2.5)
p_final_production_snow_depth

setwd(fig_dir)
tiff("Figure_1.tiff", unit = "cm", height = 4, width = 15, res = 600)
p_final_production_snow_depth
dev.off()
setwd(dir)

#--------------------------------------------------------------------------------------------------------------
#                                 FIGURE 2: ICE PHENOLOGY SCATTERPLOT
#------------------------------------------------------------------------------------------------------------
#Calculate mean chla per month and winter
Chla_after_ice_out_mean <- Chla_after_ice_out %>% select(`Chl a_µg/l`, Winter_year, Month) %>%
  group_by(Month, Winter_year) %>%
  reframe(mean_chla = mean(`Chl a_µg/l`))
#Join mean chla with ice cover period data again
Chla_after_ice_out_mean <- left_join(Chla_after_ice_out_mean, Ice_cover_period, by = c("Winter_year"), relationship = "many-to-many")
#Abbreviate months from numbers to name
Chla_after_ice_out_mean$Month_abb <- month.abb[Chla_after_ice_out_mean$Month]
#Order months for plot
Chla_after_ice_out_mean$Month_abb <- factor(Chla_after_ice_out_mean$Month_abb, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
#Make into DOY
Chla_after_ice_out_mean$Ice_off_DOY <- yday(Chla_after_ice_out_mean$`Ice break up`)
Chla_after_ice_out_mean$Ice_on_DOY <- yday(Chla_after_ice_out_mean$`Beginning of ice cover`)
Chla_after_ice_out_mean <- Chla_after_ice_out_mean %>% separate(`Ice break up`, into = c("Thaw_Year", "Thaw_Month", "Thaw_Day"))
Chla_after_ice_out_mean <- Chla_after_ice_out_mean %>% separate(`Beginning of ice cover`, into = c("Freeze_Year", "Freeze_Month", "Freeze_Day"))
Chla_after_ice_out_mean <- Chla_after_ice_out_mean %>% mutate(Ice_on_DOY = case_when(Freeze_Month == "12" | Freeze_Month == "11" ~ Ice_on_DOY-365, T~Ice_on_DOY))
Chla_after_ice_out_mean <- Chla_after_ice_out_mean %>% filter(Month_abb == "Apr" | Month_abb == "May" | Month_abb == "Jun")

#Regression plots on ice on and ice off DOY
p_ice_phenology <- ggplot() + geom_point(Chla_after_ice_out_mean, mapping = aes(x = Ice_off_DOY, y = mean_chla), col = "#660000", size = 2, shape = 16) + 
  geom_point(Chla_after_ice_out_mean, mapping = aes(x = Ice_on_DOY, y = mean_chla), col = "#000099", size = 2, shape = 17) +
  geom_smooth(Chla_after_ice_out_mean, mapping = aes(x = Ice_off_DOY, y = mean_chla), col = "#660000", method = "lm", se = F, alpha = 0.5, linewidth = 0.3, linetype="dashed") + 
  geom_smooth(Chla_after_ice_out_mean, mapping = aes(x = Ice_on_DOY, y = mean_chla), col = "#000099", method = "lm", se = F, alpha = 0.5, linewidth = 0.3, linetype="dashed") + 
  facet_wrap(~Month_abb, scales = "free_y") + theme_bw() + 
  labs(y = expression(paste("Chl-"*alpha*" ("*mu*"g "*L^-1*")"))) + 
  theme(strip.text = element_text(face = "bold")) + 
  labs(x = "Day-of-the-year")

p_Ice_period <- ggplot() + geom_smooth(Chla_after_ice_out_mean, mapping = aes(x = `Length Ice Cover (days)`, y = mean_chla), col = "black",method = "lm", se = F, linewidth = 0.5, alpha = 0.5, linetype="dashed") +
  geom_point(Chla_after_ice_out_mean, mapping = aes(x = `Length Ice Cover (days)`, y = mean_chla), col = "black", size = 2) +
  facet_wrap(~Month_abb, scales = "free_y") + theme_bw() +
  labs(y = expression(paste("Chl-"*alpha*" ("*mu*"g "*L^-1*")"))) + 
  theme(strip.text = element_text(face = "bold"))

p_final_DOY_box <- p_Ice_period/p_ice_phenology + plot_layout(guides = "collect") & theme(legend.position = "top") &
  plot_annotation(tag_levels = 'a',
                  tag_sep = '',
                  tag_prefix = '',
                  tag_suffix = '') &
  theme(plot.tag.position = "topleft",
        plot.tag = element_text(size = 10, hjust = 0, vjust = 0)) &
  theme(plot.tag = element_text(face = 'bold'))

setwd(fig_dir)
pdf("Figure_2.pdf", width=2*2.54, height=2*2.54)
p_final_DOY_box 
dev.off()

setwd(fig_dir)
tiff("Figure_2.tiff", width=12, height=9, unit = "cm", res = 600)
p_final_DOY_box 
dev.off()

#Regressions on ice phenology
# Rsqr values
Ice_period_table_r <- Chla_after_ice_out_mean %>% group_by(Month) %>% lm_table(mean_chla ~ `Length Ice Cover (days)`) %>%
  select(Month, Rsqr_adj) %>% rename("R2" = "Rsqr_adj") %>%
  mutate(R2 = round(R2, 2))
# p values
Ice_period_table_p <- Chla_after_ice_out_mean  %>% nest(data = -c(Month)) %>%
  mutate(model = map(data, ~lm(mean_chla ~ `Length Ice Cover (days)`, data = .)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>% filter(term == "`Length Ice Cover (days)`") %>%
  select(Month, p.value) %>%
  rename("p" = "p.value") %>% mutate(p = round(p, 2))
#---------------------------ICE ON DOY
# Rsqr values
Ice_on_table_r <- Chla_after_ice_out_mean %>% group_by(Month) %>% lm_table(mean_chla ~ Ice_on_DOY) %>%
  select(Month, Rsqr_adj) %>% rename("R2" = "Rsqr_adj") %>%
  mutate(R2 = round(R2, 2))
# p values
Ice_on_table_p <- Chla_after_ice_out_mean  %>% nest(data = -c(Month)) %>%
  mutate(model = map(data, ~lm(mean_chla ~ Ice_on_DOY, data = .)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>% filter(term == "Ice_on_DOY") %>% 
  select(Month, p.value) %>%
  rename("p" = "p.value") %>% mutate(p = round(p, 2))
#---------------------------ICE OFF DOY
# Rsqr values
Ice_off_table_r <- Chla_after_ice_out_mean %>% group_by(Month) %>% lm_table(mean_chla ~ Ice_off_DOY) %>%
  select(Month, Rsqr_adj) %>% rename("R2" = "Rsqr_adj") %>%
  mutate(R2 = round(R2, 2))
# p values
Ice_off_table_p <- Chla_after_ice_out_mean  %>% nest(data = -c(Month)) %>%
  mutate(model = map(data, ~lm(mean_chla ~ Ice_off_DOY, data = .)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>% filter(term == "Ice_off_DOY") %>% 
  select(Month, p.value) %>%
  rename("p" = "p.value") %>% mutate(p = round(p, 2))
#----------------------------COMBINE TABLES
# combine
Ice_period_table <- left_join(Ice_period_table_r, Ice_period_table_p)
Ice_on_table <- left_join(Ice_on_table_r, Ice_on_table_p)
Ice_off_table <- left_join(Ice_off_table_r, Ice_off_table_p)
Ice_off_table$Explanatory <- "Ice off (DOY)"
Ice_on_table$Explanatory <- "Ice on (DOY)"
Ice_period_table$Explanatory <- "Ice period"
Ice_period_table <- Ice_period_table %>% select(Explanatory, Month, R2, p)
Ice_on_table <- Ice_on_table %>% select(Explanatory, Month, R2, p)
Ice_off_table <- Ice_off_table %>% select(Explanatory, Month, R2, p)
Ice_phenology_lm_table <- rbind(Ice_period_table, Ice_on_table, Ice_off_table)
Ice_phenology_lm_table$Explanatory <- factor(Ice_phenology_lm_table$Explanatory, levels = c("Ice period", "Ice on (DOY)", "Ice off (DOY)"))


#--------------------------------------------------------------------------------------------------------------
#                                 FIGURE 3: TROPHIC MODES
#------------------------------------------------------------------------------------------------------------
#CALCULATE FUNCTIONAL GROUP - JOIN MEAN WITH MEAN
#Calculate sum of biomass per dat and genus
Phyto_data <- Phytoplankton_below_ice %>% group_by(Genus, Winter_year, Month, Day, Functional_group) %>% dplyr::summarise(across(biomass_ug_l, list(sum)))
#Calculate mean biomass per genera and month within the ice cover period
Mean_phytoplankton_genera <- Phyto_data %>% group_by(Genus, Winter_year, Month, Functional_group) %>% 
  reframe(Mean_biomass = mean(biomass_ug_l_1))
#Left join by monthly mean
Phyto_snow <- left_join(Mean_phytoplankton_genera, Mean_snow, by = c("Month", "Winter_year"))
Phyto_snow <- Phyto_snow %>% mutate(Genus = case_when(Genus == "<NA>" ~ "Other", T~Genus))
#Remove NA of snow
Phyto_snow <- Phyto_snow %>% drop_na(Mean_snow)
#Make new column called taxa from genus and remove unknown genera
Phyto_snow$Taxa <- Phyto_snow$Genus
Phyto_snow <- Phyto_snow %>% ungroup() %>% select(-Genus)
Perc <- Phyto_snow %>% filter(Taxa != "Other")
#Calculate sum of biomass per functional group and snow depth as well as total of all three groups
Perc <- Perc %>% select(-Taxa) %>% group_by(Winter_year, Month, Functional_group, Mean_snow) %>%
  mutate(sum_feeding = sum(Mean_biomass)) %>%
  group_by(Winter_year, Month, Mean_snow) %>%
  mutate(Sum_total = sum(Mean_biomass)) %>%
  filter(Sum_total > 0)
#Divide sum biomass per group with sum of total biomass
Perc <- Perc %>% select(-Mean_biomass)%>% group_by(Winter_year, Month, Mean_snow, Functional_group) %>%
  mutate(Perc = sum_feeding/Sum_total)
#Clean up data file by claulating percentage from proportin and remove unused columns
Perc <- add_column(Perc, Percentage = Perc$Perc*100, .before = 1)
Perc <- Perc %>% select(-Perc, -Sum_total, -sum_feeding)
#Claculate mean percange per month
Perc <- Perc %>% group_by(Winter_year, Month, Mean_snow, Functional_group) %>% reframe(Percentage = mean(Percentage))
#Select only months of Jan-Mar
Perc <- Perc %>% filter(Month == 01 | Month == 02 | Month == 03) %>% mutate(Month = as.factor(Month)) %>%
  mutate("Trophic identity" = Functional_group)

#Make a new column with month abbreviated
Perc$Month_abb <- month.abb[Perc$Month]
#Put correct order of groups for plotting
Perc <- Perc %>% mutate(Functional_group = case_when(Functional_group == "Autotroph" ~ "Autotrophs",
                                                                        Functional_group == "Mixotroph" ~ "Mixotrophs",
                                                                        Functional_group == "Heterotroph" ~ "Heterotrophs"))
Perc$`Trophic identity` <- factor(Perc$Functional_group, levels = c("Autotrophs", "Mixotrophs", "Heterotrophs"))
Perc$Month_abb <- factor(Perc$Month_abb, levels = c("Jan", "Feb", "Mar"))
Perc2 <- Perc
p_proportion_mixotroph <- ggplot() + 
  geom_point(Perc, mapping = aes(x = Mean_snow, y = Percentage), size = 2) + 
  theme_bw() + facet_grid(Month_abb~`Trophic identity`) + 
  labs(x = "Snow depth (cm)", y = "(%)", size = 10) + lims(y = c(0, 100)) + 
  theme(strip.text = element_text(face = "bold")) + 
  geom_smooth(Perc, mapping = aes(x = Mean_snow, y = Percentage), se = F, linewidth = 0.3, linetype = "dashed", method = "lm", col = "black") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

#GROWTH

#Calculate sum of biomass per genera and day
Phyto_data <- Phytoplankton_below_ice %>% group_by(Genus, Winter_year, Month, Day, Functional_group) %>% dplyr::summarise(across(biomass_ug_l, list(sum)))
#Join phytoplankton data and snow data
Mean_phytoplankton_genera <- Phyto_data %>% group_by(Genus, Winter_year, Month, Functional_group) %>% 
  reframe(Mean_biomass = mean(biomass_ug_l_1))
#Left join by monthly mean
Phyto_snow <- left_join(Mean_phytoplankton_genera, Mean_snow, by = c("Month", "Winter_year"))
Phyto_snow <- Phyto_snow %>% mutate(Genus = case_when(Genus == "<NA>" ~ "Other", T~Genus))
#Remove NA of snow 
Phyto_snow <- Phyto_snow %>% drop_na(Mean_snow)
#Make new column called taxa from genus and remove unknown genera
Phyto_snow$Taxa <- Phyto_snow$Genus
Phyto_snow <- Phyto_snow %>% ungroup() %>% select(-Genus)
#Remove unknown taxa, choose months of jan-mar and calculate mean biomass per month and snow depth
Perc <- Phyto_snow %>% filter(Taxa != "Other") %>%  filter(Month == 01 | Month == 02 | Month == 03) %>%
  select(-Taxa) %>% ungroup() %>%
  group_by(Winter_year, Functional_group, Month, Mean_snow) %>% reframe(Mean_feeding = mean(Mean_biomass)) %>%
  mutate("Trophic identity" = Functional_group)

#Make a new column with month abbreviated
Perc$Month_abb <- month.abb[Perc$Month]

#Put correct order for plotting
Perc <- Perc %>% mutate(Functional_group = case_when(Functional_group == "Autotroph" ~ "Autotrophs",
                                                     Functional_group == "Mixotroph" ~ "Mixotrophs",
                                                     Functional_group == "Heterotroph" ~ "Heterotrophs"))

Perc$`Trophic identity` <- factor(Perc$Functional_group, levels = c("Autotrophs", "Mixotrophs", "Heterotrophs"))

Perc$Month_abb <- factor(Perc$Month_abb, levels = c("Jan", "Feb", "Mar"))

#Plot absolute biomass in relation to snow depth
p_absolute_mixotroph <- ggplot() + geom_point(Perc, mapping = aes(x = Mean_snow, y = Mean_feeding), size = 2) +
  theme_bw() + facet_grid(Month_abb~`Trophic identity`) + labs(x = "Snow depth (cm)", y = expression(paste("Biomass  ("*mu*"g "*L^-1*")")), size = 10) +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.x.top = element_blank()) +
  theme(strip.text = element_text(face = "bold")) + 
  geom_smooth(Perc, mapping = aes(x = Mean_snow, y = Mean_feeding), se = F, linewidth = 0.3, linetype = "dashed", method = "lm", col = "black")

#Combine plots and write final version
p_mixotrophy <- ggarrange(p_proportion_mixotroph, p_absolute_mixotroph, nrow = 2, align = "v")
p_mixotrophy <- p_mixotrophy + theme(legend.position = "right") +
  plot_layout(guides = "collect")

setwd(fig_dir)
pdf("Figure_3.pdf", height = 4, width=6)
p_mixotrophy
dev.off()

setwd(fig_dir)
tiff("Figure_3.tiff", height = 15, width=12, units = "cm", res = 600)
p_mixotrophy
dev.off()

# #-------------------------------Linear regressions on trophic identities
# #PERCENTAGE

# #Non-log transformed
# Mixo_lm <- lm(Percentage~Mean_snow, data = filter(Perc2, Functional_group == "Mixotroph" & Month == 1))
# shapiro.test(Mixo_lm$residuals)
# p_val <- round(summary(Mixo_lm)$coefficients[,4][2], 3)
# rsqrd <- round(summary(Mixo_lm)$r.squared, 2)
# coeff <- as.data.frame(cbind(p_val, rsqrd))
# 
# #p=0.20, R2=0.082
# Mixo_lm <- lm(Percentage~Mean_snow, data = filter(Perc2, Functional_group == "Mixotroph" & Month == 2))
# shapiro.test(Mixo_lm$residuals)
# p_val <- round(summary(Mixo_lm)$coefficients[,4][2], 3)
# rsqrd <- round(summary(Mixo_lm)$r.squared, 2)
# coeff <- as.data.frame(cbind(p_val, rsqrd))
# 
# #p=0.03, R2=0.30
# Mixo_lm <- lm(Percentage~Mean_snow, data = filter(Perc2, Functional_group == "Mixotroph" & Month == 3))
# shapiro.test(Mixo_lm$residuals)
# p_val <- round(summary(Mixo_lm)$coefficients[,4][2], 3)
# rsqrd <- round(summary(Mixo_lm)$r.squared, 2)
# coeff <- as.data.frame(cbind(p_val, rsqrd))
# #p=0.07, R2=0.20

# #Autotrophs
# Mixo_lm <- lm(Percentage~Mean_snow, data = filter(Perc2, Functional_group == "Autotroph" & Month == 1))
# shapiro.test(Mixo_lm$residuals)
# Mixo_lm <- summary(Mixo_lm)
# Mixo_lm
# #p=0.21, R2=0.07
# Mixo_lm <- lm(Percentage~Mean_snow, data = filter(Perc2, Functional_group == "Autotroph" & Month == 2))
# shapiro.test(Mixo_lm$residuals)
# Mixo_lm <- summary(Mixo_lm)
# Mixo_lm
# #p=0.03, R2=0.30
# Mixo_lm <- lm(Percentage~Mean_snow, data = filter(Perc2, Functional_group == "Autotroph" & Month == 3))
# shapiro.test(Mixo_lm$residuals)
# Mixo_lm <- summary(Mixo_lm)
# Mixo_lm
# #p=0.08, R2=0.18

# #Heterotroph
# Mixo_lm <- lm(log(Percentage+1)~Mean_snow, data = filter(Perc2, Functional_group == "Heterotroph" & Month == 1))
# shapiro.test(Mixo_lm$residuals)
# Mixo_lm <- summary(Mixo_lm)
# Mixo_lm
# #p=0.36, R2=-0.007
# Mixo_lm <- lm(Percentage~Mean_snow, data = filter(Perc2, Functional_group == "Heterotroph" & Month == 2))
# shapiro.test(Mixo_lm$residuals)
# Mixo_lm <- summary(Mixo_lm)
# Mixo_lm
# #p=0.18, R2=0.08
# Mixo_lm <- lm(log(Percentage+1)~Mean_snow, data = filter(Perc2, Functional_group == "Heterotroph" & Month == 3))
# shapiro.test(Mixo_lm$residuals)
# Mixo_lm <- summary(Mixo_lm)
# Mixo_lm
# #p=0.56, R2=-0.06

# #GROWTH

# Perc$Functional_group <- Perc$`Trophic identity`
# Mixo_lm <- lm(log(Mean_feeding+1)~Mean_snow, data = filter(Perc, Functional_group == "Mixotroph" & Month == 1))
# shapiro.test(Mixo_lm$residuals)
# Mixo_lm <- summary(Mixo_lm)
# Mixo_lm
# #p=0.38, R2=-0.02. Log transformed.
# Hetero_lm <- lm(log(Mean_feeding+1)~Mean_snow, data = filter(Perc, Functional_group == "Heterotroph" & Month == 1))
# shapiro.test(Hetero_lm$residuals)
# Hetero_lm <- summary(Hetero_lm)
# Hetero_lm
# #p=0.94, R2=-0.11. Log transformed. 
# Auto_lm <- lm(log(Mean_feeding+1)~Mean_snow, data = filter(Perc, Functional_group == "Autotroph" & Month == 1))
# shapiro.test(Auto_lm$residuals)
# Auto_lm <- summary(Auto_lm)
# Auto_lm
# #p=0.35, R2=-0.003. Log transformed. 

# Perc$Functional_group <- Perc$`Trophic identity`
# Mixo_lm <- lm(log(Mean_feeding+1)~Mean_snow, data = filter(Perc, Functional_group == "Mixotroph" & Month == 2))
# shapiro.test(Mixo_lm$residuals)
# Mixo_lm <- summary(Mixo_lm)
# Mixo_lm
# #Log transformed. Add constant to avoid loosing true zeros in log-transformation
# Hetero_lm <- lm(log(Mean_feeding+1)~Mean_snow, data = filter(Perc, Functional_group == "Heterotroph" & Month == 2))
# shapiro.test(Hetero_lm$residuals)
# Hetero_lm <- summary(Hetero_lm)
# Hetero_lm
# 
# Auto_lm <- lm(Mean_feeding~Mean_snow, data = filter(Perc, Functional_group == "Autotroph" & Month == 2))
# shapiro.test(Auto_lm$residuals)
# Auto_lm <- summary(Auto_lm)
# Auto_lm

# #Log transformed 
# Perc$Functional_group <- Perc$`Trophic identity`
# Mixo_lm <- lm(log(Mean_feeding+1)~Mean_snow, data = filter(Perc, Functional_group == "Mixotroph" & Month == 3))
# shapiro.test(Mixo_lm$residuals)
# Mixo_lm <- summary(Mixo_lm)
# Mixo_lm
# #Log transformed. Add constant to avoid loosing true zeros in log-transformation
# Hetero_lm <- lm(log(Mean_feeding+1)~Mean_snow, data = filter(Perc, Functional_group == "Heterotroph" & Month == 3))
# shapiro.test(Hetero_lm$residuals)
# Hetero_lm <- summary(Hetero_lm)
# Hetero_lm
# #Log transformed
# Auto_lm <- lm(log(Mean_feeding+1)~Mean_snow, data = filter(Perc, Functional_group == "Autotroph" & Month == 3))
# shapiro.test(Auto_lm$residuals)
# Auto_lm <- summary(Auto_lm)
# Auto_lm

#--------------------------------------------------------------------------------------------------------------
#                                 FIGURE 4: AFTER ICE OUT NMDS, DOMINANT TAXA IN APRIL
#------------------------------------------------------------------------------------------------------------
#Calculate dominant taxa in April after ice off
Dominant_taxa <- Phytoplankton_after_ice_out_mean %>% filter(Month_abb == "Apr") %>% group_by(Genus, Month_abb) %>% reframe(Mean_phyto_biomass = mean(mean_biomass))
# NMDS ice phenology after ice-out
#Replace NA with 0
Phytoplankton_after_ice_out_mean <- hybrd.rplc_if(Phytoplankton_after_ice_out_mean)
#Order months
Phytoplankton_after_ice_out_mean$Month_abb <- factor(Phytoplankton_after_ice_out_mean$Month_abb, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
NMDS_data_filter <- Phytoplankton_after_ice_out_mean %>% spread(Genus, mean_biomass)
NMDS_data_filter <- NMDS_data_filter %>% filter(Month_abb == "Apr" | Month_abb == "May" | Month_abb == "Jun")

#Make 0 to NA and add column for month abbreviation
NMDS_data <- NMDS_data_filter
NMDS_data[NMDS_data == 0] <- NA
df <- NMDS_data
df_filter <- df
df_filter <- df_filter %>% select(-Other) #Remove unknown taxa
df_filter <- df_filter[,colSums(!is.na(df_filter)) > 0*nrow(df_filter)] #Filter data based on percentage of non-NA cases, none used
Phyto_snow <- hybrd.rplc_if(df_filter) #Replace NA with 0
Phyto_snow_filter <- Phyto_snow[rowSums(Phyto_snow[,c(14:110)])>0,] #Remove samplings without any data
NMDS_spec <- Phyto_snow_filter[, c(14:110)]
NMDS_env <- Phyto_snow_filter[,c(10:13)]
sp.hel <- decostand(NMDS_spec, method="hellinger")
NMDS_run <- metaMDS(sp.hel, autotransform = F, k=2, trymax = 100)
ef <- envfit(NMDS_run, NMDS_env[1:4], permu = 999)
ef

# PLOTTING
data.scores <- as.data.frame(scores(NMDS_run)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$Ice_period <- NMDS_env$`Length Ice Cover (days)`
data.scores$Ice_off_DOY <- NMDS_env$Ice_off_DOY
data.scores$Ice_on_DOY <- NMDS_env$Ice_on_DOY
data.scores$Month <- NMDS_env$Month_abb
head(data.scores)
species.scores <- as.data.frame(scores(NMDS_run, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores) 
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))
vec.df$variables <- c("Ice-period", "Ice-off", "Ice-on")

p9 <- ggplot() + theme_bw() + geom_hline(yintercept=0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept=0, linetype = "dashed", color = "gray50") +
  geom_point(data=data.scores, mapping = aes(x=NMDS1, y=NMDS2, col = Ice_off_DOY, shape = Month), size=4) +
  scale_colour_gradient2(low = "#FDB863",
                         mid = "#B2ABD2",
                         high = "#542788",
                         midpoint = 115) + 
  labs(shape = NULL, colour = "Ice-off") +
  geom_segment(data = vec.df,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")),
               colour="black",
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = NMDS1, y=NMDS2, label = variables),
            size=3.5) +
  coord_equal() + 
  theme(axis.text.x = element_blank())

p10 <- ggplot() + geom_hline(yintercept=0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept=0, linetype = "dashed", color = "gray50") +
  geom_text(data=species.scores, aes(x=NMDS1, y=NMDS2, label=species), alpha=0.6, size = 2, position = position_dodge(width = .9)) +
  coord_equal() +
  theme_bw() +
  scale_colour_gradient2(low = "#FDB863",
                         mid = "#B2ABD2",
                         high = "#542788",
                         midpoint = 25)

p9 <- p9 +  
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top") 

p10 <- p10 +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top")

p_final_snow_depth_NMDS <- (p9/p10) +
  plot_annotation(tag_levels = list(c("")),
                  tag_sep = '',
                  tag_prefix = '',
                  tag_suffix = '') &
  theme(plot.tag.position = "topleft",
        plot.tag = element_text(size = 8, hjust = 0, vjust = 0)) &
  theme(plot.tag = element_text(face = 'bold'),
        legend.position = "right")

#Write figures as pdf and tiff
setwd(fig_dir)
pdf("Figure_4.pdf", width=2*2.54, height=2*2.54)
p_final_snow_depth_NMDS 
dev.off()

setwd(fig_dir)
tiff("Figure_4.tiff", width=14, height=16, unit = "cm", res = 600)
p_final_snow_depth_NMDS 
dev.off()

# PERMANOVA
Permanova_env <- NMDS_env %>% mutate(Month = as.numeric(Month_abb))
sp.hel <- decostand(NMDS_spec, method="hellinger")
Permanova_ice_phenology <- adonis2(sp.hel ~ Ice_on_DOY*Ice_off_DOY*Month, data=Permanova_env, permutations=99)
Explanatory <- rownames(Permanova_ice_phenology)
Permanova_ice_phenology <- as.data.frame(Permanova_ice_phenology[1:5]) %>% mutate(SumOfSqs = round(SumOfSqs, 3),
                                                                                  R2 = round(R2, 3),
                                                                                  F = round(F, 3),
                                                                                  `Pr(>F)` = round(`Pr(>F)`, 3))


#--------------------------------------------------------------------------------------------------------------
#                                 SUPPLEMENTARY FIGURE 1: PRECIPITATION
#------------------------------------------------------------------------------------------------------------
#Plot precipitation between the meterological sites
precip_plot <- ggplot() + geom_point(precip_data, mapping = aes(x = Rain.y, y = Rain.x), alpha = 0.5) + labs(y = "", x = "Vallnora prec. (mm)") + 
  theme_bw() + geom_abline(mapping = aes(intercept = 0, slope = 1), linetype = "dashed") + 
  stat_cor(precip_data, mapping = aes(x = Rain.y, y = Rain.x), cor.coef.name = "r", method="pearson") + theme(axis.text.y = element_blank())

precip_plot2 <- ggplot() + geom_point(precip_data2, mapping = aes(x = Rain.y, y = Rain.x), alpha = 0.5) + labs(y = "Svanberga prec. (mm)", x = "Norrveda prec. (mm)") + 
  theme_bw() + geom_abline(mapping = aes(intercept = 0, slope = 1), linetype = "dashed") + 
  stat_cor(precip_data2, mapping = aes(x = Rain.y, y = Rain.x), cor.coef.name = "r", method="pearson")
#Final plot
p_final_precip <- (precip_plot2|precip_plot) & plot_annotation(tag_levels = 'a') &
  theme(plot.tag.position = "topleft",
        plot.tag = element_text(size = 10, hjust = 0, vjust = 0)) &
  theme(plot.tag = element_text(face = 'bold')) + 
  theme(plot.margin = unit(c(0.1,0.9,0.1,0.1), 'lines'))

p_final_precip

# Write plot into folder
setwd(fig_dir)
pdf("Figure_S1.pdf", width = 5, height = 4)
p_final_precip
dev.off()

setwd(fig_dir)
tiff("Figure_S1.tiff", unit = "cm", height = 7, width = 15, res = 600)
p_final_precip
dev.off()
setwd(dir)

#--------------------------------------------------------------------------------------------------------------
#                                 SUPPLEMENTARY FIGURE 2: ICE AND SNOW THICKNESS RELATION
#------------------------------------------------------------------------------------------------------------
#Calculate mean ice and snow thickness
#Join snow and ice thickness
Snow_ice_data <- left_join(Mean_snow, Mean_ice, by = c("Month", "Winter_year")) %>% drop_na()

Snow_ice_data <- Snow_ice_data %>% group_by(Winter_year) %>% reframe(Mean_snow = mean(Mean_snow),
                                                                     Mean_ice = mean(Mean_ice))

p_snow_ice_cor <- ggplot(Snow_ice_data, aes(Mean_ice, y=Mean_snow)) + 
  geom_smooth(method = "lm", se = F, col = "black", linetype = "dashed") +
  geom_point(size = 3) + theme_bw() + 
  labs(x = "Ice thickness (cm)", y = "Snow thickness (cm)") + 
  lims(x = c(0,NA),
       y = c(0,NA))

# Write plot into folder
setwd(fig_dir)
pdf("Figure_S2.pdf", width = 6, height = 2.5)
p_snow_ice_cor
dev.off()

setwd(fig_dir)
tiff("Figure_S2.tiff", unit = "cm", height = 4, width = 15, res = 600)
p_snow_ice_cor
dev.off()
setwd(dir)
#--------------------------------------------------------------------------------------------------------------
#                                 SUPPLEMENTARY FIGURE 3: ICE AND CHLA
#------------------------------------------------------------------------------------------------------------
#ICE THICKNESS AND CHLA
Ice_chla_data <- left_join(Mean_chla, Mean_ice, by = c("Winter_year", "Month")) %>% drop_na()
Ice_chla_data$Month <- as.numeric(Ice_chla_data$Month)
Ice_chla_data$Month_abbrev <- month.abb[Ice_chla_data$Month]
Ice_chla_data$Month_abbrev <- factor(Ice_chla_data$Month_abbrev, levels = c("Dec", "Jan", "Feb", "Mar", "Apr"))
Ice_chla_data$Winter_year <- as.numeric(Ice_chla_data$Winter_year)

#Plotting
#Jan
data = filter(Ice_chla_data, Month_abbrev == "Jan")
x <- data$Mean_ice
y <- data$Mean_chla
data <- data.frame(x, y) %>% drop_na()
m1 <- lm(y~x, data=data)
#m2 <- nls(y ~ SSasymp(x, Asym, R0, lrc))
m3 <- lm(log(y)~x, data=data)
AIC(m1,m3) 
BIC(m1,m3)
shapiro.test(m1$residuals) #Keep model 1, aka. linear model (non-log)

data = filter(Ice_chla_data, Month_abbrev == "Jan")
m_1 <- m1
summary(m_1)
p23 <- ggplot(data, aes(Mean_ice, y = Mean_chla)) +
  #geom_smooth(mapping = aes(x = Mean_ice, y = Mean_chla), method = "lm", linetype = "dashed", se = F, col = "black", linewidth = 0.3) + 
  geom_point(mapping = aes(x = Mean_ice, y = Mean_chla), size = 2) + 
  lims(x = c(0,55), y=c(0,NA)) + theme_bw() + ylab("January") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), se = F, color = "black", linewidth = 0.3, linetype = "dashed") 

#Feb
data = filter(Ice_chla_data, Month_abbrev == "Feb")
x <- data$Mean_ice
y <- data$Mean_chla
data <- data.frame(x, y) %>% drop_na()
m2 <- lm(y~x, data=data)
#m3 <- nls(y ~ SSasymp(x, Asym, R0, lrc))
m1 <- lm(log(y)~x, data=data)
AIC(m1,m2)
BIC(m1,m2)
shapiro.test(m2$residuals) #Keep model 1, aka. linear model (non-log)
m_2 <- m2
#Test if model is significant
summary(m_2)
data = filter(Ice_chla_data, Month_abbrev == "Feb")
#Plot the data
p <- ggplot(data, aes(Mean_ice, Mean_chla))
p <- p + #geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), se = F, color = "black", linewidth = 0.3, linetype = "dashed") + 
  geom_point(size = 2)
p24 <- p + lims(x=c(0,55)) + scale_y_continuous(name = "February") + theme_bw() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), se = F, color = "black", linewidth = 0.3, linetype = "dashed")

#Mar
data = filter(Ice_chla_data, Month_abbrev == "Mar")
#Find function coefficients for exponential decay
x <- data$Mean_ice
y <- data$Mean_chla
data <- data.frame(x, y) %>% drop_na()
#Test goodness of fit
m1 <- lm(y~x, data=data)
m2 <- lm(log(y)~x, data=data)
#m2 <- nls(y ~ SSasymp(x, Asym, R0, lrc))
AIC(m1,m2)
BIC(m1,m2)
shapiro.test(m1$residuals) #Keep model 2, linear model (log)
m_3 <- m2
summary(m_3)

data = filter(Ice_chla_data, Month_abbrev == "Mar")
#Plot everything
p <- ggplot(data, aes(Mean_ice, Mean_chla))
p <- p + geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), se = F, color = "black", linewidth = 0.3, linetype = "dashed") +
  geom_point(size = 2)
p25 <- p + #geom_function(fun = function(x) (24.905*exp(-0.07*x)), colour = "black", linetype=2) +  
  scale_y_continuous(trans = "log10", name = "March") + #,
  #sec.axis = sec_axis(trans=~.*1)) + 
  theme_bw() +
  annotation_logticks(sides = "l") + lims(x=c(0,55))


p23 <- p23 + theme(axis.title = element_blank()) + ggtitle("Jan") + 
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5), size = 8, face = "bold"),
        plot.title.position = "panel") + lims(x = c(0,55))
p24 <- p24 + theme(axis.title = element_blank()) + ggtitle("Feb") + 
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5), size = 8, face = "bold"),
        plot.title.position = "panel") + lims(x = c(0,55))
p25 <- p25 + theme(axis.title = element_blank()) + ggtitle("Mar") + 
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5), size = 8, face = "bold"),
        plot.title.position = "panel") + lims(x = c(0,55))

p23 <- p23 + theme(plot.margin = unit(c(0,0.1,0.2,0.1), 'lines'), axis.title = element_blank()) + lims(x = c(0,55))
p24 <- p24 + theme(plot.margin = unit(c(0,0.1,0.2,0.1), 'lines'), axis.title = element_blank()) + lims(x = c(0,55))
p25 <- p25 + theme(plot.margin = unit(c(0,0.1,0.2,0.1), 'lines'), axis.title = element_blank()) + lims(x = c(0,55))

p_chla_ice <- ggarrange(p23, p24, p25, ncol = 3, nrow = 1, align = "hv")
p_final_production_ice_thickness <- annotate_figure(p_chla_ice, left = text_grob(expression(paste("Chl-"*alpha*" ("*mu*"g "*L^-1*")")), rot = 90, hjust = 0.5),
                                                    bottom = text_grob("Ice thickness (cm)"))
# Write plot into folder
setwd(fig_dir)
pdf("Figure_S3.pdf", width = 6, height = 2.5)
p_final_production_ice_thickness
dev.off()

setwd(fig_dir)
tiff("Figure_S3.tiff", unit = "cm", height = 4, width = 15, res = 600)
p_final_production_ice_thickness
dev.off()
setwd(dir)

#--------------------------------------------------------------------------------------------------------------
#                                 SUPPLEMENTARY FIGURE 4: ZOOPLANKTON + NUTRIENTS
#------------------------------------------------------------------------------------------------------------
#Plot all variables
p_zoops <- ggplot(filter(Big_data_summary, Taxa == "Zooplankton"), mapping = aes(x = Month_abbrev, y = Biomass)) +
  geom_boxplot(fill = "white", outlier.shape = NA) +
  geom_point(size = 1) + 
  theme_bw() + theme(strip.text = element_text(face = "bold"),
                     axis.text = element_text(size = 8),
                     axis.title = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) +
  facet_wrap(~Taxa, scales = "free_y") + labs(x = "", y = expression("Biovolume (mm"^3*" L"^-1*")")) + 
  scale_x_discrete(labels=c("J", "F", "M", "A", "M", "J"))

p_phyto <- ggplot(filter(Big_data_summary, Taxa == "Phytoplankton"), mapping = aes(x = Month_abbrev, y = Biomass)) +
  geom_boxplot(fill = "white", outlier.shape = NA) +
  geom_point(size = 1) + 
  theme_bw() + theme(strip.text = element_text(face = "bold"),
                     axis.text = element_text(size = 8),
                     axis.title = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) +
  facet_wrap(~Taxa, scales = "free_y") + labs(x = "", y = expression("Chl-"*alpha*" ("*mu*"g L"^-1*")")) + 
  scale_x_discrete(labels=c("J", "F", "M", "A", "M", "J"))

p_nutrients <- 
  ggplot() +
  geom_boxplot(filter(Big_data_summary, Chemistry == "TN" | Chemistry == "TP" | Chemistry == "Si"), 
               mapping = aes(x = Month_abbrev, y = Concentration),
               fill = "white", outlier.shape = NA) + 
  theme_bw() + theme(strip.text = element_text(face = "bold"),
                     axis.text = element_text(size = 8),
                     axis.title = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) + geom_point(filter(Big_data_summary, Chemistry == "TN" | Chemistry == "TP" |
                                                                                  Chemistry == "Si"), 
                                                                         mapping = aes(x = Month_abbrev, y = Concentration), size = 1) +
  facet_wrap(~Chemistry, scales = "free_y") + 
  labs(x="", y = expression(paste("mg "*L^-1))) + 
  scale_x_discrete(labels=c("J", "F", "M", "A", "M", "J"))

p_temp <- ggplot() +
  geom_boxplot(filter(Big_data_summary, Chemistry == "Temperature"), 
               mapping = aes(x = Month_abbrev, y = Concentration), fill = "white", outlier.shape = NA) +
  geom_point(filter(Big_data_summary, Chemistry == "Temperature"), 
             mapping = aes(x = Month_abbrev, y = Concentration), size = 1) + 
  theme_bw() + theme(strip.text = element_text(face = "bold"),
                     axis.text = element_text(size = 8),
                     axis.title = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) +
  facet_wrap(~Chemistry, scales = "free_y") + labs(y = expression(""*~degree*C*""), x = "") + 
  scale_x_discrete(labels=c("J", "F", "M", "A", "M", "J"))

p_light <- ggplot() +
  geom_boxplot(filter(Big_data_summary, Chemistry == "SW radiation"), 
               mapping = aes(x = Month_abbrev, y = Concentration), fill = "white", outlier.shape = NA) +
  geom_point(filter(Big_data_summary, Chemistry == "SW radiation"), 
             mapping = aes(x = Month_abbrev, y = Concentration), size = 1) + 
  theme_bw() + theme(strip.text = element_text(face = "bold"),
                     axis.text = element_text(size = 8),
                     axis.title = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) +
  facet_wrap(~Chemistry, scales = "free_y") + labs(y = expression("W m"^-2), x = "") + 
  scale_x_discrete(labels=c("J", "F", "M", "A", "M", "J"))

p_zoops_phyto <- ggarrange(p_phyto, p_zoops, align = "hv")
p_temp_light <- ggarrange(p_light, p_temp, align ="hv")

p_zoops_phyto <- p_zoops_phyto + patchwork::plot_annotation(tag_levels = 'a') &
  theme(plot.tag.position = "topleft",
        plot.tag = element_text(size = 10, hjust = 0, vjust = 0)) &
  theme(plot.tag = element_text(face = 'bold')) + 
  theme(plot.margin = unit(c(0.1,0.9,0.1,0.1), 'lines'))

p_nutrients <- p_nutrients + plot_annotation(tag_levels = list(c("c"))) &
  theme(plot.tag.position = "topleft",
        plot.tag = element_text(size = 10, hjust = 0, vjust = 0)) &
  theme(plot.tag = element_text(face = 'bold')) + 
  theme(plot.margin = unit(c(0.1,0.9,0.1,0.1), 'lines'))

p_temp_light <- p_temp_light + plot_annotation(tag_levels = list(c("b"))) &
  theme(plot.tag.position = "topleft",
        plot.tag = element_text(size = 10, hjust = 0, vjust = 0)) &
  theme(plot.tag = element_text(face = 'bold')) + 
  theme(plot.margin = unit(c(0.,0.9,0.1,0.1), 'lines'))

#Arrange all plots
p_zoops_phyto_final <- ggarrange(p_zoops_phyto, p_temp_light, p_nutrients, ncol = 1)

#Write the figure as pdf and tiff
setwd(fig_dir)
pdf("Figure_S4.pdf", width = 6, height = 6) #Set new windo size and replot whatever plot you just made. 
p_zoops_phyto_final
dev.off()

tiff("Figure_S4.tiff", width = 13, height = 13, unit = "cm", res = 600) #Set new windo size and replot whatever plot you just made. 
p_zoops_phyto_final
dev.off()
setwd(dir)


#--------------------------------------------------------------------------------------------------------------
#                                 SUPPLEMENTARY FIGURE 5: ICE PHENOLOGY TAXA MASS CORRELATIONS
#------------------------------------------------------------------------------------------------------------
# TAXA MASS REGRESSIONS #
Regression_data <- Phytoplankton_after_ice_out_mean
Dominant_taxa <- Regression_data %>% filter(Month == 4) %>% group_by(Genus, Month) %>% reframe(mean_biomass = mean(mean_biomass)) %>% arrange(-mean_biomass)
#Replace NA with 0 to calculate mean
Regression_data$mean_biomass[is.na(Regression_data$mean_biomass)] <- 0
#Count each occurence of taxa to check that they are the same no. to start with
Regression_data_count <- Regression_data %>% group_by(Genus, Month) %>% mutate(n = n())
remove(Regression_data_count)
#When biomass is not 0, count no occasions
Regression_data <- Regression_data %>% group_by(Genus, Month) %>% mutate(Count = case_when(mean_biomass > 0 ~ "Yes", T ~ "No")) %>% group_by(Genus, Count, Month) %>%
  mutate(n = n())

#Calculate number of taxa without enough data
data_test <- Regression_data %>% filter(mean_biomass > 0) %>% filter(Month == 4) 
data_test <- unique(data_test$Genus)

#Filter when biomass > 0 are more than 5 occasions to run correlations
Regression_data <- Regression_data %>% group_by(Genus, Month) %>% filter(n >= 5 & Count == "Yes") %>%
  mutate(Ice_off_DOY = as.numeric(Ice_off_DOY),
         Ice_on_DOY = as.numeric(Ice_on_DOY)) %>%
  select(Month, Genus, Winter_year,Ice_on_DOY, Ice_off_DOY, `Length Ice Cover (days)`,n, mean_biomass)
#Make into long format and remove unknown genera
Regression_data <- Regression_data %>% gather(key = "Ice_phenology", value = "DOY_days", c(Ice_on_DOY:`Length Ice Cover (days)`)) %>%
  filter(Genus != "Other")
#Run correlations per taxa. Ice period is correlated to both ice-on and ice-off, only including the latter
correlation_results <- Regression_data %>%
  nest(data = -c(Genus, Ice_phenology, Month)) %>%
  mutate(cor=map(data,~cor.test(.x$DOY_days, .x$mean_biomass, method = "pearson"))) %>%
  mutate(tidied = map(cor, tidy)) %>%
  unnest(tidied) %>%
  select(-data, -cor)

Data_for_plotting <- left_join(correlation_results, Regression_data)
Data_for_plotting$Ice_phenology <- factor(Data_for_plotting$Ice_phenology, levels = c("Length Ice Cover (days)", "Ice_on_DOY", "Ice_off_DOY"))
Data_for_plotting <- Data_for_plotting %>% group_by(Genus, Month, Ice_phenology) %>% mutate(Max_biomass = max(mean_biomass),
                                                                                            Max_DOY = max(DOY_days),
                                                                                            Min_DOY = min(DOY_days))
Data_for_plotting$mean_biomass# <- Data_for_plotting$mean_biomass/1000
library(ggpubr)
Data_for_plotting <- Data_for_plotting %>% group_by(Month, Ice_phenology, Genus) %>% arrange(estimate)
Taxa_regressions <- ggplot(filter(Data_for_plotting, Ice_phenology != "Length Ice Cover (days)" & Month == 4), mapping = aes(x = DOY_days, y = mean_biomass, col = Ice_phenology, shape = Ice_phenology)) + 
  facet_wrap(~Genus, scales = "free_y") + lims(y = c(0,NA)) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", se = F, linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c("#000099", "#660000"), labels = c("Ice-on", "Ice-off")) + 
  scale_shape_manual(values = c(17,19), labels = c("Ice-on", "Ice-off")) +
  stat_cor(cor.coef.name = "r", method = "pearson", size = 3, show.legend = F) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"),
        legend.text = element_text(size = 12)) +
  labs(x = "Day-of-the-year (DOY)", y = expression("Biomass ("*mu*"g l"^-1*")"), col = "", shape = "", size = 14) +
  guides(color = guide_legend(override.aes = list(size = 3)))

setwd(fig_dir)
tiff("Figure_S5.tiff", width = 30, height = 17, units = "cm", res = 600)
Taxa_regressions
dev.off()
setwd(dir)

#Make a list with all unique taxa below snow
Phytoplankton_list <- Phytoplankton_below_ice %>% select(Phylum, Genus, Functional_group) %>% distinct() %>%
  filter(Genus != "Other")

Phytoplankton_list$"Taxonomic group" <- Phytoplankton_list$Phylum

#Make table with unique taxa and their correlation to snow depth
correlation_results_table <- correlation_results %>% ungroup() %>% select(Ice_phenology, Genus, estimate, p.value) 

#join list with unique taxa (and phylum) to the correlation results
Correlation_table_taxa <- left_join(correlation_results_table, Phytoplankton_list, by = "Genus")

#Round to 3 sign. figures and rename the columns
Correlation_table_taxa <- Correlation_table_taxa %>% 
  group_by(Genus) %>% mutate() %>%
  mutate(p.value = round(p.value, 3),
         estimate = round(estimate, 3)) %>% select(`Taxonomic group`, Genus, Ice_phenology, estimate, p.value) %>%
  arrange(`Taxonomic group`, estimate)
names(Correlation_table_taxa) <- c("Taxonomic group", "Genus", "Ice_phenology", "Correlation coefficient", "p-value")

#Write table
# setwd(table_dir)
# write_xlsx(Correlation_table_taxa, "Ice_phenology_Genus_correlations.txt")
# setwd(dir)

#--------------------------------------------------------------------------------------------------------------
#                                 SUPPLEMENTARY FIGURE 6: SNOW NMDS AND PERMANOVA
#------------------------------------------------------------------------------------------------------------
#Join phytoplankton data and snow data
Mean_phytoplankton_genera_filter <- Mean_phytoplankton_genera %>% select(-Functional_group)

Phyto_snow <- left_join(Mean_phytoplankton_genera_filter, Mean_snow, by = c("Winter_year", "Month"))
#Remove NA
Phyto_snow <- Phyto_snow %>% drop_na(Mean_snow)
Phyto_snow$Taxa <- Phyto_snow$Genus
Phyto_snow <- Phyto_snow %>% ungroup() %>% select(-Genus)

#Make 0 to NA and add column for month abbreviation
NMDS_data <- Phyto_snow
NMDS_data[NMDS_data == 0] <- NA
df <- NMDS_data
df <- df %>% filter(Month == 2 | Month == 3 | Month == 1)
df$Month_new <- month.abb[df$Month]
df$Month_new <- factor(df$Month_new, levels = c("Jan", "Feb", "Mar"))
df <- add_column(df, Month_abb = df$Month_new, .before = 1)
df <- df %>% select(-Month_new)

# #Calculate percentage of unknown taxa
df <- df %>% spread(Taxa, Mean_biomass) %>% drop_na(Mean_snow)
Phyto_snow <- hybrd.rplc_if(df)
PERMANOVA_data <- Phyto_snow
df_filter <- df
df_filter <- df_filter %>% select(-Other)
df_filter <- df_filter[,colSums(!is.na(df_filter)) > 0*nrow(df_filter)] #Filter data based on percentage of non-NA cases (>38% data)
Phyto_snow <- hybrd.rplc_if(df_filter)
Phyto_snow_filter <- Phyto_snow[rowSums(Phyto_snow[,c(5:ncol(Phyto_snow))])>0,] #Remove samplings without any data
NMDS_spec <- Phyto_snow_filter[, c(5:ncol(Phyto_snow_filter))]
NMDS_env <- Phyto_snow_filter[,c(1:2,4)]
sp.hel <- decostand(NMDS_spec, method="hellinger")
NMDS <- metaMDS(sp.hel, autotransform = F, k=2, trymax = 100)
NMDS
#Fitting snow depth as env vector
ef <- envfit(NMDS, NMDS_env[c(1,3)], permu = 9999)
ef

# PLOTTING
data.scores <- as.data.frame(scores(NMDS)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$Snow_depth <- NMDS_env$Mean_snow
data.scores$Year <- NMDS_env$Winter_year
data.scores$Month <- NMDS_env$Month_abb
head(data.scores)
species.scores <- as.data.frame(scores(NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))
vec.df$variables <- rownames(vec.df)

p9 <- ggplot() + geom_hline(yintercept=0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept=0, linetype = "dashed", color = "gray50") + 
  geom_point(data=data.scores, mapping = aes(x=NMDS1, y=NMDS2, col = Snow_depth, shape = Month), size=3.5) + 
  coord_equal() + 
  theme_bw() + 
  scale_colour_gradient2(low = "#FDB863",
                         mid = "#B2ABD2",
                         high = "#542788",
                         midpoint = 22) + labs(
                           shape = NULL, colour = "Snow depth (cm)",
                           breaks=c(0,10,20,30,40,50,60)) + 
  geom_segment(data = vec.df,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.1, "cm")),
               colour="black",
               inherit.aes = FALSE) +
  geom_text(data = vec.df,
            aes(x = NMDS1, y=NMDS2, label = "Snow depth"),
            size=3.5)

p10 <- ggplot(data=species.scores, aes(x=NMDS1, y=NMDS2)) + 
  geom_text(mapping = aes(label=species), alpha=0.6, size = 2.5) + 
  geom_vline(xintercept=0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept=0, linetype = "dashed", color = "gray50") +
  coord_equal() + theme_bw()

p9 <- p9 +
  theme(axis.title = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top") 
p10 <- p10 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top")

p_final_snow_depth_NMDS <- (p9/p10) +
  plot_annotation(tag_levels = list(c("")),
                  tag_sep = '',
                  tag_prefix = '',
                  tag_suffix = '') &
  theme(plot.tag.position = "topleft",
        plot.tag = element_text(size = 8, hjust = 0, vjust = 0)) &
  theme(plot.tag = element_text(face = 'bold'),
        legend.position = "top")


setwd(fig_dir)
tiff("Figure_S6.tiff", units="cm", width=15, height=15, res=600)
p_final_snow_depth_NMDS
dev.off()
setwd(dir)

# PERMANOVA
Permanova_env <- NMDS_env %>% mutate(Month = as.numeric(Month_abb))
sp.hel <- decostand(NMDS_spec, method="hellinger")
Permanova_snow_depth <- adonis2(sp.hel ~ Mean_snow*Month, data=Permanova_env, permutations=99)
Explanatory <- rownames(Permanova_snow_depth)
Permanova_snow_depth <- as.data.frame(Permanova_snow_depth[1:5]) %>% mutate(SumOfSqs = round(SumOfSqs, 3),
                                                                            R2 = round(R2, 3),
                                                                            `F-statistics` = round(F, 3),
                                                                            `p-value` = round(`Pr(>F)`, 3))
Permanova_snow_depth <- cbind(Explanatory,Permanova_snow_depth)
# write_xlsx(Permanova_snow_depth, paste(table_dir, "Permanova_snow_depth.xlsx", sep = ""))

#--------------------------------------------------------------------------------------------------------------
#                                 RESULTS: CALCUATE MEANS AND CORRELATIONS
#------------------------------------------------------------------------------------------------------------
#Calculate mean ice phenology
Mean_ice_period <- Ice_cover_period %>% mutate(DOY_ice_on = yday(`Beginning of ice cover`),
                                               DOY_ice_off = mean(yday(`Ice break up`)),
                                               Mean_ice_period = mean(`Length Ice Cover (days)`),
                                               SD_ice_period = sd(`Length Ice Cover (days)`)) %>%
  mutate(DOY_ice_on_mean = case_when(DOY_ice_on > 300 ~ DOY_ice_on-365, T~DOY_ice_on),
         DOY_ice_on_mean_mean = mean(DOY_ice_on_mean),
         SD_ice_on = sd(DOY_ice_on_mean),
         SD_ice_off = sd(yday(`Ice break up`)))

as.Date(Mean_ice_period$DOY_ice_on_mean)
as.Date(Mean_ice_period$DOY_ice_off)

Ice_phenology_data <- Ice_cover_period %>% mutate(DOY_ice_on = yday(`Beginning of ice cover`),
                                                  DOY_ice_off = yday(`Ice break up`),
                                                  Mean_ice_period = mean(`Length Ice Cover (days)`)) %>%
  mutate(DOY_ice_on_new = case_when(DOY_ice_on > 300 ~ DOY_ice_on-365, T~DOY_ice_on))

cor.test(Ice_phenology_data$DOY_ice_on_new, Ice_phenology_data$DOY_ice_off)


#Calculate mean ice and snow thickness
#Join snow and ice thickness
Snow_ice_data <- left_join(Mean_snow, Mean_ice, by = c("Month", "Winter_year")) %>% drop_na()

Snow_ice_data <- Snow_ice_data %>% group_by(Winter_year) %>% reframe(Mean_snow = mean(Mean_snow),
                                                                     Mean_ice = mean(Mean_ice))

cor.test(Snow_ice_data$Mean_snow, Snow_ice_data$Mean_ice)

Mean_sd <- Snow_ice_data %>%
  reframe(Mean_ice_new = mean(Mean_ice),
          Mean_snow_new = mean(Mean_snow),
          sd_snow = sd(Mean_snow),
          sd_ice = sd(Mean_ice),
          n = n())

p_snow_ice_cor <- ggplot(Snow_ice_data, aes(Mean_ice, y=Mean_snow)) + 
  geom_smooth(method = "lm", se = F, col = "black", linetype = "dashed") +
  geom_point(size = 3) + theme_bw() + 
  labs(x = "Ice thickness (cm)", y = "Snow thickness (cm)") + 
  lims(x = c(0,NA),
       y = c(0,NA))

#Max ice thickness and correaltion to ice duration
Max_ice <- Ice_thickness %>% group_by(Winter_year) %>% reframe(Max_ice = max(Ice_thickness))

Cor_data <- left_join(Max_ice, Ice_cover_period, by = "Winter_year")
cor.test(Cor_data$Max_ice, Cor_data$`Length Ice Cover (days)`)

setwd(fig_dir)
tiff("Figure_S2.tiff", width = 6, height = 5, units = "cm", res = 500)
p_snow_ice_cor
dev.off()
setwd(dir)

setwd(fig_dir)
pdf("Figure_S2.pdf", width = 5, height = 4)
p_snow_ice_cor
dev.off()
setwd(dir)

#Calculate percentage unknowns
#Snow and ice data in Jan-Mar
Phyto_data <- Phytoplankton_below_ice %>% group_by(Genus, Winter_year, Month, Day, Functional_group) %>% dplyr::summarise(across(biomass_ug_l, list(sum))) %>%
  mutate(Functional_group = case_when(Functional_group == "Autrotroph" ~ "Autotroph", T~Functional_group))
#Calculate mean biomass per genera and month within the ice cover period
Mean_phytoplankton_genera <- Phyto_data %>% group_by(Genus, Winter_year, Month, Functional_group) %>% 
  reframe(Mean_biomass = mean(biomass_ug_l_1)) %>% ungroup()
Mean_phytoplankton_genera <- Mean_phytoplankton_genera %>% filter(Month == 1 | Month == 2 | Month == 3) %>% select(-Functional_group)

Unknown_genera <- Mean_phytoplankton_genera %>% filter(Genus == "Other") %>% reframe(sum(Mean_biomass))
Known_genera <- Mean_phytoplankton_genera %>% filter(Genus != "Other") %>% reframe(sum(Mean_biomass))
Perc <- (Unknown_genera/(Known_genera+Unknown_genera))*100

#Caulculate dominant taxa in April after ice-off
Dominant_taxa <- Phytoplankton_after_ice_out_mean %>% filter(Month_abb == "Apr") %>% group_by(Genus, Month_abb) %>% reframe(Mean_phyto_biomass = mean(mean_biomass))

# Calculate unknown taxa in biomass after ice-off data (Apr-Jun)
#Replace NA with 0
Phytoplankton_after_ice_out_mean <- hybrd.rplc_if(Phytoplankton_after_ice_out_mean)
#Order months
Phytoplankton_after_ice_out_mean$Month_abb <- factor(Phytoplankton_after_ice_out_mean$Month_abb, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
#Calculate percentage
Known_genera <- Phytoplankton_after_ice_out_mean %>% filter(Genus != "Other") %>% reframe(Biomass = sum(mean_biomass))
Unknown_genera <- Phytoplankton_after_ice_out_mean %>% filter(Genus == "Other") %>% reframe(Biomass = sum(mean_biomass))
Perc <- (Unknown_genera/(Known_genera+Unknown_genera))*100

#Calculate samplings with 0-7.5 m instead of 0-20 m.
Mean_chla <- Chla_below_ice %>% select(-Date, -Day) %>% group_by(Winter_year, Month, max_depth_m) %>%
  drop_na(`Chl a_µg/l`) %>% reframe(Mean_chla = mean(`Chl a_µg/l`)) %>% filter(Month == 1 | Month == 2 | Month == 3)
#--------------------------------------------------------------------------------------------------------------
#                                 RESULTS: ICE THICKNESS PERMANOVA
#------------------------------------------------------------------------------------------------------------
Mean_phytoplankton_genera_filter <- Mean_phytoplankton_genera
Phyto_ice <- left_join(Mean_ice, Mean_phytoplankton_genera_filter, by = c("Winter_year", "Month"))
#Remove NA
Phyto_ice <- Phyto_ice %>% drop_na(Mean_ice)
Phyto_ice$Taxa <- Phyto_ice$Genus
Phyto_ice <- Phyto_ice %>% ungroup() %>% select(-Genus)
#Make 0 to NA and add column for month abbreviation
NMDS_data <- Phyto_ice
NMDS_data[NMDS_data == 0] <- NA
df <- NMDS_data
df <- df %>% filter(Month == 2 | Month == 3 | Month == 1)
df$Month_new <- month.abb[df$Month]
df$Month_new <- factor(df$Month_new, levels = c("Jan", "Feb", "Mar"))
df <- add_column(df, Month_abb = df$Month_new, .before = 1)
df <- df %>% select(-Month_new)
df <- df %>% spread(Taxa, Mean_biomass) %>% drop_na(Mean_ice)
#NMDS (removing unknowns and samplings with no obs)
df <- hybrd.rplc_if(df)
df_filter <- df
df_filter <- df_filter %>% select(-Other)
df_filter <- df_filter[,colSums(!is.na(df_filter)) > 0*nrow(df_filter)] #Filter data based on percentage of non-NA cases (>38% data)
Phyto_ice <- hybrd.rplc_if(df_filter)
Phyto_ice_filter <- Phyto_ice[rowSums(Phyto_ice[,c(5:ncol(Phyto_ice))])>0,] #Remove samplings without any data
NMDS_spec <- Phyto_ice_filter[, c(5:ncol(Phyto_ice_filter))]
NMDS_env <- Phyto_ice_filter[,c(1:2,4)]

# PERMANOVA
Permanova_env <- NMDS_env %>% mutate(Month = as.numeric(Month_abb))
sp.hel <- decostand(NMDS_spec, method="hellinger")
Permanova_ice_thickness <- adonis2(NMDS_spec ~ Mean_ice*Month, data=Permanova_env, permutations=99)
Explanatory <- rownames(Permanova_ice_thickness)
Permanova_ice_thickness <- as.data.frame(Permanova_ice_thickness[1:5]) %>% mutate(SumOfSqs = round(SumOfSqs, 3),
                                                                                  R2 = round(R2, 3),
                                                                                  `F-statistics` = round(F, 3),
                                                                                  `p-value` = round(`Pr(>F)`, 3))
Permanova_ice_thickness <- cbind(Explanatory,Permanova_ice_thickness)
#write_xlsx(Permanova_snow_depth, paste(table_dir, "Permanova_snow_depth.xlsx", sep = ""))

#--------------------------------------------------------------------------------------------------------------
#                                 RESULTS: DOMINANT TAXA BELOW SNOW
#------------------------------------------------------------------------------------------------------------
#Replace NA with 0
Phytoplankton_data_snow <- hybrd.rplc_if(Phytoplankton_below_ice)
#Calculate sum of biomass per taxa and day
Phyto_data <- Phytoplankton_data_snow %>% group_by(Genus, Winter_year, Month, Day, Date) %>% dplyr::summarise(across(biomass_ug_l, list(sum))) %>% ungroup()
#Find unique taxa in the whole data set and assign to make sure all taxa are present (0) when calculating mean
Unique_data <- Phyto_data %>% select(Date, Winter_year, Month, Day, Genus)
Unique_data$count <- NA
Unique_data <- Unique_data %>% complete(Date, Genus, fill = list(count = NA))
#Bind unique list of taxa for each sampling dat with the phytoplankton data (every taxa not present gets NA)
Phyto_data_final <- left_join(Unique_data, Phyto_data, by = c("Genus", "Date", "Winter_year", "Month", "Day")) %>%
  mutate(Winter_year = as.numeric(Winter_year),
         Month = as.numeric(Month))
remove(Phyto_data, Unique_data)
#Choose only taxa in the months Jan-Mar
Phyto_data_final <- Phyto_data_final %>% filter(Month == 1 | Month == 2 | Month == 3)
#Separate dates
Snow_data <- separate(Snow_on_ice, col = "Datum", into = c("Year", "Month", "Day"))
#Calculate mean snow depth per month and winter year
Snow_data <- Snow_data %>% group_by(Winter_year, Month) %>% reframe(Snow_depth = mean(Snow_depth)) %>%
  mutate(Month = as.numeric(Month))
#Replace NA with 0 to calculate mean
Phyto_data_final$biomass_ug_l_1[is.na(Phyto_data_final$biomass_ug_l_1)] <- 0
#Calculate mean biomass per month and winter year
Phyto_data_final <- Phyto_data_final %>% group_by(Genus, Winter_year, Month) %>% reframe(mean_biomass = mean(biomass_ug_l_1)) 
#Join phytoplankton data and snow data
Phyto_snow <- left_join(Phyto_data_final, Snow_data, by = c("Winter_year", "Month")) %>% drop_na(Snow_depth)
remove(Phyto_data_final, Snow_data)
#Make new column called taxa from genus
Phyto_snow$Taxa <- Phyto_snow$Genus
Phyto_snow <- Phyto_snow %>% select(-Genus)
#Find which taxa are present when snow depth is not 0
Phyto_snow <- Phyto_snow %>% group_by(Taxa) %>% mutate(Keep = case_when(Snow_depth > 0 ~ "Yes", T~"No"))
#Filter the taxa present with a biomass bigger than 0 and make a vector of the taxa list
vector_snow_taxa <- Phyto_snow %>% filter(mean_biomass > 0 & Keep == "Yes") %>% distinct(Taxa)
vector_snow_taxa <- as.vector(vector_snow_taxa$Taxa)
#Filter the taxa which are present below snow in the data
Phyto_snow <- Phyto_snow %>% filter(Taxa %in% vector_snow_taxa)
#Calculate mean biomass and snow depth per year (based on Jan-Mar)
Phyto_snow_annual_mean <- Phyto_snow %>% ungroup() %>% group_by(Winter_year, Taxa) %>% reframe(Biomass = mean(mean_biomass),
                                                                                               Snow_depth = mean(Snow_depth)) %>% select(Winter_year, Taxa, Snow_depth, Biomass)
remove(Phyto_snow)
#Assign to new name
Regression_data <- Phyto_snow_annual_mean
remove(Phyto_snow_annual_mean)
#Find dominant taxa below ice
Dominant_taxa <- Regression_data %>% group_by(Taxa) %>% reframe(mean_biomass = mean(Biomass)) %>% arrange(-mean_biomass)

