#--------------------------------------------------------------------------------------------------------------
#                                 META DATA
#------------------------------------------------------------------------------------------------------------
# Script by: Ellinor Jakobsson, 2023-2024
# Script for plotting all figures and doing all analyses in manuscript:
#"Effects of changing snow- and ice cover conditions on phytoplankton biomass and community composition in a mesotrophic lake"

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
Snow_on_ice <- read_excel(paste(dir, "Snow_data_on_ice.xlsx", sep = ""), sheet = 2)

#Ice thickness data 1999-2019
Ice_thickness <- read_excel(paste(dir, "Ice_thickness.xlsx", sep = ""), sheet = 2)

#Ice cover period data 1997-2019
Ice_cover_period <- read_excel(paste(dir, "/Ice_cover_period.xlsx", sep = ""), sheet = 2)

#Phytoplankton data below ice 1997-2019
Phytoplankton_below_ice <- read_excel(paste(dir, "/Phytoplankton_below_ice.xlsx", sep = ""), sheet = 2)

#Phytoplankton biomass after ice-out 1997-2019
Phytoplankton_after_ice_out <- read_excel(paste0(dir, "/Phytoplankton_after_ice_out.xlsx"), sheet = 2)

#Chla data after ice-out 1997-2019
Chla_after_ice_out <- read_excel(paste0(dir, "/Chla_after_ice_out.xlsx"), sheet = 2)

#Write chla within ice-period
Chla_below_ice <- read_excel(paste0(dir, "/Chla_below_ice.xlsx"), sheet = 2)

#Write summary data
Zooplankton_nutrient_light <- read_excel(paste0(dir, "/Zooplankton_nutrient_light.xlsx"), sheet = 2) 

#Precipitation data from all stations
#Data from Svanberga precip and Vallnora precip
Vallnora_precip <- read_excel(paste0(dir, "/Vallnora_precip.xlsx"), sheet = 2)

#Data from Svanberga precip and Norrveda precip
Norrveda_precip <- read_excel(paste0(dir, "/Norrveda_precip.xlsx"), sheet = 2) 

#Useful functions and defined parameters for the script
#Read function to replace NA with 0. Used later on.
hybrd.rplc_if <- function(x) { mutate_if(x, is.numeric, ~replace(., is.na(.), 0)) }

#Set within which dates to filter i.e., ice on and ice off dates from ice cover period data
Lower_range <- as.Date(Ice_cover_period$`Beginning of ice cover`)
Upper_range <- as.Date(Ice_cover_period$`Ice break up`)

#Format to date format
Snow_on_ice$Date <- as.Date(Snow_on_ice$Date)

#Calculate monthly mean of snow and chla per year and transform snow depth to cm
Mean_snow <- Snow_on_ice %>% select(-Day, -Date) %>% group_by(Winter_year, Month) %>% drop_na(Snow_depth) %>%
  reframe(Mean_snow = mean(Snow_depth)*100) #snow depth in cm

Mean_chla <- Chla_below_ice %>% select(-Date, -Day) %>% group_by(Winter_year, Month) %>%
  drop_na(`Chl a_µg/l`) %>% reframe(Mean_chla = mean(`Chl a_µg/l`)) #chla in ug L-1

Mean_ice <- Ice_thickness %>% group_by(Winter_year, Month) %>% drop_na(Ice_thickness) %>%
  reframe(Mean_ice = mean(Ice_thickness)) #Ice thickness in cm

#--------------------------------------------------------------------------------------------------------------
#                                 FIGURE 1: PHYTOPLANKTON CHL-a ~ SNOW DEPTH
#------------------------------------------------------------------------------------------------------------
#Join mean snow and mean chla
Snow_chla_data <- left_join(Mean_chla, Mean_snow, by = c("Winter_year", "Month"))
Snow_chla_data$Month <- as.numeric(Snow_chla_data$Month)
Snow_chla_data$Month_abbrev <- month.abb[Snow_chla_data$Month]
Snow_chla_data$Month_abbrev <- factor(Snow_chla_data$Month_abbrev, levels = c("Dec", "Jan", "Feb", "Mar", "Apr"))
Snow_chla_data$Winter_year <- as.numeric(Snow_chla_data$Winter_year)

#------------------------------------------------------Plot snow and chl-a
#January
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

#Extract model params
p_val <- round((summary(m_1)$coefficients[,4][2]*2), 3) #Multiply by 2 for Holms correction and round with three decimals
rsqrd <- round(summary(m_1)$adj.r.squared, 2)
coeff <- as.data.frame(cbind(p_val, rsqrd))

# generate position of labs for plot
labpos <- data %>% 
  ungroup() %>%
  reframe(maxheight = max(Mean_chla), 
            minheight = min(Mean_chla),
            xpos = mean(Mean_snow))

#bind the position of labs and the lm stats to be plotted
plotting_labs <- cbind(labpos, coeff)

p20 <- ggplot(data, aes(Mean_snow, y = Mean_chla)) +
  geom_smooth(mapping = aes(x = Mean_snow, y = Mean_chla), method = "lm", linetype = "dashed", se = F, col = "black", linewidth = 0.3) + 
  geom_point(mapping = aes(x = Mean_snow, y = Mean_chla), size = 2) + 
  lims(x = c(0,55), y=c(0,NA)) + theme_bw() + ylab("January") +
  geom_text(plotting_labs,
            mapping = aes(x = 20, y = minheight + 0.85 * (maxheight - minheight), label = paste0("p = ", format(round(p_val, 4), scientific = F), ",\nR² = ", round(rsqrd, 2))),
            family = "serif", fontface = "italic", size = 2.5)


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
#Extract model params
p_val <- round((summary(m_2)$coefficients[,4][2]*1), 3) #Multiply by 1 for Holms correction and round with three decimals
rsqrd <- round(summary(m_2)$adj.r.squared, 2)
coeff <- as.data.frame(cbind(p_val, rsqrd))

#Test if model is significant
summary(m_2)
data = filter(Snow_chla_data, Month_abbrev == "Feb")

# generate position of labs for plot
labpos <- data %>% 
  ungroup() %>%
  reframe(maxheight = max(Mean_chla), 
          minheight = min(Mean_chla),
          xpos = mean(Mean_snow))

#bind the position of labs and the lm stats to be plotted
plotting_labs <- cbind(labpos, coeff)

#Plot the data
p <- ggplot(data, aes(Mean_snow, Mean_chla))
p <- p + geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), se = F, color = "black", linewidth = 0.3, linetype = "dashed") + 
  geom_point(size = 2)
p21 <- p + lims(x=c(0,55)) + scale_y_continuous(trans = "log10", name = "February") +
  annotation_logticks(sides = "l") + theme_bw() + 
  geom_text(plotting_labs,
            mapping = aes(x = 20, y = minheight + 0.75 * (maxheight - minheight), label = paste0("p = ", format(round(p_val, 4), scientific = F), ",\nR² = ", round(rsqrd, 2))),
            family = "serif", fontface = "italic", size = 2.5)

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

# generate position of labs for plot
labpos <- data %>% 
  ungroup() %>%
  reframe(maxheight = max(Mean_chla), 
          minheight = min(Mean_chla),
          xpos = mean(Mean_snow))

#Extract model params
p_val <- (summary(m_3)$coefficients[,4][2])*3 #Multiply by 3 for Holms correction and round with three decimals
rsqrd <- round(summary(m_3)$adj.r.squared, 2)
coeff <- as.data.frame(cbind(p_val, rsqrd))

#bind the position of labs and the lm stats to be plotted
plotting_labs <- cbind(labpos, coeff)

#Plot everything
p <- ggplot(data, aes(Mean_snow, Mean_chla))
p <- p + geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), se = F, color = "black", linewidth = 0.3, linetype = "dashed") +
  geom_point(size = 2)
p22 <- p + #geom_function(fun = function(x) (24.905*exp(-0.07*x)), colour = "black", linetype=2) +  
  scale_y_continuous(trans = "log10", name = "March") + #,
                     #sec.axis = sec_axis(trans=~.*1)) + 
                     theme_bw() +
  annotation_logticks(sides = "l") + lims(x=c(0,55)) +
  geom_text(plotting_labs,
            mapping = aes(x = 20, y = minheight + 0.68 * (maxheight - minheight), label = paste0("p < ", format(round(p_val, 4), scientific = F), ",\nR² = ", round(rsqrd, 2))),
            family = "serif", fontface = "italic", size = 2.5)

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
tiff("Figure_1.tiff", unit = "cm", height = 4, width = 15, res = 1200)
p_final_production_snow_depth
dev.off()
setwd(dir)

#--------------------------------------------------------------------------------------------------------------
#                                 FIGURE 2: CHL-a ~ ICE PHENOLOGY SCATTERPLOT 
#------------------------------------------------------------------------------------------------------------
#Set the correct order of months
Chla_after_ice_out$Month_abb <- factor(Chla_after_ice_out$Month_abb, levels = c("Apr", "May", "Jun"))
 
#Regressions on ice phenology
# Rsqr values
Ice_period_table_r <- Chla_after_ice_out %>% group_by(Month) %>% lm_table(Mean_chla ~ `Length Ice Cover (days)`) %>%
  select(Month, Rsqr_adj) %>% rename("R2" = "Rsqr_adj") %>%
  mutate(R2 = round(R2, 4))
# p values
Ice_period_table_p <- Chla_after_ice_out  %>% nest(data = -c(Month)) %>%
  mutate(model = map(data, ~lm(Mean_chla ~ `Length Ice Cover (days)`, data = .)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>% filter(term == "`Length Ice Cover (days)`") %>%
  select(Month, p.value) %>%
  rename("p" = "p.value")
#---------------------------ICE ON DOY
# Rsqr values
Ice_on_table_r <- Chla_after_ice_out %>% group_by(Month) %>% lm_table(Mean_chla ~ Ice_on_DOY) %>%
  select(Month, Rsqr_adj) %>% rename("R2" = "Rsqr_adj") %>%
  mutate(R2 = round(R2, 4))
# p values
Ice_on_table_p <- Chla_after_ice_out  %>% nest(data = -c(Month)) %>%
  mutate(model = map(data, ~lm(Mean_chla ~ Ice_on_DOY, data = .)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>% filter(term == "Ice_on_DOY") %>% 
  select(Month, p.value) %>%
  rename("p" = "p.value") %>% mutate(p = round(p, 2))
#---------------------------ICE OFF DOY
# Rsqr values
Ice_off_table_r <- Chla_after_ice_out %>% group_by(Month) %>% lm_table(Mean_chla ~ Ice_off_DOY) %>%
  select(Month, Rsqr_adj) %>% rename("R2" = "Rsqr_adj") %>%
  mutate(R2 = round(R2, 4))
# p values
Ice_off_table_p <- Chla_after_ice_out  %>% nest(data = -c(Month)) %>%
  mutate(model = map(data, ~lm(Mean_chla ~ Ice_off_DOY, data = .)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>% filter(term == "Ice_off_DOY") %>% 
  select(Month, p.value) %>%
  rename("p" = "p.value") %>% mutate(p = round(p, 2))
#----------------------------COMBINE TABLES
# combine
Ice_period_table <- left_join(Ice_period_table_r, Ice_period_table_p)
Ice_on_table <- left_join(Ice_on_table_r, Ice_on_table_p)
Ice_off_table <- left_join(Ice_off_table_r, Ice_off_table_p)
Ice_off_table$Explanatory <- "Ice_off_DOY"
Ice_on_table$Explanatory <- "Ice_on_DOY"
Ice_period_table$Explanatory <- "Ice_period"
Ice_period_table <- Ice_period_table %>% select(Explanatory, Month, R2, p)
Ice_on_table <- Ice_on_table %>% select(Explanatory, Month, R2, p)
Ice_off_table <- Ice_off_table %>% select(Explanatory, Month, R2, p)
Ice_phenology_lm_table <- rbind(Ice_period_table, Ice_on_table, Ice_off_table)
Ice_on_table$Month_abb <- month.abb[Ice_on_table$Month]
Ice_phenology_lm_table$Month_abb <- month.abb[Ice_phenology_lm_table$Month]

Chla_after_ice_out_long <- Chla_after_ice_out %>%
  select(Winter_year, Month_abb, Mean_chla, Ice_on_DOY, Ice_off_DOY, `Length Ice Cover (days)`) %>%
  group_by(Winter_year, Month_abb) %>% 
  gather(Parameter, Days_doy, -Mean_chla, -Winter_year, -Month_abb)

# Labs for plot
labpos <- Chla_after_ice_out_long %>% 
  group_by(Parameter, Month_abb) %>% 
  summarise(maxheight = max(Mean_chla), 
            minheight = min(Mean_chla),
            xpos = mean(Days_doy*0.90)) %>%
  select(Parameter, Month_abb, maxheight, minheight, xpos) %>%
  mutate(Parameter = case_when(Parameter == "Length Ice Cover (days)"~"Ice_period",
                               T~Parameter))

#Change name to prepare for joing of datasets
Ice_phenology_lm_table$Parameter <- Ice_phenology_lm_table$Explanatory
Ice_phenology_lm_table <- Ice_phenology_lm_table %>% select(Parameter, Month_abb, R2, p)

#bind the position of labs and the lm stats to be plotted
plotting_labs <- left_join(labpos, Ice_phenology_lm_table)

plotting_labs <- plotting_labs
plotting_labs <- plotting_labs

# Put months in correct order
plotting_labs$Month_abb <- factor(plotting_labs$Month_abb, levels = c("Apr", "May", "Jun"))
Chla_after_ice_out$Month_abb <- factor(Chla_after_ice_out$Month_abb, levels = c("Apr", "May", "Jun"))

#Regression plots on ice on and ice off DOY

p_ice_phenology <- ggplot() + geom_point(Chla_after_ice_out, mapping = aes(x = Ice_off_DOY, y = Mean_chla), col = "#660000", size = 1, shape = 16) + 
  geom_point(Chla_after_ice_out, mapping = aes(x = Ice_on_DOY, y = Mean_chla), col = "#000099", size = 1, shape = 17) +
  geom_smooth(filter(Chla_after_ice_out, Month_abb == "Apr"), mapping = aes(x = Ice_off_DOY, y = Mean_chla), col = "#660000", method = "lm", se = F, alpha = 0.5, linewidth = 0.3, linetype="dashed") + 
#  geom_smooth(filter(Chla_after_ice_out), mapping = aes(x = Ice_on_DOY, y = Mean_chla), col = "#000099", method = "lm", se = F, alpha = 0.5, linewidth = 0.3, linetype="dashed") + 
  facet_wrap(~Month_abb, scales = "free_y") + 
  theme_bw() + 
  labs(y = expression(paste("Chl-"*alpha*" ("*mu*"g "*L^-1*")"))) + 
  theme(strip.text = element_text(face = "bold")) + 
  labs(x = "Day-of-the-year") +
  geom_text(filter(plotting_labs, Parameter != "Ice_period"),
            mapping = aes(x = xpos, y = minheight + 0.85 * (maxheight - minheight), label = paste0("p = ", round(p, 3), ",\nR² = ", format(round(R2, 2), scientific = F ))),
            family = "serif", fontface = "italic", size = 2.5)


p_Ice_period <- ggplot() + 
  geom_smooth(filter(Chla_after_ice_out, Month_abb == "Apr"), mapping = aes(x = `Length Ice Cover (days)`, y = Mean_chla), col = "black",method = "lm", se = F, linewidth = 0.5, alpha = 0.5, linetype="dashed") +
  geom_point(Chla_after_ice_out, mapping = aes(x = `Length Ice Cover (days)`, y = Mean_chla), col = "black", size = 1) +
  facet_wrap(~Month_abb, scales = "free_y") + theme_bw() +
  labs(y = expression(paste("Chl-"*alpha*" ("*mu*"g "*L^-1*")")), x = "Ice cover period (days)") + 
  theme(strip.text = element_text(face = "bold"))  +
  geom_text(filter(plotting_labs, Parameter == "Ice_period"),
            mapping = aes(x = xpos*0.5, y = minheight + 0.85 * (maxheight - minheight), label = paste0("p = ", round(p, 2), ",\nR² = ", format(round(R2, 2), scientific = F))),
            family = "serif", fontface = "italic", size = 2.5)

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
tiff("Figure_2.tiff", width=12, height=9, unit = "cm", res = 1200)
p_final_DOY_box 
dev.off()


#--------------------------------------------------------------------------------------------------------------
#                                 FIGURE 3: TROPHIC MODES ~ SNOW DEPTH
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
#Clean up data file by calculating percentage from proportion and remove unused columns
Perc <- add_column(Perc, Percentage = Perc$Perc*100, .before = 1)
Perc <- Perc %>% select(-Perc, -Sum_total, -sum_feeding)
#Calculate mean percentage per month
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

# #-------------------------------Beta regression on trophic modes
# #PERCENTAGE
library(betareg)

#Calculate the proportion for the test
Perc2$Proportion <- Perc2$Percentage/100

# Beta regression on proportion of trophic modes
#January
betareg_trophicmode <- betareg(Proportion ~ Functional_group*Mean_snow, link = "logit", data = filter(Perc2, Functional_group != "Heterotrophs" & Month == 1))
betareg_trophicmode1 <- as.data.frame(summary(betareg_trophicmode)$coefficients$mean[-1,])
betareg_trophicmode1$Month <- "Jan"
betareg_trophicmode1$row <- rownames(betareg_trophicmode1)

#February
betareg_trophicmode <- betareg(Proportion ~ Functional_group*Mean_snow, link = "logit", data = filter(Perc2, Functional_group != "Heterotrophs" & Month == 2))
betareg_trophicmode2 <- as.data.frame(summary(betareg_trophicmode)$coefficients$mean[-1,])
betareg_trophicmode2$Month <- "Feb"
betareg_trophicmode2$row <- rownames(betareg_trophicmode2)

#March
betareg_trophicmode <- betareg(Proportion ~ Functional_group*Mean_snow, link = "logit", data = filter(Perc2, Functional_group != "Heterotrophs" & Month == 3))
betareg_trophicmode3 <- as.data.frame(summary(betareg_trophicmode)$coefficients$mean[-1,])
betareg_trophicmode3$Month <- "Mar"
betareg_trophicmode3$row <- rownames(betareg_trophicmode3)


#Adjust for multiple testing of interaction effect
final_results_trophic_modes <- rbind(betareg_trophicmode1, betareg_trophicmode2, betareg_trophicmode3)
final_results_trophic_modes_interaction <- final_results_trophic_modes %>% filter(row == "Functional_groupMixotrophs:Mean_snow")
final_results_trophic_modes_interaction$p_adj <- p.adjust(final_results_trophic_modes_interaction$`Pr(>|z|)`, method = "holm") #Adjust multiple testing for the interaction

final_results_trophic_modes <- final_results_trophic_modes %>% filter(row != "Functional_groupMixotrophs:Mean_snow")
final_results_trophic_modes
final_results_trophic_modes_interaction

# ABSOLUTE BIOMASS
#---------------Jan
Mixo_lm <- lm(log(Mean_feeding)~Mean_snow, data = filter(Perc, Functional_group == "Mixotrophs" & Month == 1))
#shapiro.test(Mixo_lm$residuals)
summary(Mixo_lm) #0.39
#Mixo_lm1 <- as.data.frame(Mixo_lm$coefficients[,4])[2,]
Mixo_lm1 <- lm(log(Mean_feeding)~Mean_snow, data = filter(Perc, Functional_group == "Mixotrophs" & Month == 1))
p_val <- round((summary(Mixo_lm1)$coefficients[,4][2]), 3)
rsqrd <- round(summary(Mixo_lm1)$adj.r.squared, 2)
coeffs1 <- as.data.frame(cbind(p_val, rsqrd))
# Log transformed

Auto_lm <- lm(log(Mean_feeding)~Mean_snow, data = filter(Perc, Functional_group == "Autotrophs" & Month == 1))
#shapiro.test(Auto_lm$residuals)
summary(Auto_lm) #0.43
#Auto_lm1 <- as.data.frame(Auto_lm$coefficients[,4])[2,]
Auto_lm1 <- lm(log(Mean_feeding)~Mean_snow, data = filter(Perc, Functional_group == "Autotrophs" & Month == 1))
p_val <- round((summary(Auto_lm1)$coefficients[,4][2]), 3)
rsqrd <- round(summary(Auto_lm1)$adj.r.squared, 2)
coeffs2 <- as.data.frame(cbind(p_val, rsqrd))
# Log transformed


#----------------Feb
Mixo_lm <- lm(log(Mean_feeding)~Mean_snow, data = filter(Perc, Functional_group == "Mixotrophs" & Month == 2))
#shapiro.test(Mixo_lm$residuals)
summary(Mixo_lm) #0.92
#Mixo_lm2 <- as.data.frame(Mixo_lm$coefficients[,4])[2,]
Mixo_lm2 <- lm(log(Mean_feeding)~Mean_snow, data = filter(Perc, Functional_group == "Mixotrophs" & Month == 2))
p_val <- round((summary(Mixo_lm2)$coefficients[,4][2]), 3)
rsqrd <- round(summary(Mixo_lm2)$adj.r.squared, 2)
coeffs3 <- as.data.frame(cbind(p_val, rsqrd))
# Log transformed

Auto_lm <- lm(Mean_feeding~Mean_snow, data = filter(Perc, Functional_group == "Autotrophs" & Month == 2))
#shapiro.test(Auto_lm$residuals)
summary(Auto_lm) #0.02
#Auto_lm2 <- as.data.frame(Auto_lm$coefficients[,4])[2,]
Auto_lm2 <- lm(Mean_feeding~Mean_snow, data = filter(Perc, Functional_group == "Autotrophs" & Month == 2))
p_val <- round((summary(Auto_lm2)$coefficients[,4][2]), 3) 
rsqrd <- round(summary(Auto_lm2)$adj.r.squared, 2)
coeffs4 <- as.data.frame(cbind(p_val, rsqrd))
# Not log transformed 

#-----------------Mar
Mixo_lm <- lm(log(Mean_feeding)~Mean_snow, data = filter(Perc, Functional_group == "Mixotrophs" & Month == 3))
#shapiro.test(Mixo_lm$residuals)
summary(Mixo_lm) #0.51
#Mixo_lm3 <- as.data.frame(Mixo_lm$coefficients[,4])[2,]
Mixo_lm3 <- lm(log(Mean_feeding)~Mean_snow, data = filter(Perc, Functional_group == "Mixotrophs" & Month == 3))
p_val <- round((summary(Mixo_lm3)$coefficients[,4][2]), 3) 
rsqrd <- round(summary(Mixo_lm3)$adj.r.squared, 2)
coeffs5 <- as.data.frame(cbind(p_val, rsqrd))
# Log transformed

Auto_lm <- lm(log(Mean_feeding)~Mean_snow, data = filter(Perc, Functional_group == "Autotrophs" & Month == 3))
#shapiro.test(Auto_lm$residuals)
summary(Auto_lm) #0.08
#Auto_lm3 <- as.data.frame(Auto_lm$coefficients[,4])[2,]
Auto_lm3 <- lm(log(Mean_feeding)~Mean_snow, data = filter(Perc, Functional_group == "Autotrophs" & Month == 3))
p_val <- round((summary(Auto_lm3)$coefficients[,4][2]), 3) #Multiply by 2 for Holms correction and round with three decimals
rsqrd <- round(summary(Auto_lm3)$adj.r.squared, 2)
coeffs6 <- as.data.frame(cbind(p_val, rsqrd))
# Log tranformed

# # Perform holms corrections
coefficients_all <- rbind(coeffs1, coeffs2, coeffs3, coeffs4, coeffs5, coeffs6)
coefficients_all$Functional_group <- c("Mixotrophs", "Autotrophs","Mixotrophs", "Autotrophs","Mixotrophs", "Autotrophs")
coefficients_all$Month <- c(1,1,2,2,3,3)
coefficients_all$p_adj <- p.adjust(coefficients_all$p_val, method = "holm")

# generate position of labs for plot
labpos <- filter(Perc, Functional_group != "Heterotrophs") %>% 
  ungroup() %>%
  group_by(Month, Functional_group) %>%
  reframe(maxheight = max(Mean_feeding), 
          minheight = min(Mean_feeding),
          xpos = 20) %>%
  mutate(Month = as.numeric(Month))

#bind the position of labs and the lm stats to be plotted
plotting_labs <- left_join(labpos, coefficients_all)
plotting_labs$Month_abb <- month.abb[plotting_labs$Month]
plotting_labs <- plotting_labs %>% mutate(p_val = case_when(p_val > 1~1, T~p_val))
#----------------------------------------------------PLOTTING
#Remove heterotrophs from plots
Perc <- Perc %>% filter(Functional_group != "Heterotrophs")
Perc$Month_abb <- factor(Perc$Month_abb, levels = c("Jan", "Feb", "Mar"))
plotting_labs$Month_abb <- factor(plotting_labs$Month_abb)
Perc <- left_join(Perc, plotting_labs)

add_logticks  <- function (base = 10, sides = "bl", scaled = TRUE, 
                           short = unit(0.1, "cm"), mid = unit(0.2, "cm"),  long = unit(0.3, "cm"), 
                           colour = "black",  size = 0.5, linetype = 1, alpha = 1, color = NULL, 
                           data =data.frame(x = NA),... )   {
  if (!is.null(color)) 
    colour <- color
  layer(geom = "logticks", params = list(base = base, 
                                         sides = sides, scaled = scaled, short = short, 
                                         mid = mid, long = long, colour = colour, size = size, 
                                         linetype = linetype, alpha = alpha, ...), 
        stat = "identity", data = data , mapping = NULL, inherit.aes = FALSE, position = "identity",
        show.legend = FALSE)
}


#Plot absolute biomass in relation to snow depth

p_absolute_mixotroph_Jan <- ggplot() + 
  geom_point(filter(Perc, Month_abb == "Jan"), mapping = aes(x = Mean_snow, y = Mean_feeding), size = 1.5) +
  theme_bw() + 
  labs(x = "", y = "", size = 10) +
  theme_bw() +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.x.top = element_blank()) +
  theme(strip.text = element_text(face = "bold")) + 
  geom_text(filter(Perc, Month_abb == "Jan"), mapping = aes(x = 40, y = 10, label = paste0("p = 1.00", ",\nR² = ", round(rsqrd, 2))),
          family = "serif", fontface = "italic", size = 2.5) +
  facet_grid2(Month_abb~`Trophic identity`, axes = "all") +
  scale_y_continuous(trans = "log10") +
  annotation_logticks(sides = "l") + lims(x=c(0,55)) + 
  theme(plot.margin = unit(c(0,0,-1,-1), 'lines'))

#T, R, B, L
p_absolute_mixotroph_Feb <- ggplot() + 
  geom_point(filter(Perc, Month_abb == "Feb"), mapping = aes(x = Mean_snow, y = Mean_feeding), size = 1.5) +
  labs(x = "", y = "", size = 10) +
  theme_bw() +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.x.top = element_blank()) +
  theme(strip.text = element_text(face = "bold")) +
  geom_text(filter(Perc, Month_abb == "Feb"), mapping = aes(x = 40, y = 35, label = paste0("p = ", format(round(p_adj, 2), scientific = F), ",\nR² = ", round(rsqrd, 2))),
          family = "serif", fontface = "italic", size = 2.5) +
  facet_grid2(Month_abb~Functional_group, scales = "free", independent = T) +
  lims(x=c(0,55)) +
  facetted_pos_scales(y = list(Functional_group == "Mixotrophs" ~ scale_y_continuous(trans = "log10"))) + 
  add_logticks(side = "l", data = data.frame(x = NA, Functional_group = "Mixotrophs")) + 
  theme(plot.margin = unit(c(-1,0,-1,-1), 'lines'))

p_absolute_mixotroph_Mar <- ggplot() + 
  geom_point(filter(Perc, Month_abb == "Mar"), mapping = aes(x = Mean_snow, y = Mean_feeding), size = 1.5) +
  theme_bw() + 
  labs(x = "Snow depth (cm)", y = "", size = 10) +
  theme_bw() +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.x.top = element_blank()) +
  theme(strip.text = element_text(face = "bold")) + 
  geom_text(filter(Perc, Month_abb == "Mar"), mapping = aes(x = 40, y = 15, label = paste0("p = ", format(round(p_adj, 2), scientific = F), ",\nR² = ", round(rsqrd, 2))),
            family = "serif", fontface = "italic", size = 2.5) +
  facet_grid2(Month_abb~`Trophic identity`, axes = "all") +
  scale_y_continuous(trans = "log10") +
  annotation_logticks(sides = "l") + lims(x=c(0,55)) + 
  theme(plot.margin = unit(c(-1,0,-1,-1), 'lines'))


#Plot proportions
Perc2 <- Perc2 %>% filter(`Trophic identity` != "Heterotrophs")

#Fit the betareg again
betareg_trophicmode1 <- betareg(Proportion ~ Functional_group*Mean_snow, link = "logit", data = filter(Perc2, Functional_group != "Heterotrophs" & Month == 1))

p_proportion_mixotroph_Jan <- ggplot() + 
  geom_point(filter(Perc2, Month_abb == "Jan"), mapping = aes(x = Mean_snow, y = Proportion), size = 1.5) +
  theme_bw() + 
  labs(x = "", y = "", size = 10) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold")) + 
  facet_grid2(Month_abb~`Trophic identity`, axes = "all")+ 
  theme(plot.margin = unit(c(-1,0,-1,-1), 'lines')) + 
  geom_line(filter(Perc2, Month == 1), mapping = aes(x = Mean_snow, y = predict(betareg_trophicmode1, filter(Perc2, Month == 1))), 
            linetype = "dashed") +
  # geom_text(filter(Perc2, Month_abb == "Jan"), mapping = aes(x = 40, y = 0.75, label = paste0("p = 0.04")),
  #           family = "serif", fontface = "italic", size = 2.5) + 
  lims(x = c(0, 55))
 

betareg_trophicmode2 <- betareg(Proportion ~ Functional_group*Mean_snow, link = "logit", data = filter(Perc2, Functional_group != "Heterotrophs" & Month == 2))

p_proportion_mixotroph_Feb <- ggplot() + 
  geom_point(filter(Perc2, Month_abb == "Feb"), mapping = aes(x = Mean_snow, y = Proportion), size = 1.5) +
  theme_bw() + 
  labs(x = "", y = "", size = 10) +
  theme_bw() +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.x.top = element_blank()) +
  theme(strip.text = element_text(face = "bold")) + 
  facet_grid2(Month_abb~`Trophic identity`, axes = "all")+ 
  theme(plot.margin = unit(c(-1,0,-1,-1), 'lines')) + 
  geom_line(filter(Perc2, Month == 2), mapping = aes(x = Mean_snow, y = predict(betareg_trophicmode2, filter(Perc2, Month == 2))), 
            linetype = "dashed") +
  # geom_text(filter(Perc2, Month_abb == "Feb"), mapping = aes(x = 40, y = 0.75, label = paste0("p < 0.001")),
  #           family = "serif", fontface = "italic", size = 2.5) + 
  lims(x = c(0, 55))

betareg_trophicmode3 <- betareg(Proportion ~ Functional_group*Mean_snow, link = "logit", data = filter(Perc2, Functional_group != "Heterotrophs" & Month == 3))


p_proportion_mixotroph_Mar <- ggplot() + 
  geom_point(filter(Perc2, Month_abb == "Mar"), mapping = aes(x = Mean_snow, y = Proportion), size = 1.5) +
  theme_bw() + 
  labs(x = "Snow depth (cm)", y = "", size = 10) +
  theme_bw() +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.x.top = element_blank()) +
  theme(strip.text = element_text(face = "bold")) + 
  facet_grid2(Month_abb~`Trophic identity`, axes = "all") + 
  theme(plot.margin = unit(c(-1,0,-1,-1), 'lines')) + 
  geom_line(filter(Perc2, Month == 3), mapping = aes(x = Mean_snow, y = predict(betareg_trophicmode3, filter(Perc2, Month == 3))), 
                linetype = "dashed") +
  # geom_text(filter(Perc2, Month_abb == "Mar"), mapping = aes(x = 40, y = 0.75, label = paste0("p = 0.016")),
  #           family = "serif", fontface = "italic", size = 2.5) + 
  lims(x = c(0, 55))



#Arrange plots
p_trophic_mode_final <- ggarrange(p_proportion_mixotroph_Jan,
          p_proportion_mixotroph_Feb,
          p_proportion_mixotroph_Mar,
          p_absolute_mixotroph_Jan,
          p_absolute_mixotroph_Feb,
          p_absolute_mixotroph_Mar,align = "hv", nrow = 6)

p_trophic_mode_final <- annotate_figure(p_trophic_mode_final, 
                                        left = text_grob(expression(paste("      Biomass  ("*mu*"g "*L^-1*")                ", "                                    Proportional biomass")), rot = 90, hjust = 0.5))

setwd(fig_dir)
pdf("Figure_3.pdf", height = 4, width=6)
p_trophic_mode_final
dev.off()

setwd(fig_dir)
tiff("Figure_3.tiff", height = 20, width=8, units = "cm", res = 1200)
p_trophic_mode_final
dev.off()


#--------------------------------------------------------------------------------------------------------------
#                                 FIGURE 4: AFTER ICE OUT CCA, DOMINANT TAXA IN APRIL
#------------------------------------------------------------------------------------------------------------
#Calculate dominant taxa in April after ice off
Dominant_taxa <- Phytoplankton_after_ice_out %>% filter(Month_abb == "Apr") %>% group_by(Genus, Month_abb, Month) %>% reframe(Mean_phyto_biomass = mean(mean_biomass))
# NMDS ice phenology after ice-out
#Replace NA with 0
Phytoplankton_after_ice_out <- hybrd.rplc_if(Phytoplankton_after_ice_out)
#Order months
Phytoplankton_after_ice_out$Month_abb <- factor(Phytoplankton_after_ice_out$Month_abb, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
NMDS_data_filter <- Phytoplankton_after_ice_out %>% spread(Genus, mean_biomass)
NMDS_data_filter <- NMDS_data_filter %>% filter(Month_abb == "Apr" | Month_abb == "May" | Month_abb == "Jun") #Select only Apr-Jun (after ice-out)

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
NMDS_env <- Phyto_snow_filter[,c(3:4, 11:13)]

library(ggrepel)
# CCA 
cca_env <- NMDS_env %>% mutate(Month = as.numeric(Month))
cca_ice_phenology <- cca(NMDS_spec ~ Ice_on_DOY*Ice_off_DOY+Month+Ice_on_DOY*Month+Ice_off_DOY*Month, data=cca_env)
cca_ice_phenology_anova_terms <- anova.cca(cca_ice_phenology, by = "term", permutations = 9999)
cca_ice_phenology_anova_model <- anova.cca(cca_ice_phenology, permutations = 9999)
cca_ice_phenology_anova_model
cca_ice_phenology_anova_terms

#CCA1: 0.3599504,
#CCA2: 0.2454282

#extracting the data as data frame; env data
veg_1 <- as.data.frame(cca_ice_phenology$CCA$biplot)
veg_1["env"] <- row.names(veg_1)
veg_1 <- veg_1 %>% mutate(env = case_when(env == "Ice_on_DOY"~"Ice-on",
                                          env == "Ice_off_DOY"~"Ice-off",
                                          T~env))

# PLOTTING
data.scores <- as.data.frame(scores(cca_ice_phenology)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$Ice_period <- NMDS_env$`Length Ice Cover (days)`
data.scores$Ice_off_DOY <- NMDS_env$Ice_off_DOY
data.scores$Ice_on_DOY <- NMDS_env$Ice_on_DOY
data.scores$Month <- NMDS_env$Month
data.scores$Month <- month.abb[data.scores$Month]
data.scores$Month <- factor(data.scores$Month, levels = c("Apr", "May", "Jun"))

species.scores <- as.data.frame(scores(cca_ice_phenology, display = "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores) 


p9 <- ggplot() + 
  theme_bw() + 
  geom_hline(yintercept=0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept=0, linetype = "dashed", color = "gray50") +
  geom_point(data=data.scores, mapping = aes(x=CCA1, y=CCA2, col = Ice_off_DOY, shape = Month), size=3) +
  scale_colour_gradient2(low = "#FDB863",
                         mid = "#B2ABD2",
                         high = "#542788",
                         midpoint = 115) + 
  labs(shape = NULL, colour = "Ice-off", x = "CCA1 (36.0 %)", y = "CCA2 (24.5 %)") +
  geom_segment(data = filter(veg_1, env != "Ice_on_DOY:Ice_off_DOY" & env != "Ice_on_DOY:Month" & env != "Ice_off_DOY:Month"), aes(x = 0, y = 0, xend = CCA1, yend = CCA2), color = "purple", arrow = arrow(length = unit(0.1, "cm"))) +
  geom_text_repel(data = filter(veg_1, env != "Ice_on_DOY:Ice_off_DOY" & env != "Ice_on_DOY:Month" & env != "Ice_off_DOY:Month"), aes(x = CCA1, y = CCA2, label = filter(veg_1, env != "Ice_on_DOY:Ice_off_DOY" & env != "Ice_on_DOY:Month" & env != "Ice_off_DOY:Month")$env), color = "purple", nudge_y = -0.05, size = 3)


p10 <- ggplot() + geom_hline(yintercept=0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept=0, linetype = "dashed", color = "gray50") +
  coord_equal() +
  lims(x = c(-3.3, 1.2)) +
  theme_bw() +
  geom_text(data = species.scores, aes(x = CCA1, y = CCA2, label = species.scores$species), size = 1.5) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18)) + 
  theme_bw() +
  labs (x = "CCA1 (36.0 %)", y = "") + 
  geom_segment(data = filter(veg_1, env != "Ice_on_DOY:Ice_off_DOY" & env != "Ice_on_DOY:Month" & env != "Ice_off_DOY:Month"), aes(x = 0, y = 0, xend = CCA1, yend = CCA2), color = "purple", arrow = arrow(length = unit(0.1, "cm"))) +
  geom_text_repel(data = filter(veg_1, env != "Ice_on_DOY:Ice_off_DOY" & env != "Ice_on_DOY:Month" & env != "Ice_off_DOY:Month"), aes(x = CCA1, y = CCA2, label = filter(veg_1, env != "Ice_on_DOY:Ice_off_DOY" & env != "Ice_on_DOY:Month" & env != "Ice_off_DOY:Month")$env), color = "purple", nudge_y = -0.05, size = 3)


p9 <- p9 +  
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position = "top") 

p10 <- p10 +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position = "top")

# p_final_snow_depth_NMDS <- (p9/p10) +
#   plot_annotation(tag_levels = list(c("")),
#                   tag_sep = '',
#                   tag_prefix = '',
#                   tag_suffix = '') &
#   theme(plot.tag.position = "topleft",
#         plot.tag = element_text(size = 8, hjust = 0, vjust = 0)) &
#   theme(plot.tag = element_text(face = 'bold'),
#         legend.position = "right")

p_final_ice_phenology_NMDS <- p9 + theme(plot.margin = unit(c(0,0,-1,0), "pt")) +
  p10 + theme(plot.margin = unit(c(-1,0,0,0), "pt")) + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "top")
p_final_ice_phenology_NMDS


#Write figures as pdf and tiff
setwd(fig_dir)
pdf("Figure_4.pdf", width=2*2.54, height=2*2.54)
p_final_ice_phenology_NMDS
dev.off()

setwd(fig_dir)
tiff("Figure_4.tiff", width=20, height = 10, unit = "cm", res = 600)
p_final_ice_phenology_NMDS 
dev.off()

#--------------------------------------------------------------------------------------------------------------
#                                 SUPPLEMENTARY FIGURE 1: PRECIPITATION
#------------------------------------------------------------------------------------------------------------
#Plot precipitation between the meteorological sites
p_Vallnora_precip <- ggplot() + 
  geom_point(filter(Vallnora_precip), mapping = aes(x = Rain_mm_Vallnora, y = Rain_mm_Svanberga), alpha = 0.5) + labs(y = "", x = "Vallnora prec. (mm)") + 
  theme_bw() + geom_abline(mapping = aes(intercept = 0, slope = 1), linetype = "dashed") + 
  stat_cor(Vallnora_precip, mapping = aes(x = Rain_mm_Svanberga, y = Rain_mm_Vallnora), cor.coef.name = "r", method="pearson") + 
  theme(axis.text.y = element_blank())

p_Norrveda_precip <- ggplot() + 
  geom_point(Norrveda_precip, mapping = aes(x = Rain_mm_Svanberga, y = Rain_mm_Norrveda), alpha = 0.5) + 
  labs(y = "Svanberga prec. (mm)", x = "Norrveda prec. (mm)") + 
  theme_bw() + geom_abline(mapping = aes(intercept = 0, slope = 1), linetype = "dashed") + 
  stat_cor(Norrveda_precip, mapping = aes(x = Rain_mm_Svanberga, y = Rain_mm_Norrveda), cor.coef.name = "r", method="pearson")

#Final plot
p_final_precip <- (p_Norrveda_precip|p_Vallnora_precip) & plot_annotation(tag_levels = 'a') &
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

m_1 <- cor.test(Snow_ice_data$Mean_snow, Snow_ice_data$Mean_ice, method = "pearson")
# #Extract model params
# p_val <- round(m_1[3]$p.value, 4) 
# rsqrd <- round(m_1[4]$estimate, 2)
# coeff <- as.data.frame(cbind(p_val, rsqrd))
# 
# # generate position of labs for plot
# labpos <- Snow_ice_data %>% 
#   ungroup() %>%
#   reframe(maxheight = max(Mean_snow), 
#           minheight = min(Mean_snow),
#           xpos = mean(Mean_ice))
# 
# #bind the position of labs and the lm stats to be plotted
# plotting_labs <- cbind(labpos, coeff)

p_snow_ice_cor <- ggplot(Snow_ice_data, aes(Mean_ice, y=Mean_snow)) + 
  geom_abline(intercept = 0 , slope = 0.7889886, linetype = "dashed") + 
  geom_point(size = 3) + theme_bw() + 
  labs(x = "Ice thickness (cm)", y = "Snow thickness (cm)") + 
  lims(x = c(0,NA),
       y = c(0,NA)) +
  # geom_text(plotting_labs,
  #           mapping = aes(x = 15, y =30, label = paste0("p = ", format(round(p_val, 3), scientific = F), ",\nR = ", round(rsqrd, 2))),
  #           family = "serif", fontface = "italic", size = 3) + 
  stat_cor(Snow_ice_data, mapping = aes(x = Mean_ice, y = Mean_snow), cor.coef.name = "R", method="pearson", p.accuracy = 0.001)


# Write plot into folder
setwd(fig_dir)
pdf("Figure_S2.pdf", width = 6, height = 2.5)
p_snow_ice_cor
dev.off()

setwd(fig_dir)
tiff("Figure_S2.tiff", unit = "cm", height = 5, width = 6, res = 600)
p_snow_ice_cor
dev.off()
setwd(dir)
#--------------------------------------------------------------------------------------------------------------
#                                 SUPPLEMENTARY FIGURE 3: ICE AND CHL-a
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

# generate position of labs for plot
labpos <- data %>% 
  ungroup() %>%
  reframe(maxheight = max(Mean_chla), 
          minheight = min(Mean_chla),
          xpos = mean(Mean_ice))

#Extract model params
p_val <- round(summary(m_1)$coefficients[,4][2]*2,2) #Multiply by 2 for Holms correction and round with three decimals
rsqrd <- round(summary(m_1)$adj.r.squared, 2)
coeff <- as.data.frame(cbind(p_val, rsqrd))

#bind the position of labs and the lm stats to be plotted
plotting_labs <- cbind(labpos, coeff)

p23 <- ggplot(data, aes(Mean_ice, y = Mean_chla)) +
  #geom_smooth(mapping = aes(x = Mean_ice, y = Mean_chla), method = "lm", linetype = "dashed", se = F, col = "black", linewidth = 0.3) + 
  geom_point(mapping = aes(x = Mean_ice, y = Mean_chla), size = 2) + 
  lims(x = c(0,55), y=c(0,NA)) + theme_bw() + ylab("January") +
  geom_text(plotting_labs,
            mapping = aes(x = 45, y = minheight + 0.68 * (maxheight - minheight), label = paste0("p = ", format(round(p_val, 4), scientific = F), ",\nR² = ", round(rsqrd, 2))),
            family = "serif", fontface = "italic", size = 2.5)


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

# generate position of labs for plot
labpos <- data %>% 
  ungroup() %>%
  reframe(maxheight = max(Mean_chla), 
          minheight = min(Mean_chla),
          xpos = mean(Mean_ice))

#Extract model params
p_val <- round(summary(m_2)$coefficients[,4][2]*1, 2) #Multiply by 1 for Holms correction and round with three decimals
rsqrd <- round(summary(m_2)$adj.r.squared, 2)
coeff <- as.data.frame(cbind(p_val, rsqrd))

#bind the position of labs and the lm stats to be plotted
plotting_labs <- cbind(labpos, coeff)

#Plot the data
p <- ggplot(data, aes(Mean_ice, Mean_chla))
p <- p + #geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), se = F, color = "black", linewidth = 0.3, linetype = "dashed") + 
  geom_point(size = 2)
p24 <- p + lims(x=c(0,55)) + scale_y_continuous(name = "February") + theme_bw() +   
  geom_text(plotting_labs, mapping = aes(x = 45, y = minheight + 0.71 * (maxheight - minheight), label = paste0("p = ", format(round(p_val, 4), scientific = F), ",\nR² = ", round(rsqrd, 2))), family = "serif", fontface = "italic", size = 2.5)

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

# generate position of labs for plot
labpos <- data %>% 
  ungroup() %>%
  reframe(maxheight = max(Mean_chla), 
          minheight = min(Mean_chla),
          xpos = mean(Mean_ice))

#Extract model params
p_val <- round(summary(m_3)$coefficients[,4][2]*3, 3) #Multiply by 3 for Holms correction and round with three decimals
rsqrd <- round(summary(m_3)$adj.r.squared, 2)
coeff <- as.data.frame(cbind(p_val, rsqrd))

#bind the position of labs and the lm stats to be plotted
plotting_labs <- cbind(labpos, coeff)

#Plot everything
p <- ggplot(data, aes(Mean_ice, Mean_chla))
p <- p + geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), se = F, color = "black", linewidth = 0.3, linetype = "dashed") +
  geom_point(size = 2)
p25 <- p + #geom_function(fun = function(x) (24.905*exp(-0.07*x)), colour = "black", linetype=2) +  
  scale_y_continuous(trans = "log10", name = "March") + #,
  #sec.axis = sec_axis(trans=~.*1)) + 
  theme_bw() +
  annotation_logticks(sides = "l") + lims(x=c(0,55)) +
  geom_text(plotting_labs,
            mapping = aes(x = 45, y = minheight + 0.35 * (maxheight - minheight), label = paste0("p = ", format(round(p_val, 4), scientific = F), ",\nR² = ", round(rsqrd, 2))),
            family = "serif", fontface = "italic", size = 2.5)


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
p_final_production_ice_thickness
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

Zooplankton_nutrient_light$Month_abbrev <- factor(Zooplankton_nutrient_light$Month_abbrev, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun"))
#Plot all variables
p_zoops <- ggplot(filter(Zooplankton_nutrient_light, Taxa == "Zooplankton"), mapping = aes(x = Month_abbrev, y = Biomass)) +
  geom_boxplot(fill = "white", outlier.shape = NA) +
  geom_point(size = 1) + 
  theme_bw() + theme(strip.text = element_text(face = "bold"),
                     axis.text = element_text(size = 8),
                     axis.title = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) +
  facet_wrap(~Taxa, scales = "free_y") + labs(x = "", y = expression("Biovolume (mm"^3*" L"^-1*")")) + 
  scale_x_discrete(labels=c("J", "F", "M", "A", "M", "J"))

p_phyto <- ggplot(filter(Zooplankton_nutrient_light, Taxa == "Phytoplankton"), mapping = aes(x = Month_abbrev, y = Biomass)) +
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
  geom_boxplot(filter(Zooplankton_nutrient_light, Chemistry == "TN" | Chemistry == "TP" | Chemistry == "Si"), 
               mapping = aes(x = Month_abbrev, y = Concentration),
               fill = "white", outlier.shape = NA) + 
  theme_bw() + theme(strip.text = element_text(face = "bold"),
                     axis.text = element_text(size = 8),
                     axis.title = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) + geom_point(filter(Zooplankton_nutrient_light, Chemistry == "TN" | Chemistry == "TP" |
                                                                                  Chemistry == "Si"), 
                                                                         mapping = aes(x = Month_abbrev, y = Concentration), size = 1) +
  facet_wrap(~Chemistry, scales = "free_y") + 
  labs(x="", y = expression(paste("mg "*L^-1))) + 
  scale_x_discrete(labels=c("J", "F", "M", "A", "M", "J"))

p_temp <- ggplot() +
  geom_boxplot(filter(Zooplankton_nutrient_light, Chemistry == "Temperature"), 
               mapping = aes(x = Month_abbrev, y = Concentration), fill = "white", outlier.shape = NA) +
  geom_point(filter(Zooplankton_nutrient_light, Chemistry == "Temperature"), 
             mapping = aes(x = Month_abbrev, y = Concentration), size = 1) + 
  theme_bw() + theme(strip.text = element_text(face = "bold"),
                     axis.text = element_text(size = 8),
                     axis.title = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) +
  facet_wrap(~Chemistry, scales = "free_y") + labs(y = expression(""*~degree*C*""), x = "") + 
  scale_x_discrete(labels=c("J", "F", "M", "A", "M", "J"))

p_light <- ggplot() +
  geom_boxplot(filter(Zooplankton_nutrient_light, Chemistry == "SW radiation"), 
               mapping = aes(x = Month_abbrev, y = Concentration), fill = "white", outlier.shape = NA) +
  geom_point(filter(Zooplankton_nutrient_light, Chemistry == "SW radiation"), 
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
Regression_data <- Phytoplankton_after_ice_out
#Calculate dominant taxa
Dominant_taxa <- Regression_data %>% filter(Month == 4) %>% group_by(Genus, Month) %>% reframe(mean_biomass = mean(mean_biomass)) %>% arrange(-mean_biomass)
#Replace NA with 0 to calculate mean
Regression_data$mean_biomass[is.na(Regression_data$mean_biomass)] <- 0
#Count each occurrence of taxa to check that they are the same no. to start with
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
  select(Month, Genus, Winter_year,Ice_on_DOY, Ice_off_DOY, `Length Ice Cover (days)`, n, mean_biomass)

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

correlation_results_iceon <- correlation_results %>% filter(Ice_phenology == "Ice_on_DOY" & Month == 4)
correlation_results_iceoff <- correlation_results %>% filter(Ice_phenology == "Ice_off_DOY" & Month == 4)

p_adj_ice_off <- p.adjust(correlation_results_iceoff$p.value, method = "holm")
correlation_results_iceoff <- cbind(correlation_results_iceoff, p_adj_ice_off)
p_adj_ice_on <- p.adjust(correlation_results_iceon$p.value, method = "holm")
correlation_results_iceon <- cbind(correlation_results_iceon, p_adj_ice_on)

correlation_results <- rbind(correlation_results_iceon, correlation_results_iceoff)

Data_for_plotting <- left_join(correlation_results, Regression_data)

Data_for_plotting$Ice_phenology <- factor(Data_for_plotting$Ice_phenology, levels = c("Length Ice Cover (days)", "Ice_on_DOY", "Ice_off_DOY"))
Data_for_plotting <- Data_for_plotting %>% group_by(Genus, Month, Ice_phenology) %>% mutate(Max_biomass = max(mean_biomass),
                                                                                            Max_DOY = max(DOY_days),
                                                                                            Min_DOY = min(DOY_days))
library(ggpubr)
Data_for_plotting <- Data_for_plotting %>% group_by(Month, Ice_phenology, Genus) %>% arrange(estimate)
Data_for_plotting <- Data_for_plotting %>% mutate(p_adj = case_when(Genus == "Aulacoseira" & Ice_phenology == "Ice_on_DOY"~ 0.760,
                                                                    Genus == "Chlamydomonas" & Ice_phenology == "Ice_on_DOY"~0.872, T~1))
coefficients_all <- Data_for_plotting %>% group_by(Month, Ice_phenology, Genus) %>% reframe(p_adj = mean(...12),
                                                                                            estimate = mean(estimate))

# generate position of labs for plot
labpos <- Data_for_plotting %>% 
  ungroup() %>%
  group_by(Month, Ice_phenology, Genus) %>%
  reframe(maxheight = max(mean_biomass), 
          minheight = min(mean_biomass),
          xpos = mean(DOY_days)) %>%
  mutate(Month = as.numeric(Month))

#bind the position of labs and the lm stats to be plotted
plotting_labs <- left_join(labpos, coefficients_all)
#----------------------------------------------------PLOTTING
#Remove heterotrophs from plots
Data_for_plotting <- left_join(Data_for_plotting, plotting_labs)
Data_for_plotting <- Data_for_plotting %>% filter(Month == 4)

Taxa_regressions <- ggplot(filter(Data_for_plotting, Ice_phenology != "Length Ice Cover (days)" & Month == 4), mapping = aes(x = DOY_days, y = mean_biomass, col = Ice_phenology, shape = Ice_phenology)) + 
  facet_wrap(~Genus, scales = "free_y", ncol = 4, nrow = 7) + lims(y = c(0,NA)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("#000099", "#660000"), labels = c("Ice-on", "Ice-off")) + 
  scale_shape_manual(values = c(17,19), labels = c("Ice-on", "Ice-off")) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"),
        legend.text = element_text(size = 12)) +
  labs(x = "Day-of-the-year (DOY)", y = expression("Biomass ("*mu*"g l"^-1*")"), col = "", shape = "", size = 14) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  geom_text(Data_for_plotting, mapping = aes(x = Max_DOY*0.8, y = Max_biomass* 0.8, label = paste0("p = ", format(round(p_adj, 2), scientific = F), ",\nr = ", round(estimate, 3))), family = "serif", fontface = "italic", size = 2.5)


setwd(fig_dir)
tiff("Figure_S5.tiff", width = 25, height = 25, units = "cm", res = 1200)
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
#                                 SUPPLEMENTARY FIGURE 6: SNOW AND ICE CCA
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
df_filter <- df
df_filter <- df_filter %>% select(-Other)
df_filter <- df_filter[,colSums(!is.na(df_filter)) > 0*nrow(df_filter)] #Filter data based on percentage of non-NA cases (not used)
Phyto_snow <- hybrd.rplc_if(df_filter)
Phyto_snow_filter <- Phyto_snow[rowSums(Phyto_snow[,c(5:ncol(Phyto_snow))])>0,] #Remove samplings without any data
NMDS_spec <- Phyto_snow_filter[, c(5:ncol(Phyto_snow_filter))]
NMDS_env <- Phyto_snow_filter[,c(2:3,4)]
#sp.hel <- decostand(NMDS_spec, method="hellinger") #tranform data for unconstrained

library(ggrepel)
# CCA 
cca_env <- NMDS_env %>% mutate(Month = as.numeric(Month))
cca_snow_depth <- cca(NMDS_spec ~ Mean_snow*Month, data=cca_env)
cca_snow_depth_anova <- anova.cca(cca_snow_depth, by = "term", permutations = 9999)
cca_snow_depth_anova_model <- anova.cca(cca_snow_depth, permutations = 9999)
#CCA1 = 0.6465045 ~ 64.7%
#CCA2 = 0.2076479 ~ 20.8%
#extracting the data as data frame; env data
veg_1_snow <- as.data.frame(cca_snow_depth$CCA$biplot)
veg_1_snow["env"] <- row.names(veg_1_snow)

#extracting the data; genusv
veg_2_snow <- as.data.frame(cca_snow_depth$CCA$v)
veg_2_snow["genus"] <- row.names(veg_2_snow)
veg_1_snow <- veg_1_snow %>% mutate(env = case_when(env == "Mean_snow"~"Snow", T~env))

# PLOTTING
data.scores_snow <- as.data.frame(scores(cca_snow_depth)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores_snow$site <- rownames(data.scores_snow)  # create a column of site names, from the rownames of data.scores
data.scores_snow$Snow <- NMDS_env$Mean_snow
data.scores_snow$Month <- NMDS_env$Month
data.scores_snow$Month <- month.abb[data.scores_snow$Month]
data.scores_snow$Month <- factor(data.scores_snow$Month, levels = c("Jan", "Feb", "Mar"))

species.scores_snow <- as.data.frame(scores(cca_snow_depth, display = "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores_snow$species <- rownames(species.scores_snow)  # create a column of species, from the rownames of species.scores
head(species.scores) 


p11 <- ggplot() + 
  theme_bw() + 
  geom_hline(yintercept=0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept=0, linetype = "dashed", color = "gray50") +
  geom_point(data=data.scores_snow, mapping = aes(x=CCA1, y=CCA2, col = Snow, shape = Month), size=4) +
  guides(shape  = guide_legend(order = 1)) +
  scale_colour_gradient2(low = "#FDB863",
                         mid = "#B2ABD2",
                         high = "#542788",
                         midpoint = 20) + 
  labs(shape = NULL, colour = "Snow depth", x = "", y = "CCA2 (20.8 %)") +
  geom_segment(data = filter(veg_1_snow, env != "Mean_snow:Month"), aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel(data = filter(veg_1_snow, env != "Mean_snow:Month"), aes(x = CCA1, y = CCA2, label = filter(veg_1_snow, env != "Mean_snow:Month")$env), nudge_y = -0.05, size = 3)


p12 <- ggplot() + 
  geom_hline(yintercept=0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept=0, linetype = "dashed", color = "gray50") +
  coord_equal() +
  theme_bw() +
  geom_text(data = species.scores_snow, aes(x = CCA1, y = CCA2, label = species.scores_snow$species), size = 2.5) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18)) + 
  theme_bw() +
  labs (x = "CCA1 (64.7 %)", y = "CCA2 (20.8 %)") + 
  geom_segment(data = filter(veg_1_snow, env != "Mean_snow:Month"), aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel(data = filter(veg_1_snow, env != "Mean_snow:Month"), aes(x = CCA1, y = CCA2, label = filter(veg_1_snow, env != "Mean_snow:Month")$env), nudge_y = -0.05, size = 3) +
  lims(x = c(-1.6, 1.8))


p11 <- p11 +  
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "right") 

p12 <- p12 +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "right")



#ICE THICKNESS CCA
#Join phytoplankton data and ice data
Mean_phytoplankton_genera_filter <- Mean_phytoplankton_genera %>% select(-Functional_group)
Phyto_ice <- left_join(Mean_phytoplankton_genera_filter, Mean_ice, by = c("Winter_year", "Month"))

#Remove NA of ice thickness
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

# #Calculate percentage of unknown taxa
df <- df %>% spread(Taxa, Mean_biomass) %>% drop_na(Mean_ice)
Phyto_ice <- hybrd.rplc_if(df)
df_filter <- df
df_filter <- df_filter %>% select(-Other)
df_filter <- df_filter[,colSums(!is.na(df_filter)) > 0*nrow(df_filter)] #Filter data based on percentage of non-NA cases (>38% data)
Phyto_ice <- hybrd.rplc_if(df_filter)
Phyto_ice_filter <- Phyto_ice[rowSums(Phyto_ice[,c(5:ncol(Phyto_ice))])>0,] #Remove samplings without any data
NMDS_spec <- Phyto_ice_filter[, c(5:ncol(Phyto_ice_filter))]
NMDS_env <- Phyto_ice_filter[,c(2:3,4)]
#sp.hel <- decostand(NMDS_spec, method="hellinger")

# CCA 
cca_env <- NMDS_env %>% mutate(Month = as.numeric(Month))
cca_ice_thickness <- cca(NMDS_spec ~ Mean_ice*Month, data=cca_env)
cca_ice_thickness_anova <- anova.cca(cca_ice_thickness, by = "term", permutations = 9999)
cca_ice_thickness_anova_model <- anova.cca(cca_ice_thickness, permutations = 9999)

#extracting the data as data frame; env data
veg_1 <- as.data.frame(cca_ice_thickness$CCA$biplot)
veg_1["env"] <- row.names(veg_1)

#extracting the data; genusv
veg_2 <- as.data.frame(cca_ice_thickness$CCA$v)
veg_2["genus"] <- row.names(veg_2)
veg_1 <- veg_1 %>% mutate(env = case_when(env == "Mean_ice"~"Ice", T~env))

# PLOTTING
data.scores <- as.data.frame(scores(cca_ice_thickness)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$Ice <- NMDS_env$Mean_ice
data.scores$Month <- NMDS_env$Month
data.scores$Month <- month.abb[data.scores$Month]
data.scores$Month <- factor(data.scores$Month, levels = c("Jan", "Feb", "Mar"))

species.scores <- as.data.frame(scores(cca_ice_thickness, display = "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores) 


p13 <- ggplot() + 
  theme_bw() + 
  geom_hline(yintercept=0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept=0, linetype = "dashed", color = "gray50") +
  geom_point(data=data.scores, mapping = aes(x=CCA1, y=CCA2, col = Ice, shape = Month), size=4) +
  guides(shape  = guide_legend(order = 1)) +
  scale_colour_gradient2(low = "#FDB863",
                         mid = "#B2ABD2",
                         high = "#542788",
                         midpoint = 35) + 
  labs(shape = NULL, colour = "Ice thickness", x = "", y = "CCA2 (39.2 %)") +
  geom_segment(data = filter(veg_1, env != "Mean_ice:Month"), aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel(data = filter(veg_1, env != "Mean_ice:Month"), aes(x = CCA1, y = CCA2, label = filter(veg_1, env != "Mean_ice:Month")$env), nudge_y = -0.05, size = 3)


p14 <- ggplot() + geom_hline(yintercept=0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept=0, linetype = "dashed", color = "gray50") +
  coord_equal() +
  theme_bw() +
  geom_text(data = species.scores, aes(x = CCA1, y = CCA2, label = species.scores$species), size = 2.5) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18)) + 
  theme_bw() +
  labs (x = "CCA1 (53.5 %)", y = "CCA2 (39.2 %)") + 
  geom_segment(data = filter(veg_1, env != "Mean_ice:Month"), aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel(data = filter(veg_1, env != "Mean_ice:Month"), aes(x = CCA1, y = CCA2, label = filter(veg_1, env != "Mean_ice:Month")$env), nudge_y = -0.05, size = 3) + 
  lims(x = c(-1.9,1.95))


p13 <- p13 +  
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "left")

p14 <- p14 +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "left")


# Combine plots
library(patchwork)

p_cca_snow_depth <- p11 +
  p12 + 
  plot_layout(ncol = 1, guides = "collect") & 
  theme(legend.position = "left") + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
p_cca_snow_depth <- wrap_elements(p_cca_snow_depth)

p_cca_ice_thickness <- p13 + 
  p14 +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "right") +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
p_cca_ice_thickness <- wrap_elements(p_cca_ice_thickness)

p_ice_snow_cca_suppl <- (p_cca_snow_depth+p_cca_ice_thickness) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = c('a'))
#p_ice_snow_cca_suppl <- wrap_elements(p_ice_snow_cca_suppl)

setwd(fig_dir)
tiff("Figure_S6.tiff", units="cm", width=30, height=20, res=1200)
p_ice_snow_cca_suppl
dev.off()
setwd(dir)


#--------------------------------------------------------------------------------------------------------------
#                                 SUPPLEMENTARY FIGURE 7: Chl-a: Biomass
#------------------------------------------------------------------------------------------------------------
# Replace Nas with 0s
Phytoplankton_below_ice <- hybrd.rplc_if(Phytoplankton_below_ice)

#Sum all phytoplankton biomass per sampling day across all taxa (based on phylum)
Phytoplankton_data_sum_day <- Phytoplankton_below_ice %>% group_by(Day, Winter_year, Month) %>% reframe(sum_biomass = sum(biomass_ug_l))

Chla_snow <- left_join(Chla_below_ice, Snow_on_ice, by = c("Day", "Month", "Year", "Winter_year"))
Phyto_snow <- left_join(Phytoplankton_data_sum_day, Snow_on_ice, by = c("Day", "Winter_year", "Month"))
Snow_ratio_biomass_chla_data <- left_join(Chla_snow, Phyto_snow)

#Calculate mean
Snow_ratio_biomass_chla_data <- Snow_ratio_biomass_chla_data %>%
  drop_na(sum_biomass, `Chl a_µg/l`, Snow_depth) %>%
  mutate(Snow_depth = Snow_depth*100) %>%
  filter(Month == 1 | Month == 2 | Month == 3) %>%
  group_by(Winter_year, Month) %>%
  reframe(
    N_tot = n(),
    Snow_mean = mean(Snow_depth),
    Snow_sd = sd(Snow_depth),
    Biomass_mean = mean(sum_biomass),
    Biomass_sd = sd(sum_biomass),
    Chla_mean = mean(`Chl a_µg/l`),
    Chla_sd = sd(`Chl a_µg/l`))

Snow_ratio_biomass_chla_data$Ratio <- Snow_ratio_biomass_chla_data$Chla_mean/Snow_ratio_biomass_chla_data$Biomass_mean
Snow_ratio_biomass_chla_data$Month <- as.numeric(Snow_ratio_biomass_chla_data$Month)
Snow_ratio_biomass_chla_data$Month_abbrev <- month.abb[Snow_ratio_biomass_chla_data$Month] 
Snow_ratio_biomass_chla_data$Month_abbrev <- factor(Snow_ratio_biomass_chla_data$Month_abbrev, levels = c("Jan", "Feb", "Mar"))
Snow_ratio_biomass_chla_data$Winter_year <- as.numeric(Snow_ratio_biomass_chla_data$Winter_year)

lmJan <- lm(Ratio ~ Snow_mean, filter(Snow_ratio_biomass_chla_data, Month == 1))
summary(lmJan)
lmFeb <- lm(Ratio ~ Snow_mean, filter(Snow_ratio_biomass_chla_data, Month == 2))
summary(lmFeb)
lmMar <- lm(Ratio ~ Snow_mean, filter(Snow_ratio_biomass_chla_data, Month == 3))
summary(lmMar)

#Extract model params
p_val <- round((summary(lmJan)$coefficients[,4][2]), 2) #Multiply by 2 for Holms correction and round with three decimals
rsqrd <- round(summary(lmJan)$adj.r.squared, 2)
coeff1 <- as.data.frame(cbind(p_val, rsqrd))

#Extract model params
p_val <- round((summary(lmFeb)$coefficients[,4][2]), 2) #Multiply by 2 for Holms correction and round with three decimals
rsqrd <- round(summary(lmFeb)$adj.r.squared, 2)
coeff2 <- as.data.frame(cbind(p_val, rsqrd))

#Extract model params
p_val <- round((summary(lmMar)$coefficients[,4][2]), 2) #Multiply by 2 for Holms correction and round with three decimals
rsqrd <- round(summary(lmMar)$adj.r.squared, 2)
coeff3 <- as.data.frame(cbind(p_val, rsqrd))

coeffs <- rbind(coeff1, coeff2, coeff3)
Month_abbrev <- rbind("Jan", "Feb", "Mar")
coeffs <- cbind(coeffs, Month_abbrev)

# generate position of labs for plot
labpos <- Snow_ratio_biomass_chla_data %>% 
  ungroup() %>%
  reframe(maxheight = max(Ratio), 
          minheight = min(Ratio),
          xpos = mean(Snow_mean))

#bind the position of labs and the lm stats to be plotted
plotting_labs <- cbind(labpos, coeffs)
Snow_ratio_biomass_chla_data$Month_abbrev <- factor(Snow_ratio_biomass_chla_data$Month_abbrev, levels = c("Jan", "Feb", "Mar"))
plotting_labs$Month_abbrev <- factor(plotting_labs$Month_abbrev, levels = c("Jan", "Feb", "Mar"))
#plot the chla conc. per biomass conc.
p_chla_per_biomass <- ggplot(Snow_ratio_biomass_chla_data, aes(x=Snow_mean, y = Ratio)) + geom_point(size = 3) + 
  facet_wrap(~Month_abbrev) + 
  theme_bw() + 
  labs(x = "Snow depth (cm)", y = expression(paste("Chl-"*alpha*"/Biomass"))) +
  geom_text(plotting_labs,
            mapping = aes(x = 20, y = minheight + 0.85 * (maxheight - minheight), label = paste0("p = ", format(round(p_val, 4), scientific = F), ",\nR² = ", round(rsqrd, 2))),
            family = "serif", fontface = "italic", size = 2.5) + 
  theme(strip.text = element_text(face = "bold"))



#Write figure into fig dir.
setwd(fig_dir)
tiff("Figure_S7.tiff", units="cm", width=15, height=5, res=500)
p_chla_per_biomass
dev.off()

setwd(fig_dir)
pdf("Figure_S7.pdf", width=4.5, height=2)
p_chla_per_biomass
dev.off()
setwd(dir)

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
Dominant_taxa <- Phytoplankton_after_ice_out %>% filter(Month_abb == "Apr") %>% group_by(Genus, Month_abb) %>% reframe(Mean_phyto_biomass = mean(mean_biomass))

# Calculate unknown taxa in biomass after ice-off data (Apr-Jun)
#Replace NA with 0
Phytoplankton_after_ice_out <- hybrd.rplc_if(Phytoplankton_after_ice_out)
#Order months
Phytoplankton_after_ice_out$Month_abb <- factor(Phytoplankton_after_ice_out$Month_abb, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
#Calculate percentage
Known_genera <- Phytoplankton_after_ice_out %>% filter(Genus != "Other") %>% reframe(Biomass = sum(mean_biomass))
Unknown_genera <- Phytoplankton_after_ice_out %>% filter(Genus == "Other") %>% reframe(Biomass = sum(mean_biomass))
Perc <- (Unknown_genera/(Known_genera+Unknown_genera))*100

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
Snow_data <- separate(Snow_on_ice, col = "Date", into = c("Year", "Month", "Day"))
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

