# Data preparation
library(tidyverse)
library(readxl)
library(epitools)

# New updated data with meningitis (25.04.2024)
# All df are stored in raw_data
dat <- readxl::read_excel("raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_withdeath_meningitis_clean.xlsx") #%>% 
# glimpse()

dat <- dat %>% 
  dplyr::rename(Earliest.specimen.date = Earliestspecimendate,
         current.region.name = currentregionname)

dat_G <- dat %>% 
  dplyr::mutate(AGEYR = ifelse(AGEYR >= 90, 90, as.numeric(AGEYR)), # For incidence calculation, data grouped for people aged 90+
         year = year(Earliest.specimen.date),
         month = month(Earliest.specimen.date),
         vacc = case_when(
           year < 2006 ~ "Pre-PCV7",
           year >= 2006 & year < 2011 ~ "PCV7",
           year >= 2011 ~ "PCV13",
           TRUE ~ NA_character_
         ),
         ageGroup = case_when( # edit 5 age bands
           AGEYR < 5 ~ "<5",
           AGEYR >= 5 & AGEYR < 19 ~ "5-18",
           AGEYR >= 19 & AGEYR < 31 ~ "19-30",
           AGEYR >= 31 & AGEYR < 65 ~ "31-64",
           AGEYR >= 65 ~ "65+",
           is.na(AGEYR) ~ "Unknown" # 16 IDs have no AGEYR
           # TRUE ~ "Unknown" 
         ),
         ageGroup2 = case_when(
           AGEYR < 15 ~ "children",
           AGEYR >= 15 ~ "adults",
           is.na(AGEYR) ~ "Unknown" # 16 IDs have no AGEYR
         ),
         current.region.name = ifelse(current.region.name == "EASTERN", "EAST", current.region.name), # Wrong perception of "EASTERN" that should be "EAST"
         current.region.name = case_when(
           current.region.name == "E MIDS" ~ "East Midlands",
           current.region.name == "EAST" ~ "East of England",
           current.region.name == "LONDON" ~ "London",
           current.region.name == "N EAST" ~ "North East",
           current.region.name == "N WEST" ~ "North West",
           current.region.name == "S EAST" ~ "South East",
           current.region.name == "S WEST" ~ "South West",
           current.region.name == "W MIDS" ~ "West Midlands",
           current.region.name == "YORK&HUM" ~ "Yorkshire and The Humber",
           TRUE ~ current.region.name
         ),
         ageLabel = ifelse(AGEYR >= 90, 90, as.numeric(AGEYR)), # For incidence calculation, data grouped for people aged 90+
  ) #%>% 
  # glimpse()

# EpiDescription based on incidences and CI
# Total population data by age, year for each region
# SOURCE: https://www.nomisweb.co.uk/
pop <- readxl::read_excel("raw_data/nomis_2024_04_15_124553_DCedit.xlsx") #%>% 
  # glimpse()

pop_l <- pop %>% 
  tidyr::pivot_longer(cols = `2001`:`2022`,
               names_to = "Year",
               values_to = "PopSize") %>% 
  dplyr::mutate(Age = gsub("Age ", "", Age),
         Age = ifelse(Age == "Aged 90+", 90, as.numeric(Age)), # For incidence calculation, data grouped for people aged 90+
         ageGroup = case_when( # edit 5 age bands
           Age < 5 ~ "<5",
           Age >= 5 & Age < 19 ~ "5-18",
           Age >= 19 & Age < 31 ~ "19-30",
           Age >= 31 & Age < 65 ~ "31-64",
           Age >= 65 ~ "65+",
           is.na(Age) ~ "Unknown" # 16 IDs have no Age
           # TRUE ~ "Unknown" 
         ),
         ageGroup2 = case_when(
           Age < 15 ~ "children",
           Age >= 15 ~ "adults",
           is.na(Age) ~ "Unknown" # 16 IDs have no AGEYR
         ),
         Year = as.numeric(Year)) %>% 
  glimpse()

# Vaccination programme:
# SOURCE: https://www.gov.uk/government/publications/pneumococcal-the-green-book-chapter-25
vaccine_UK <- data.frame(
  year = c(2006, 2011),
  vaccine = c("PCV7", "PCV13")
)
# Colour names:
# https://www.datanovia.com/en/blog/awesome-list-of-657-r-color-names/
col_map <- c("<5" = "indianred4",
             "5-18" = "orange",
             "19-30" = "seagreen4",
             "31-64" = "steelblue",
             "65+" = "purple3",
             "Unknown" = "black",
             "children" = "darkred",
             "adults" = "darkblue"
)

vacc_map <- c("PCV7" = "gray80",
              "PCV13" = "gray20")

col_imD <- c(incid_Ser1 = "deepskyblue3",
             incid_m = "green",
             incid_D = "maroon")

# CI calculations for children-adults
ageGroup2 <- dat_G %>% 
  dplyr::group_by(year, ageGroup2) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup()

pop_ageGroup2 <- pop_l %>% 
  dplyr::group_by(Year, ageGroup2) %>% 
  dplyr::summarise(PopSize = sum(PopSize)) %>% 
  dplyr::ungroup()

all_ageGroup2 <- merge(ageGroup2, pop_ageGroup2,
               by.x = c("year","ageGroup2"),
               by.y = c("Year", "ageGroup2")) %>%
  dplyr::mutate(Conf_Int = epitools::binom.exact(counts, PopSize),
         incid_Ser1 = Conf_Int$proportion) # per-100,000 population

write.csv(all_ageGroup2, "raw_data/incidence_CI_per_year_2_ageGroup.csv", row.names = FALSE)

# Viz counts
png("pictures/counts_2ageGroups.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(all_ageGroup2, aes(x = year, y = counts, group = ageGroup2,
                color = ageGroup2)) +
  geom_line(size = 1.5) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("children", "adults"),
                     labels = c("Children (< 15)", "Adults")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  geom_label(aes(x = 2006, y = 150, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 150, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("The Counts of Serotype 1 in England by Demographic Groups") +
  xlab("Year") +
  ylab("Serotype 1 Cases")
dev.off()

# Viz incidence
png("pictures/incidence_2ageGroups.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(all_ageGroup2, aes(x = year, y = Conf_Int$proportion*100000, group = ageGroup2,
                  color = ageGroup2)) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = Conf_Int$lower*100000, ymax = Conf_Int$upper*100000), # It doesn't matter whether I add the CI or not because the Pop data is quite huge, I suppose (?)
                width = .1) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("children", "adults"),
                     labels = c("Children (< 15)", "Adults")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  scale_linetype_manual(values = c(vacc_map),
                        name = "Vaccine",
                        labels = c("PCV7", "PCV13")) +
  geom_label(aes(x = 2006, y = 0.15, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 0.15, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("The Incidence of Serotype 1 in England \nby Demographic Groups (per 100,000)") +
  xlab("Year") +
  ylab("Serotype 1 Incidence")
dev.off()

# CI calculations for 5 ageGroups
ageGroup5 <- dat_G %>% 
  dplyr::group_by(year, ageGroup) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup()

pop_ageGroup5 <- pop_l %>% 
  dplyr::group_by(Year, ageGroup) %>% 
  dplyr::summarise(PopSize = sum(PopSize)) %>% 
  dplyr::ungroup()

all_ageGroup5 <- merge(ageGroup5, pop_ageGroup5,
                       by.x = c("year","ageGroup"),
                       by.y = c("Year", "ageGroup")) %>%
  dplyr::mutate(Conf_Int = epitools::binom.exact(counts, PopSize),
                incid_Ser1 = Conf_Int$proportion) # per-100,000 population

write.csv(all_ageGroup5, "raw_data/incidence_CI_per_year_5_ageGroup.csv", row.names = FALSE)

# Viz counts
png("pictures/counts_5ageGroups.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(all_ageGroup5, aes(x = year, y = counts, group = ageGroup,
                          color = ageGroup)) +
  geom_line(size = 1.5) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("<5", "5-18", "19-30", "31-64", "65+", "Unknown"),
                     labels = c("<5", "5-18", "19-30", "31-64", "65+", "Unknown")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  geom_label(aes(x = 2006, y = 175, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 250, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("The Counts of Serotype 1 in England by Demographic Groups") +
  xlab("Year") +
  ylab("Serotype 1 Cases")
dev.off()

# Viz incidence
png("pictures/incidence_5ageGroups.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(all_ageGroup5, aes(x = year, y = Conf_Int$proportion*100000, group = ageGroup,
                  color = ageGroup)) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = Conf_Int$lower*100000, ymax = Conf_Int$upper*100000), # It doesn't matter whether I add the CI or not because the Pop data is quite huge, I suppose (?)
                width = .1) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("<5", "5-18", "19-30", "31-64", "65+", "Unknown"),
                     labels = c("<5", "5-18", "19-30", "31-64", "65+", "Unknown")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  scale_linetype_manual(values = c(vacc_map),
                        name = "Vaccine",
                        labels = c("PCV7", "PCV13")) +
  geom_label(aes(x = 2006, y = 0.15, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 0.15, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("The Incidence of Serotype 1 in England \nby Demographic Groups (per 100,000)") +
  xlab("Year") +
  ylab("Serotype 1 Incidence")
dev.off()

# Basic case count data without age structure or regions
# Create all hypothetical recorded disease date
dat_G$Earliest.specimen.date <- as.Date(dat_G$Earliest.specimen.date)
all_date <- data.frame(allDate = seq.Date(from = min(dat_G$Earliest.specimen.date),
                                          to = max(dat_G$Earliest.specimen.date), 
                                          by = 1))
all_date$day <- 1:nrow(all_date)

# Coz the incidence only requires 2 columns called "counts" and "Day" in NUMBERS
# The counts (but in 0 counts the date are not recorded)
Natm_ni <- dat_G %>% 
  dplyr::group_by(Earliest.specimen.date) %>% 
  dplyr::summarise(counts_Ser1 = n()) %>% 
  dplyr::ungroup() #%>% 
# glimpse()

Natm_nmeningitis <- dat_G %>% 
  dplyr::filter(MeningitisFlag == "Y") %>% 
  dplyr::group_by(Earliest.specimen.date) %>% 
  dplyr::summarise(counts_meningitis = n()) %>% 
  dplyr::ungroup() #%>% 
# glimpse()

Natm_n30DDeath <- dat_G %>% 
  dplyr::filter(`30daydeath` == "D") %>% 
  dplyr::group_by(Earliest.specimen.date) %>% 
  dplyr::summarise(counts_30DDeath = n()) %>% 
  dplyr::ungroup() #%>% 
# glimpse()


# Create a new df based on counts per day for Serotype 1, meningitis, and 30 days death
Natm_n_i <- dplyr::full_join(all_date, Natm_ni,
                      by = c("allDate" = "Earliest.specimen.date"))

Natm_n_im <- dplyr::full_join(Natm_n_i, Natm_nmeningitis,
                       by = c("allDate" = "Earliest.specimen.date"))

Natm_n_imD <- dplyr::full_join(Natm_n_im, Natm_n30DDeath,
                        by = c("allDate" = "Earliest.specimen.date")) %>% 
  replace(is.na(.), 0) #%>% # NA means no data of meningitis or 30 days death, changed them to 0
  # glimpse()


# Examples on https://github.com/mrc-ide/mcstate/blob/master/inst/sir_incidence.csv
# Requires case count per aligned day only

# Viz per-day counts by base R plot
png("pictures/daily_cases.png", width = 17, height = 12, unit = "cm", res = 1200)
par(bty = "n", mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0))
col_imD <- c(counts_Ser1 = "deepskyblue3",
             counts_meningitis = "green",
             counts_30DDeath = "maroon")
plot(Natm_n_imD$allDate, Natm_n_imD$counts_Ser1, type = "b",
     xlab = "Date (year)", ylab = "Counts",
     ylim = c(0, max(Natm_n_imD$counts_Ser1)+2),
     col = col_imD[1], pch = 20)

lines(Natm_n_imD$allDate, Natm_n_imD$counts_meningitis,
      type = "b", col = col_imD[2], pch = 20)
lines(Natm_n_imD$allDate, Natm_n_imD$counts_30DDeath,
      type = "b", col = col_imD[3], pch = 20)
legend("topleft", names(col_imD), fill = col_imD, bty = "n")
dev.off()

# Viz per-week counts by base R plot
# https://stackoverflow.com/questions/30431444/plotting-by-week-with-ggplot-in-r
# https://stackoverflow.com/questions/3777174/plotting-two-variables-as-lines-using-ggplot2-on-the-same-graph
Natm_n_imD$weeks <- cut(Natm_n_imD[,"allDate"], breaks="week")
Nat_weekly <- Natm_n_imD %>% 
  dplyr::group_by(weeks) %>% 
  dplyr::summarise(counts_Ser1_weekly = sum(counts_Ser1),
                   counts_meningitis_weekly = sum(counts_meningitis),
                   counts_30DDeath_weekly = sum(counts_30DDeath)) %>% 
  dplyr::ungroup() #%>% 

dir.create("inputs")

incidence_weekly <- Nat_weekly %>% 
  dplyr::select(weeks, counts_Ser1_weekly) %>% 
  dplyr::rename(cases = counts_Ser1_weekly) # That annoying name

write.csv(incidence_weekly, "inputs/incidence_weekly.csv", row.names = FALSE)

png("pictures/weekly_cases.png", width = 17, height = 12, unit = "cm", res = 1200)
col_imD_weekly <- c(counts_Ser1_weekly = "deepskyblue3",
             counts_meningitis_weekly = "green",
             counts_30DDeath_weekly = "maroon")
ggplot(Nat_weekly, aes(as.Date(weeks))) +
  geom_line(aes(y = counts_Ser1_weekly, colour = "counts_Ser1_weekly")) +
  geom_line(aes(y = counts_meningitis_weekly, colour = "counts_meningitis_weekly")) +
  geom_line(aes(y = counts_30DDeath_weekly, colour = "counts_30DDeath_weekly")) +
  scale_x_date() +
  scale_color_manual(values = col_imD_weekly,
                     name = "Cases",
                     breaks = c("counts_Ser1_weekly", "counts_meningitis_weekly", "counts_30DDeath_weekly"),
                     labels = c("Serotype 1", "Meningitis", "30 Days Death")
  ) +
  ggtitle("The Counts of Serotype 1 in England") +
  xlab("Year") +
  ylab("Serotype 1 Cases (Aggregated by Week)")
dev.off()

## 2. Data Fitting #############################################################
# The anatomy of an mcstate particle filter, as noted above, consists of three main components: \n 
# 1. A set of observations to fit the model to, generated using mcstate::particle_filter_data(). \n 
# 2. A model to fit, which must be a dust generator, either dust::dust() or odin.dust::odin_dust(). \n 
# 3. A comparison function, which is an R function which calculates the likelihood of the state given the data at one time point.

# There is a calibration function in mcstate to fit our model to data.
# https://mrc-ide.github.io/mcstate/articles/sir_models.html
incidence <- Natm_n_imD %>% 
  dplyr::select(day, counts_Ser1) %>% 
  dplyr::rename(cases = counts_Ser1) # That annoying name

dir.create("inputs")
write.csv(incidence, "inputs/incidence.csv", row.names = FALSE)

png("pictures/hist_daily_cases.png", width = 17, height = 12, unit = "cm", res = 1200)
hist(incidence$cases,
     main = "Histogram of Daily Cases",
     xlab = "Daily Incidence") # huge zero daily cases occur
dev.off()

## 3. AMR Analysis #############################################################
# Load dat_G first.
link_ID <- readxl::read_excel("raw_data/gubbins/ukhsa_assemblies_02_07_24.xlsx")
link_ID$assembly_name <- substr(link_ID$assembly_name, 1, (nchar(link_ID$assembly_name)-6)) # That annoying last 6 chara of ".fasta"

amr_smx <- read.csv("raw_data/gubbins/resistance_folp_smx.csv")
amr_tmp <- read.csv("raw_data/gubbins/resistance_dhfr_tmp.csv")

dat_G_amr <- dplyr::left_join(dat_G, link_ID, by = "ID")
dat_G_amr <- dplyr::left_join(dat_G_amr, amr_smx, by = c("assembly_name" = "isolate_id"))
dat_G_amr <- dplyr::left_join(dat_G_amr, amr_tmp, by = c("assembly_name" = "isolate_id"))
dat_G_amr <- dat_G_amr %>% 
  dplyr::rename(resistance_smx = Resistance.x,
                resistance_tmp = Resistance.y)

