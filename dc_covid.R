# *** DC COVID-19 Vulnerability & Epidemiological Data

# ----------
# Setup
# ----------

library(tidyverse)
library(data.table)
library(readxl)
library(lubridate)
library(viridis)
install.packages("boxr")
library(boxr)


# ----------
# 1: Data Importing and Cleaning
# ----------

make_path <- function(file_name) {
  return(file.path(getwd(), "data", file_name))}

pop_tract_raw <- read.csv(make_path("DC_Population_by_Tract.csv"))
neighborhood_to_tract_raw <- read.csv(make_path("DC_Health_Planning_Neighborhoods_to_Census_Tracts.csv"))

# From: https://dcgov.app.box.com/v/DCHealthStatisticsData
pos_neighborhood_raw <- read_xlsx(make_path("COVID19_DCHealthStatisticsDataV3 (NewFileStructure).xlsx"), 
                                  "Total Positives by Neighborhood")
test_neighborhood_raw <- read_xlsx(make_path("COVID19_DCHealthStatisticsDataV3 (NewFileStructure).xlsx"), 
                                   "Total Tests by Neighborhood")



# ----------
# 2. Data Processing: Join Data, Aggregate CCVI, & Calculate Epi Metrics
# ----------

# Step 1: Join Population & Aggregate to Neighborhood Level

# Join tract population and Neighborhood-to-Tract
neighborhood_to_tract_joined <- full_join(neighborhood_to_tract_raw, pop_tract_raw, 
                               by = "GEOID")

#  Clean Columns
neighborhood_to_tract_joined <- neighborhood_to_tract_joined %>% select(-c("OBJECTID.x", "OBJECTID.y",
                                                     "STATE", "COUNTY"))
neighborhood_to_tract_joined <- neighborhood_to_tract_joined %>% rename_all(tolower)

neighborhood_pop <- neighborhood_to_tract_joined %>% 
  group_by(dc_hpn_name, code, hpn_label) %>% 
  summarize(neighborhood_pop = sum(population))

neighborhood_pop <- neighborhood_pop %>% arrange(hpn_label)



# Step 3: Aggregate Epi (test & case) data into weekly  tables
# Due to variations in the hours of COVID testing sites and labs, COVID tests
# and cases can vary dramatically based on the day of the week. Therefore, 7 day 
# running averages or weekly averages are recommended. 

process_epi_data <- function(df) {
  # TODO Future: Add option to summarize by month
  df <- rename(df, total = 3)
  
  df <- df %>% 
    mutate(epiweek = epiweek(Date)) %>% 
    arrange(Neighborhood, epiweek) %>%
    group_by(Neighborhood) %>%
    mutate(new = total - lag(total))
  
  df <- df%>%
    group_by(Neighborhood, epiweek) %>%
    summarise(days_w_data = n(), 
              new = sum(new),
              total = max(total))
  
  # Excluding weeks with less than 7 days represented
  # TODO Future: Explore in detail the dates with missing data and decide whether to include
  df <- df %>% filter(days_w_data %in% c(5, 6, 7))
  df <- df %>% select(-days_w_data)

  return(df)
}

tests_weekly <- process_epi_data(test_neighborhood_raw)
tests_weekly$new[tests_weekly$new < 0] <- NA
tests_weekly <- tests_weekly %>% rename(new_tests = new, total_tests = total)

cases_weekly <- process_epi_data(pos_neighborhood_raw)
cases_weekly$new[cases_weekly$new < 0] <- NA
cases_weekly <- cases_weekly %>% rename(new_cases = new, total_cases = total)

epi_weekly <- full_join(tests_weekly, cases_weekly, by = c("Neighborhood", "epiweek"))

epi_weekly$code <- sub(": .*","", epi_weekly$Neighborhood)

temp <- data.frame(neighborhood_pop %>% select(code, neighborhood_pop))
epi_weekly <- left_join(epi_weekly, temp,  by = "code")

epi_weekly <- epi_weekly %>% filter(Neighborhood != "Unknown")

epi_weekly <- epi_weekly %>% mutate(
  tpr = (new_cases / new_tests *100),
  test_rate = (new_tests/neighborhood_pop * 1000),
  case_rate = (new_cases/neighborhood_pop * 100000)
)

epi_weekly <- left_join(epi_weekly, epi_weekly %>% 
                          group_by(Neighborhood) %>%
                          summarise(cumulative_incidence = (max(total_cases) / max(neighborhood_pop) * 100000)))

# Pivot Long
epi_weekly <- epi_weekly %>% select(-c("code", "hpn_label"))

epi_weekly_long <- epi_weekly %>% 
  pivot_longer(cols = -c("Neighborhood", "epiweek", "ccvi_sum", "ccvi_quintile"),
               names_to = "variable", values_to = "value")
epi_weekly_long$ccvi_quintile <- as.factor(epi_weekly_long$ccvi_quintile)


# ----------
# 3. Analysis & Visualization
# ----------

## VISUALIZE 
# I chose a heat map because it is a  simple way to look at two variables 
# over a time series. 
# X axis = epidemiological weeks
# Y axis = neighborhoods, in order of vulnerability (most vulnerable first)
# Color = epi variables. I chose TPR and Testing Rate, because the lower 
  # testing rates indicated by these variables will suppress the differences 
  # in incidence rate that is likely happening in the ground truth

# TODO Future: create choropleth maps using the shapefiles on opendata.dc 

# TPR Heat Map
epi_weekly %>% 
  filter(epiweek != 35 & epiweek != 36 & 
           epiweek >= 22) %>% 
  arrange(cumulative_incidence) %>%
  ggplot(aes(x = epiweek, 
             y = reorder(Neighborhood, cumulative_incidence), 
             fill = tpr)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_viridis(discrete=FALSE, option="A") +
  labs(title = "DC COVID: Neighborhood Test Positivity Rate, May 24 to November 27") +
  ylab("Health Neighborhood, Ordered By Cumulative Incidence") +
  xlab("Epidemiological Weeks") +
  labs(title = "Test Positivity") +
  scale_x_continuous(breaks = c(25, 30, 35, 40, 45, 50)) + 
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank()) +
  ggsave("tpr_viz.png", path = "output", width = 6, height = 8)

# Test Rate Heat Map
epi_weekly %>% 
  filter(epiweek >= 22) %>% 
  arrange(cumulative_incidence) %>%
  ggplot(aes(x = epiweek, 
             y = reorder(Neighborhood, cumulative_incidence), 
             fill = test_rate)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_viridis(discrete=FALSE, option="A") +
  labs(title = "DC COVID: Neighborhood Testing Rate, May 24 to November 27") +
  ylab("Health Neighborhood, Ordered By Cumulative Incidence") +
  xlab("Epidemiological Weeks") +
  labs(title = "Test Rate / 1000 Residents") +
  scale_x_continuous(breaks = c(25, 30, 35, 40, 45, 50)) + 
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank()) +
  ggsave("testing_rate_viz.png", path = "output", width = 6, height = 8)

# Incidence Rate Map
epi_weekly %>% 
  filter(epiweek != 35 & epiweek != 36 & 
           epiweek >= 22) %>% 
  arrange(cumulative_incidence) %>%
  ggplot(aes(x = epiweek, 
             y = reorder(Neighborhood, cumulative_incidence), 
             fill = case_rate)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_viridis(discrete=FALSE, option="A") +
  labs(title = "DC COVID: Neighborhood Incidence, May 24 to November 27") +
  ylab("Health Neighborhood, Ordered By Cumulative Incidence") +
  xlab("Epidemiological Weeks") +
  labs(title = "Incidence Rate / 100,000") +
  scale_x_continuous(breaks = c(25, 30, 35, 40, 45, 50)) + 
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank()) + 
  ggsave("incidence_rate_viz.png", path = "output", width = 6, height = 8)

# Line Plot - serious spaghetti lines, but wanted to check it out just in case
# a clear pattern emerged
epi_weekly_long %>% 
  filter(variable %in% c("test_rate") & epiweek != 35 & epiweek != 36) %>%
  ggplot(aes(x = epiweek, y = value, group=Neighborhood, color = ccvi_quintile)) +
  geom_line() +
  theme_minimal() +
  theme(legend.position = "bottom")

# QUINTILES
# Less spaghetti, but also less nuanced and more difficult to interpret. 
quintile_summary <- epi_weekly %>% group_by(ccvi_quintile, epiweek) %>%
  mutate(tpr_weighted = (tpr*pop_weight),
         test_rate_weighted = (test_rate*pop_weight), 
         case_rate_weighted = (case_rate*pop_weight)) %>%
  summarise(tpr = sum(tpr_weighted), 
            test_rate = sum(test_rate_weighted), 
            case_rate = sum(case_rate_weighted))
quintile_summary <- quintile_summary %>% 
  pivot_longer(cols = -c("epiweek", "ccvi_quintile"),
               names_to = "variable", values_to = "value")
quintile_summary$ccvi_quintile <- as.factor(quintile_summary$ccvi_quintile)

quintile_summary %>% 
  filter(epiweek != 35 & epiweek != 36 & variable %in% ("case_rate")) %>%
  ggplot(aes(x = epiweek, y = value, group=ccvi_quintile, color = ccvi_quintile)) +
  geom_line() +
  scale_fill_viridis() +
  labs(title = "ADD TITLE") +
  theme_minimal() +
  theme(legend.position = "bottom") + 
  facet_grid(rows = vars(variable))

# REGRESSION
# Without identifying how testing policy and access have changed over the time
# period of the dataset, there are too many potential confounders for regression
# and machine learning techniques to be used. However, I was curious so I did a 
# simple linear regression with the three outcome variables. 
epi_weekly <- setDT(epi_weekly)

lr_tpr <- lm(case_rate ~ epiweek + theme1_sum + theme2_sum + theme3_sum + theme4_sum 
              + theme5_sum + theme6_sum + ccvi_sum,
              data = epi_weekly[is.na(epi_weekly$case_rate) == FALSE])
summary(lr_tpr)




  


