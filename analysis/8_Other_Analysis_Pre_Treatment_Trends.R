###############################################################################
# Author: FYO
# CGIC evaluations: processing outcome data for CGIC evaluation
# 
# Last Updated: 
###############################################################################

# Packages needed for analysis
library(tidyverse)
library(vtable)
library(lubridate)
library(digest)
library(purrr)
library(readr)
library(reshape)
library(R.utils)

################################ Loading data #################################
# To load, open this file AFTER loading the .Rproj file, which sets the 
# working directory to /CGIC_Evaluation_2019/analysis/ 
# The files here come from the R script titled 6_Quasi_Experimental_Analysis_Diff_in_Diff.R

daily_shotspotter <- read.delim(file = "../data/daily_shotspotter.csv", sep = ",", stringsAsFactors = FALSE) %>%
  mutate(event_date = ymd(event_date))
daily_calls <- read.delim(file = "../data/daily_calls.csv", sep = ",", stringsAsFactors = FALSE)%>%
  mutate(event_date = ymd(event_date))
daily_dcr_violent <- read.delim(file = "../data/daily_dcr_violent.csv", sep = ",", stringsAsFactors = FALSE)%>%
  mutate(event_date = ymd(event_date))
daily_dcr_gun <- read.delim(file = "../data/daily_dcr_gun.csv", sep = ",", stringsAsFactors = FALSE)%>%
  mutate(event_date = ymd(event_date))
daily_arrest_broad <- read.delim(file = "../data/daily_arrest_broad.csv", sep = ",", stringsAsFactors = FALSE)%>%
  mutate(event_date = ymd(event_date))
daily_arrest_violation <- read.delim(file = "../data/daily_arrest_violation.csv", sep = ",", stringsAsFactors = FALSE)%>%
  mutate(event_date = ymd(event_date))
daily_arrest_violent <-read.delim(file = "../data/daily_arrest_violent.csv", sep = ",", stringsAsFactors = FALSE)%>%
  mutate(event_date = ymd(event_date))


##################### Creating Pre-treatment Trend Figures ####################
# Note that these figures do not appear in the report, but were used to assess
# The potential vailidity of the parallel trends assumption prior to running DiD.

# Shotspotter
monthly_agg_shotspotter <- daily_shotspotter %>%
  filter(event_date >= ymd("2016-11-01") & event_date < ymd("2017-11-01")) %>% 
  group_by(Treated, month_year) %>%
  summarise(event_count = sum(shotspotter_count))

shotspotter_trend_plot <- ggplot(data = monthly_agg_shotspotter, 
                                 aes(x = month_year, y = event_count, group = Treated, color = Treated)) +
  geom_line() + theme_bw() + ylab("Event Count") + xlab("Date") + geom_vline(xintercept = ymd("2017-11-01"), colour = "grey60", linetype = 2) +
  ggtitle("Outcome: Shotspotter")

shotspotter_trend_plot
ggsave("../output/shotspotter_trend.pdf")

# Calls for Service 
monthly_agg_calls <- daily_calls %>%
  filter(event_date >= ymd("2016-11-01") & event_date < ymd("2017-11-01")) %>% 
  group_by(Treated, month_year) %>%
  summarise(event_count = sum(calls_count))

calls_trend_plot <- ggplot(data = monthly_agg_calls, aes(x = month_year, y = event_count, group = Treated, color = Treated)) +
  geom_line() + theme_bw() + ylab("Event Count") + xlab("Date") + geom_vline(xintercept = ymd("2017-11-01"), colour = "grey60", linetype = 2) +
  ggtitle("Outcome: Calls for Service for Sounds of Gunshots")

calls_trend_plot
ggsave("../output/calls_trend.pdf")


# DCR Violent Crime
monthly_agg_dcr_violent <- daily_dcr_violent %>%
  filter(event_date >= ymd("2016-11-01") & event_date < ymd("2017-11-01")) %>% 
  group_by(Treated, month_year) %>%
  summarise(event_count = sum(dcr_violent_count))

dcr_violent_trend_plot <- ggplot(data = monthly_agg_dcr_violent, aes(x = month_year, y = event_count, group = Treated, color = Treated)) +
  geom_line() + theme_bw() + ylab("Event Count") + xlab("Date") + geom_vline(xintercept = ymd("2017-11-01"), colour = "grey60", linetype = 2) +
  ggtitle("Outcome: Violent Crime")

dcr_violent_trend_plot
ggsave("../output/dcr_violent_trend.pdf")


# DCR Gun Crime
monthly_agg_dcr_gun <- daily_dcr_gun %>%
  filter(event_date >= ymd("2016-11-01") & event_date < ymd("2017-11-01")) %>% 
  group_by(Treated, month_year) %>%
  summarise(event_count = sum(dcr_gun_count))

dcr_gun_trend_plot <- ggplot(data = monthly_agg_dcr_gun, aes(x = month_year, y = event_count, group = Treated, color = Treated)) +
  geom_line() + theme_bw() + ylab("Event Count") + xlab("Date") + geom_vline(xintercept = ymd("2017-11-01"), colour = "grey60", linetype = 2) +
  ggtitle("Outcome: Violent Gun Crime")

dcr_gun_trend_plot
ggsave("../output/dcr_gun_trend.pdf")

# Broad arrest
monthly_agg_arrest_broad <- daily_arrest_broad %>%
  filter(event_date >= ymd("2016-11-01") & event_date < ymd("2017-11-01")) %>% 
  group_by(Treated, month_year) %>%
  summarise(event_count = sum(arrest_broad_count))

arrest_broad_trend_plot <- ggplot(data = monthly_agg_arrest_broad, aes(x = month_year, y = event_count, group = Treated, color = Treated)) +
  geom_line() + theme_bw() + ylab("Event Count") + xlab("Date") + geom_vline(xintercept = ymd("2017-11-01"), colour = "grey60", linetype = 2) +
  ggtitle("Outcome: Gun Arrests")

arrest_broad_trend_plot
ggsave("../output/arrest_broad_trend.pdf")


# Gun Violence Arrests
monthly_agg_arrest_violent <- daily_arrest_violent %>%
  filter(event_date >= ymd("2016-11-01") & event_date < ymd("2017-11-01")) %>% 
  group_by(Treated, month_year) %>%
  summarise(event_count = sum(arrest_violent_count))

arrest_violent_trend_plot <- ggplot(data = monthly_agg_arrest_violent, aes(x = month_year, y = event_count, group = Treated, color = Treated)) +
  geom_line() + theme_bw() + ylab("Event Count") + xlab("Date") + geom_vline(xintercept = ymd("2017-11-01"), colour = "grey60", linetype = 2) +
  ggtitle("Outcome: Arrests for Violent Crimes Involving Guns")

arrest_violent_trend_plot
ggsave("../output/arrest_violent_trend.pdf")

# Gun Violation Arrests
monthly_agg_arrest_violation <- daily_arrest_violation %>%
  filter(event_date >= ymd("2016-11-01") & event_date < ymd("2017-11-01")) %>% 
  group_by(Treated, month_year) %>%
  summarise(event_count = sum(arrest_weapon_violation_count))

arrest_violation_trend_plot <- ggplot(data = monthly_agg_arrest_violation, aes(x = month_year, y = event_count, group = Treated, color = Treated)) +
  geom_line() + theme_bw() + ylab("Event Count") + xlab("Date") + geom_vline(xintercept = ymd("2017-11-01"), colour = "grey60", linetype = 2) +
  ggtitle("Outcome: Gun Violation Arrests")

arrest_violation_trend_plot
ggsave("../output/arrest_violation_trend.pdf")

