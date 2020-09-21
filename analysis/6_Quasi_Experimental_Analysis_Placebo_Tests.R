###############################################################################
# CGIC evaluations: Implementing Placebo tests
# Author: FYO
###############################################################################

# Packages needed for analysis
library(tidyverse)
library(vtable)
library(lubridate)
library(digest)
library(purrr)
library(readxl)
library(reshape)
library(gridExtra)

# Saving lab colors for plotting
lab_colors <- c("#2b4888", "#de4057")

################################ Loading data #################################

# Importing rolling means
rolling_means <- read.delim(file = "../data/rolling_means.csv", 
                            sep = ",", stringsAsFactors = FALSE)

# Creating list of dates
start_date = ymd('2016-11-01')
end_date = ymd('2019-04-30')
date_list <- seq(start_date, end_date, by = "days")


############################ Processing Function ##############################

clean_placebo_district <- function(data_path = "../data/data_placebos/",
                                   real_pre_name = "1D_calls_real_pre.csv",
                                   real_post_name = "1D_calls_real_post.csv",
                                   synth_pre_name = "1D_calls_synth_pre.csv",
                                   synth_post_name = "1D_calls_synth_post.csv") {
  
  # Description: Combines and cleanssynth results generated using python. Note that this is
  #              called by the "create_synth_plots" and "create_placebo_plots" functions.
  #
  # Args: data_path = file path for input data
  #       real_pre_name = file name for .csv containing real pretreatment outcome measure
  #       real_post_name = file name for .csv containing real posttreatment outcome measure
  #       synth_pre_name = file name for .csv containing synthetic pretreatment outcome measure
  #       synth_post_name = file name for .csv containing synthetic posttreatment outcome measure
  #
  # Returns: real and synthetic control values for a given outcome and district
  
  # Creating filename strings for import
  real_pre_path <- paste(data_path, real_pre_name, sep = "")
  real_post_path <- paste(data_path, real_post_name, sep = "")
  synth_pre_path <- paste(data_path, synth_pre_name, sep = "")
  synth_post_path <- paste(data_path, synth_post_name, sep = "")
  
  # reading in data
  real_pre <- read.delim(file = real_pre_path, sep = ",", stringsAsFactors = FALSE)
  real_post <- read.delim(file = real_post_path, sep = ",", stringsAsFactors = FALSE)
  
  synth_pre <- read.delim(file = synth_pre_path, sep = ",", stringsAsFactors = FALSE)
  synth_post <- read.delim(file = synth_post_path, sep = ",", stringsAsFactors = FALSE)
  
  # Appending pre and post
  real_combined <- bind_rows(real_pre, real_post)
  synth_combined <- bind_rows(synth_pre, synth_post)
  
  
  # Converting from wide (PSAs in different columns) to long, adding date, merging synth and real
  real_long <- real_combined %>%
    gather(event_psa, real_value, names(real_combined)) %>%
    mutate(event_psa = substring(event_psa,2))
  
  real_long$event_date <- rep(date_list, length(real_combined))
  
  synth_long <- synth_combined %>%
    gather(event_psa, synth_value, names(synth_combined)) %>%
    mutate(event_psa = substring(event_psa,2))
  
  synth_long$event_date <- rep(date_list, length(synth_combined))
  
  combined <- real_long %>%
    inner_join(synth_long, by = c("event_psa", "event_date")) %>%
    mutate(prediction_error = real_value - synth_value) %>%
    mutate(squared_prediction_error = prediction_error ^ 2) %>%
    mutate(treatment_period = event_date >= ymd('2017-11-01'))
  
  return(combined)
  
}

############################### Placebo plots  ################################

create_placebo_plots <- function(data_path = "../data/data_placebos/",
                                 outcome = "calls",
                                 district_list = c("1D", "2D", "3D", "4D", "5D", "6D", "7D"),
                                 figure_title = "Average Monthly gap in Calls for Service for Sounds of Gunshots for treated and placebo PSAs",
                                 output_name = "calls_placebo_figure.pdf",
                                 output_path = "../output/") {
  
  
  # Description: Creates placebo test line plots
  #
  # Args: data_path = file path for input data
  #       outcome = name of the outcome variable to be used
  #       district_list = name of the districts to plot results for (defaults to entire district list, but can change to do subsets)
  #       Figure title = string, desired title for  figure.
  #       output_name = name for output file (string ending with desired file extension)
  #       output_path = directory where you want to save results
  #
  # Returns: cleaned and combined synth results for all districts listed for the chosen outcome. Note that
  #          The output from this function is what feeds into the calculate_mspe function.
  
  
  for (district in district_list) {
    
    # Naming files
    real_pre = paste(district, "_",outcome, "_real_pre.csv", sep = "")
    real_post = paste(district, "_", outcome, "_real_post.csv", sep = "")
    synth_pre = paste(district, "_", outcome, "_synth_pre.csv", sep = "")
    synth_post = paste(district, "_", outcome, "_synth_post.csv", sep = "")
    
    # Combining and cleaning up placebo info for each district with fxn above
    assign(paste("temp_", district, sep = ""), clean_placebo_district(data_path = data_path,
                                                                      real_pre_name = real_pre,
                                                                      real_post_name = real_post,
                                                                      synth_pre_name = synth_pre,
                                                                      synth_post_name = synth_post))
    # Combining cleaned district files
    if (district == district_list[1]) {
      placebo_combined <- eval(parse(text = paste("temp_", district, sep = "")))
      
    } else {
      placebo_combined <- rbind(placebo_combined, eval(parse(text = paste("temp_", district, sep = ""))))
    }
  }
  
  # Marking treated PSAs
  placebo_combined <- placebo_combined %>%
    mutate(`Treated PSA` = ifelse(event_psa %in% c(701:708), "Treated PSA", "Untreated PSA"))
  
  # Creating monthly averages for plotting
  placebo_monthly <- placebo_combined %>%
    mutate(year = year(event_date)) %>%
    mutate(month = month(event_date)) %>%
    mutate(month_year =  ymd(paste(year, month, "01", sep = "-"))) %>%
    group_by(month_year, event_psa, `Treated PSA`) %>%
    summarise(avg_error = mean(prediction_error))
  
  # Plotting outcome over time
  placebo_plot <- ggplot(data = placebo_monthly, aes(x = month_year, y = avg_error, group = event_psa, colour = `Treated PSA`)) +
    geom_line() + theme_bw() + scale_color_manual(values = lab_colors) +
    ylab("avg. monthly gap b/w real and synthetic values") + xlab("") + 
    geom_vline(xintercept = ymd("2017-11-01"), colour = "grey60", linetype = 2) +
    ggtitle(figure_title) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.background = element_rect(colour = "grey80")
    )
  
  # View and save plot
  placebo_plot
  
  save_string <- paste(output_path, output_name, sep = "")
  ggsave(save_string)
  
  # Returning combined and processed placebo data
  return(placebo_combined)
}


calls_placebo_data <- create_placebo_plots(data_path = "../data/data_placebos/",
                                           outcome = "calls",
                                           district_list = c("1D", "2D", "3D", "4D", "5D", "6D", "7D"),
                                           figure_title = "Average Monthly Gap in Calls for Service for Sounds of Gunshots",
                                           output_name = "calls_placebo_figure.pdf")

shotspotter_placebo_data <- create_placebo_plots(data_path = "../data/data_placebos/",
                                                 outcome = "shotspotter",
                                                 district_list = c("1D", "3D", "4D", "5D", "6D", "7D"),
                                                 figure_title = "Average Monthly Gap in Shotspotter Alerts",
                                                 output_name = "shotspotter_placebo_figure.pdf")

dcr_violent_placebo_data <- create_placebo_plots(data_path = "../data/data_placebos/",
                                                 outcome = "dcr_violent",
                                                 district_list = c("1D", "2D", "3D", "4D", "5D", "6D", "7D"),
                                                 figure_title = "Average Monthly Gap in Violent Crime",
                                                 output_name = "dcr_violent_placebo_figure.pdf")

dcr_gun_placebo_data <- create_placebo_plots(data_path = "../data/data_placebos/",
                                             outcome = "dcr_gun",
                                             district_list = c("1D", "2D", "3D", "4D", "5D", "6D", "7D"),
                                             figure_title = "Average Monthly Gap in Gun Crime",
                                             output_name = "dcr_gun_placebo_figure.pdf")

arrest_broad_placebo_data <- create_placebo_plots(data_path = "../data/data_placebos/",
                                                  outcome = "arrest_broad",
                                                  district_list = c("1D", "2D", "3D", "4D", "5D", "6D", "7D"),
                                                  figure_title = "Average Monthly Gap in Gun Arrests",
                                                  output_name = "arrest_broad_placebo_figure.pdf")

arrest_violation_placebo_data <- create_placebo_plots(data_path = "../data/data_placebos/",
                                                      outcome = "arrest_violation",
                                                      district_list = c("1D", "2D", "3D", "4D", "5D", "6D", "7D"),
                                                      figure_title = "Average Monthly Gap in Arrests for Gun Violations",
                                                      output_name = "arrest_violation_placebo_figure.pdf")

arrest_violent_placebo_data <- create_placebo_plots(data_path = "../data/data_placebos/",
                                                    outcome = "arrest_violent",
                                                    district_list = c("1D", "2D", "3D", "4D", "5D", "6D", "7D"),
                                                    figure_title = "Average Monthly Gap in Arrests for Violent Gun Crimes",
                                                    output_name = "arrest_violent_placebo_figure.pdf")

################################# MSPE Plots ##################################

calculate_mspe <- function(data = calls_placebo_data,
                           figure_title = "Mean Square Prediction Error Ratio: Calls for Service",
                           save_name = "calls_mspe",
                           output_path_figure ="../output/",
                           output_path_data = "../output/MPSE_Calculations/"
) {
  
  # Description: Creates placebo test line plots
  #
  # Args: data = name of df containing intput data
  #       figure_title = string, desired title for  figure.
  #       save_name = string, file name for output (no file extension needed, automatically saves .pdf)
  #       output_path_figure = string, file path for folder where figures should be saved
  #       output_path_data = string, file path for folder where calculations should be saved
  #
  # Returns: currently returns the plot (but can return data by commenting out the current return call and uncommenting the other)
  
  
  # Calculating MSPE before and after treatment for each
  mspe <- data %>%
    group_by(event_psa, treatment_period, `Treated PSA`) %>%
    summarise(mspe = mean(squared_prediction_error)) %>%
    spread(treatment_period, mspe) %>%
    select(event_psa, `Treated PSA`, posttreatment_mspe = `TRUE`, pretreatment_mspe = `FALSE`) %>%
    mutate(`Post MSPE / Pre MSPE` = posttreatment_mspe / pretreatment_mspe)
  
  # Making Treated PSA an ordered factor so that treated dots are on top
  mspe$`Treated PSA` <- factor(mspe$`Treated PSA`, levels = c("Treated PSA", "Untreated PSA", ordered = TRUE))
  
  mspe_plot <- ggplot(data = mspe, aes(x = `Post MSPE / Pre MSPE`, colour = `Treated PSA`, fill = `Treated PSA`)) +
    geom_histogram(binwidth = .5, alpha = 0.5) + ylab("") + theme_bw() + 
    scale_color_manual(values = lab_colors) + scale_fill_manual(values = lab_colors) +
    ggtitle(figure_title) +
    theme(plot.title = element_text(face = "bold", size = 11),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.background = element_rect(colour = "grey80"))
  
  mspe_plot
  
  # Saving results
  save_path_pdf <- paste(output_path_figure, save_name, ".pdf", sep = "")
  save_path_csv <- paste(output_path_data,save_name, ".csv", sep = "")
  ggsave(save_path_pdf)
  write.csv(mspe, file = save_path_csv)
  
  # Retuning underlying data for plot
  #return(mspe)
  
  # Returning the plot
  return(mspe_plot)
}

mspe_calls <- calculate_mspe(data = calls_placebo_data, save_name = "calls_mspe", 
                             figure_title = "Mean Square Prediction Error Ratio: Calls for Service")
mspe_shotspotter <- calculate_mspe(data = shotspotter_placebo_data, save_name = "shotspotter_mspe", 
                                   figure_title = "Mean Square Prediction Error Ratio: Shotspotter Alerts")
mspe_dcr_violent <- calculate_mspe(data = dcr_violent_placebo_data, save_name = "dcr_violent_mspe", 
                                   figure_title = "Mean Square Prediction Error Ratio: Violent Crime")
mspe_dcr_gun <- calculate_mspe(data = dcr_gun_placebo_data, save_name = "dcr_gun_mspe", 
                               figure_title = "Mean Square Prediction Error Ratio: Gun Crime")
mspe_arrest_broad <- calculate_mspe(data = arrest_broad_placebo_data, save_name = "arrest_broad_mspe", 
                                    figure_title = "Mean Square Prediction Error Ratio: Gun Arrests")
mspe_arrest_violation <- calculate_mspe(data = arrest_violation_placebo_data, save_name = "arrest_violation_mspe", 
                                        figure_title = "Mean Square Prediction Error Ratio: Gun Violation Arrests")
mspe_arrest_violent <- calculate_mspe(data = arrest_violent_placebo_data, save_name = "arrest_violent_mspe",
                                      figure_title = "Mean Square Prediction Error Ratio: Arrests for Violent Gun Crime")
