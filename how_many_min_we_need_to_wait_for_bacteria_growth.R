how_many_min_we_need_to_wait_for_bacteria_growth <- function(current_OD, doubling_time_in_min) {
  
  library(tidyverse)
  
  
  
  OD600 <- 5 * (10 ^ 8) # how many cells per mL, when OD600 = 1
  
  current_cell_number <- current_OD * OD600 # how many cells we have right now
  
  
  
  # how many cells we want
  aimmed_OD <- c(0.4, 0.5)
  target_cell_number_lowest <- aimmed_OD[1] * OD600
  target_cell_number_highest <- aimmed_OD[2] * OD600
  
  time_to_double_shortest <- round(log2(target_cell_number_lowest / current_cell_number) * doubling_time_in_min, 0)
  time_to_double_longest <- round(log2(target_cell_number_highest / current_cell_number) * doubling_time_in_min, 0)
  
  
  
  print(paste("time need to wait:", time_to_double_shortest, "~", time_to_double_longest, "min"))
}
