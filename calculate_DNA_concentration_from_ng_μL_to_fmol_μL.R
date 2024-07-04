calculate_DNA_concentration_from_ng_μL_to_fmol_μL <- function(DNA_Conc.ng_μL, DNA_length.bp) {
  
  library(tidyverse)
  
  Aver_molecular_weight_of_per_bp <- 617.96 # unit: g/mol
  molecular_weight_of_DNA <- DNA_length.bp * Aver_molecular_weight_of_per_bp # unit: g/mol
  DNA_Conc.fmol_μL <- round((DNA_Conc.ng_μL / molecular_weight_of_DNA) * 10 ^ 6, 2)
  
  
  
  
  output <- paste("DNA concentration (fmol/μL):", DNA_Conc.fmol_μL)
  
  return(output)
}
