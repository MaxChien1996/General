volume_loading_for_Nanopore <- function(DNA_Conc.ng_μL, Aver_Weight_of_library.bp){

  library(tidyverse)

  Aver_MW_of_per_bp <- 617.96 #unit: g/mol
  Aver_MW_of_library <- Aver_Weight_of_library.bp * Aver_MW_of_per_bp #unit: g/mol
  DNA_Conc.fmol_μL <- (DNA_Conc.ng_μL / Aver_MW_of_library) * 10 ^ 6

  n <- seq(0, 5, by = 0.01) #0 and 5 represent minimum and maximum volume can load
  all_df <- data.frame(volume = 0)
  for(volume in n){
    if(volume * DNA_Conc.fmol_μL < 50 & volume * DNA_Conc.fmol_μL > 20)
    {
      df <- data.frame(volume)
      all_df <- bind_rows(all_df, df)
    } else {
      volume = volume
    }
  }
  all_df <- all_df %>%
    filter(volume != 0)

  print("DNA loading should be between 20 ~ 50 fmol")
  print(paste("DNA_Conc. (fmol/μL):", DNA_Conc.fmol_μL))

  output <- paste("loading volume:", min(all_df), "~", max(all_df), "μL", "(maximum volume shouldn't be more than 5 μL)")

  return(output)
}
