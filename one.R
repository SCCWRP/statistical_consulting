
calculate_AWM_bight18 <- function(B18_station, dd_inc = F, pah_inc = T, 
                                  PCB_change = F, excluded_strata = c()){
  B18_station <- B18_station %>% 
    filter(!stratum %in% excluded_strata) %>%
    mutate(
      result = case_when(
        result == -88 ~ 0,
        TRUE ~ result
      )
    )
  
  
  # Fines -------------------------------------------------------------------
  
  fines <- B18_station %>% 
    filter(analytename %in% c('Phi 4.5', 'Phi 5.0', 'Phi 5.5', 'Phi 6.0', 'Phi 6.5'
                              , 'Phi 7.0', 'Phi 7.5', 'Phi 8.0', 'Phi 8.5'
                              , 'Phi 9.0', 'Phi 9.5', 'Phi 10.0', 'Phi 10.5'
                              , 'Phi 11.0', 'Phi 11.5', 'Phi 12.0', 'Phi 12.5'
                              , 'Phi 13.0', 'Phi 13.5', 'Phi 14.0', 'Phi 14.5')) %>% 
    group_by(stratum, analyteclass, stationid) %>% 
    summarise(
      sum_result = sum(result, na.rm =T),
      denominator = max(chemareaweight),
      numerator = sum_result * denominator
    ) %>% 
    group_by(analyteclass) %>% 
    summarise(
      mean = sum(numerator, na.rm = T)/sum(denominator, na.rm = T), 
      sd = sqrt(sum(((sum_result - mean)*denominator)^2, na.rm = T)/
                  (sum(denominator, na.rm = T))^2),
      CI_95 = 1.96*sd, 
      mean_result = mean(sum_result, na.rm = T),
      sd_result = sd(sum_result, na.rm = T),
      CI_95_result = 1.96*sd_result,
      Min = min(na.omit(sum_result)),
      Max = max(na.omit(sum_result)),
      Median = median(na.omit(sum_result)),
      Pcnt_10 = quantile(na.omit(sum_result), .1),
      Pcnt_90 = quantile(na.omit(sum_result), .9),
      total = sum(sum_result)
    ) %>% 
    mutate(
      analyteclass = rep('Fines %')
    ) %>% 
    rename(analyte = analyteclass)
  
  
  
  # Inorganic ----------------------------------------------------------------
  
  inorganics <- B18_station %>% 
    filter(analyteclass == 'Inorganics') %>% 
    group_by(analytename) %>% 
    summarize(
      mean = sum(result *chemareaweight, na.rm =T)/
        sum(chemareaweight, na.rm = T),
      sd = sqrt(sum(((result - mean)*chemareaweight)^2, na.rm = T)/
                  (sum(chemareaweight, na.rm = T))^2),
      CI_95 = 1.96*sd,
      mean_result = mean(result, na.rm = T),
      sd_result = sd(result, na.rm = T),
      CI_95_result = 1.96*sd_result,
      Min = min(na.omit(result)),
      Max = max(na.omit(result)),
      Median = median(na.omit(result)),
      Pcnt_10 = quantile(na.omit(result), .1),
      Pcnt_90 = quantile(na.omit(result), .9)
    ) %>% 
    rename(analyte = analytename)
  
  
  AWM <- fines %>% 
    bind_rows(inorganics) %>% 
    arrange(analyte)
  return(AWM)
}



area_weights = ETL.Extract$station_completion.areas[c("Region...8","Total Area")][1:12,] %>% 
  filter(`Region...8` != "Bight") %>%
  rename(stratum = `Region...8`,
         chemareaweight = `Total Area`)
  
B18_station = merge(ETL.Extract$station_completion.trawls,ETL.Extract$chemistry_bight.tbl_chemresults,by="stationid") %>%
  merge(area_weights,by="stratum")
B18_station$result = as.numeric(B18_station$result)







x=calculate_AWM_bight18(B18_station)
