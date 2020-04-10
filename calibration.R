## calibration.R: 


with(ETL.Load,{
  calibration_data = ETL.Transform$calibration_data
})


#############################################
#                                           #
# Make functions available                  #
#                                           #  
#############################################


### Methods to normalize residuals of calibration curves
with(calibration.env,{

  # creating a wrapper function in case we choose to augment with mahalanobis distance, 
  # or add unsupervised learning
    normalize_residuals = function(...,method="SCCWRP"){
      if (method=="SCCWRP"){
        function_call = normalize_residuals.sccwrp(...)
        return(function_call)
      }
    }
  
    
    
    
    
    
    
    
    
    
# SCCWRP's preferred methodology:
    
    normalize_residuals.sccwrp = function(clean_sites,rm_column,rm_value,tm_column,tm_value,maxiter=100){

      ####################################################################################
      # data expected/provided to function
      ####################################################################################
      
      #clean_data:
      # stationid      lat      long     Stratum TraceMetal  PPM ReferenceMetal   PPH
      # 1 B18-10203 34.44325 -120.4297 Inner Shelf    Arsenic 6.28           Iron 0.893
      # 2 B18-10206 34.39825 -119.8652   Mid Shelf    Arsenic 4.58           Iron 0.761
      # 3 B18-10208 34.33409 -119.4346 Inner Shelf    Arsenic 8.25           Iron 1.500
      # 4 B18-10209 34.28353 -119.3544 Inner Shelf    Arsenic 5.93           Iron 1.050
      # 5 B18-10210 34.24356 -119.3853 Inner Shelf    Arsenic 5.53           Iron 1.360
      # 6 B18-10211 34.22837 -119.3520 Inner Shelf    Arsenic 5.14           Iron 1.430
      
      #rm_column:
      # [1] "ReferenceMetal"
      
      #rm_value:
      # [1] "PPH"
      
      #tm_column:
      # [1] "TraceMetal"
      
      #tm_value:
      # [1] "PPM"
      
      ####################################################################################
      # function_code
      ####################################################################################
      
    trace_metals = unique(clean_sites[[tm_column]])
    reference_metals = unique(clean_sites[[rm_column]])
    colnames_clean_sites = colnames(clean_sites)
    return_clean_sites = clean_sites[0,]
    return_new_contaminated_sites = clean_sites[0,]
    
    for (trace_metal in trace_metals){
      for (reference_metal in reference_metals){
        #print(c(contaminant,trace_metal))
        trace_metal_rows = clean_sites[[tm_column]]==trace_metal
        reference_metal_rows = clean_sites[[rm_column]]==reference_metal
        selector = (trace_metal_rows + reference_metal_rows)==2
        
        clean_sites_tm_rm = clean_sites[selector,]
        contaminated_sites_tm_rm = clean_sites[0,]
        
        sw=0
        i=0
        while(sw<=0.05){
          formula = paste(tm_value,"~",rm_value)
          model = lm(formula ,data=clean_sites_tm_rm)
          
          tryCatch({sw = shapiro.test(resid(model))$p.val},error = function(e){i = maxiter+1})
          
          if(i>maxiter){
            #  print("Something has gone horribly wrong")
            clean_sites_tm_rm = clean_sites_tm_rm[0,]
            break
          } else if (i==0 && sw > 0.05){
            break
          }
          i=i+1
          clean_sites_tm_rm$residuals = NA
          clean_sites_tm_rm$outliers  = NA
          
          clean_sites_tm_rm$residuals = residuals(model) #Store Residuals
          SD2 = 2*sd(clean_sites_tm_rm$residuals)      #Compute 2 sd
          clean_sites_tm_rm$outliers = ifelse(abs(clean_sites_tm_rm$residuals)>SD2, 1, 0)
          contaminated_sites_tm_rm = rbind(contaminated_sites_tm_rm,clean_sites_tm_rm[clean_sites_tm_rm$outliers == 1,])
          clean_sites_tm_rm = clean_sites_tm_rm [clean_sites_tm_rm$outliers != 1,]
          
        }
        return_clean_sites = rbind(return_clean_sites,clean_sites_tm_rm[,colnames_clean_sites])
        return_new_contaminated_sites = rbind(return_new_contaminated_sites,contaminated_sites_tm_rm[,colnames_clean_sites])
        
      }}
    
    return(list(clean_sites = return_clean_sites,new_contaminated_sites = return_new_contaminated_sites))
    
  }
  
  
})







### methods to generate reference metal overview table

with(calibration.env,{
  
  
  rm_model_summary_kable = function(reference_metal,normalized_data,rm_column,rm_value,tm_column,tm_value,calibration_sites = F,slope_test = F){
    ####################################################################################
    # data expected/provided to function
    ####################################################################################
    
    # reference_metal:
    # [1] "Aluminum"
    
    #normalized_data:
    # stationid      lat      long     Stratum TraceMetal  PPM ReferenceMetal   PPH
    # 1 B18-10203 34.44325 -120.4297 Inner Shelf    Arsenic 6.28           Iron 0.893
    # 2 B18-10206 34.39825 -119.8652   Mid Shelf    Arsenic 4.58           Iron 0.761
    # 3 B18-10208 34.33409 -119.4346 Inner Shelf    Arsenic 8.25           Iron 1.500
    # 4 B18-10209 34.28353 -119.3544 Inner Shelf    Arsenic 5.93           Iron 1.050
    # 5 B18-10210 34.24356 -119.3853 Inner Shelf    Arsenic 5.53           Iron 1.360
    # 6 B18-10211 34.22837 -119.3520 Inner Shelf    Arsenic 5.14           Iron 1.430
    
    #rm_column:
    # [1] "ReferenceMetal"
    
    #rm_value:
    # [1] "PPH"
    
    #tm_column:
    # [1] "TraceMetal"
    
    #tm_value:
    # [1] "PPM"
    
    ####################################################################################
    # function_code
    ####################################################################################
    
    
    # only look at one reference metal
    normalized_data = normalized_data[normalized_data[[rm_column]]==reference_metal,]
    
    # we're going to iterate through the trace metals, and perform the same analysis on each
    # !! need to update pred_interval = predict(normal.model[[i]],newdata = normal.data[[i]], interval="prediction",level = 0.99)
    
    names = unique(normalized_data[[tm_column]])
    n = length(names)
    samplesize = rep(0,n)
    min = c(rep(0,n))
    max = c(rep(0,n))
    r2 = c(rep(0,n))
    sigma = c(rep(0,n))
    slope = c(rep(0,n))
    intercept = c(rep(0,n))
    metal_sum_of_squares = c(rep(0,n))
    mean_metal = c(rep(0,n))
    var_slope = c(rep(0,n))
    var_intercept = c(rep(0,n))
    mse = c(rep(0,n))
    fstat = c(rep(0,n))
    
    for(i in 1:n){
      subset_data = normalized_data[normalized_data[[tm_column]]==names[i],]
      samplesize[i] = nrow(subset_data)
      formula = paste(tm_value,"~",rm_value)
      normal_model = summary(lm(formula,data=subset_data))
      min[i] = min(subset_data[[tm_value]])
      max[i] = max(subset_data[[tm_value]])
      var_slope[i] = coefficients(normal_model)[2,2]^2
      var_intercept[i] = coefficients(normal_model)[1,2]^2
      metal_sum_of_squares[i] = sum((subset_data[[rm_value]] - mean(subset_data[[rm_value]]))^2)
      mean_metal[i] = mean(subset_data[[rm_value]])
      slope[i] = round(normal_model$coefficients[2,1],3)
      intercept[i] = round( normal_model$coefficients[1,1],3)
      sigma[i] = 2.576*round(normal_model$sigma,3)
      
      r2[i] = round(normal_model$r.squared,3)
      
      mse[i] = mean(normal_model$residuals^2)
      fstat[i] = normal_model$fstatistic[1]
    }
    
    df = as.data.frame(cbind(names,samplesize,min,max,r2,slope,intercept,sigma))

    names(df) = c("Reference Metal(% dry) versus","Sample Size","Minimum","Maximum",
                  "$r^2$","Slope  (m)","intercept (b)","Prediction Interval")
    
    return_kable = kable(df, format ="html",booktabs = T,escape=F,align="c") %>%
      kable_styling(position = "center") %>%
      column_spec(1:8, width = "5cm")
    
    
    df$mean_metal = mean_metal
    df$rm_sum_of_squares = metal_sum_of_squares
    df$resi_var = (sigma/2.57)^2
    df$var_slope = var_slope
    df$var_intercept = var_intercept
    table_3 <- read_csv("Data/table_3.csv")
    df = merge(df,table_3,by.x="Reference Metal(% dry) versus",by.y = "reference")

    df$se_1999 = (df$ci_1999/(2.57*sqrt(1+(1/df$n_1999))))^2
    df$sse_1999 = df$se_1999*(df$n_1999-2)
    df$sst_1999 = df$sse_1999/(1-df$r2_1999)
    df$r_1999 = sqrt(df$r2_1999)
    df$s_y_1999 = sqrt(df$sst/(df$n_1999 - 1 ))
    df$s_x_1999 = (df$r_1999 * df$s_y_1999)/df$slope_1999
    df$ss_1999 = (df$s_x_1999)^2*(df$n_1999-1)
    #######################
    # I believe this should be square
    #######################
    df$var_slope_1999 = (sqrt(df$se_1999/df$ss_1999))^2
    df$var_intercept_1999 = (sqrt(df$se_1999)*(1/df$n_1999 + df$mean_metal^2/df$ss_1999))^2
    df$slope = as.numeric(as.character(df[["Slope  (m)"]]))
    df$intercept = as.numeric(as.character(df[["intercept (b)"]]))
    df$slope_equal_z_value = (df$slope - df$slope_1999)/sqrt(df$var_slope+df$var_slope_1999)
    df$intercept_equal_z_value = (df$intercept - df$intercept_1999)/sqrt(df$var_intercept+df$var_intercept_1999)
    
    df$slope_p_val = 2*(1-pnorm(abs(df$slope_equal_z_value)))
    df$intercept_p_val = 2*(1-pnorm(abs(df$intercept_equal_z_value)))
    
    
    
    
     df = df[c("Reference Metal(% dry) versus","intercept_1999", "intercept", "var_intercept_1999","var_intercept",
              "slope_1999","slope", "var_slope_1999","var_slope","intercept_p_val","slope_p_val")]
    
  names(df) = c("Trace Metal","Intercept (1999)","Intercept (2018)", "Intercept Variance (1999)", "Intercept variance (2018)",
   "Slope (1999)","Slope (2018)","Slope Variance (1999)","Slope Variance (2018)",  "P-Value of Intercepts","P-Value of Slopes")
    
    
      
    
    slope_test_kable = kable(df, format ="html",booktabs = T,escape=F,align="c",digits=4) %>%
      kable_styling(position = "center") %>%
      column_spec(1:11, width = "5cm")

          
    if(slope_test){
      if(reference_metal == "Iron"){

      return(slope_test_kable)
      }
      else{
        df=as.data.frame("No P-value available")
        names(df) = reference_metal
        return(kable(df))
      }
    }
    else{
    return(return_kable)
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  normalizer_summary_kable = function(normalized_data,rm_column,rm_value,tm_column,tm_value,calibration_sites = F,return_kable_and_df = F){
    names = c(unique(normalized_data[[tm_column]]))
    n = length(names)
    iron_sample_size = rep(0,n)
    iron_r2 = rep(0,n)
    iron_pval = rep(0,n)
    
    aluminum_sample_size = rep(0,n)
    aluminum_r2 = rep(0,n)
    aluminum_pval = rep(0,n)
    
    grains_sample_size = rep(0,n)
    grains_r2 = rep(0,n)
    grains_pval = rep(0,n)
    
    for (i in 1:n){
      subset_data = normalized_data[normalized_data[[tm_column]]==names[i],]
      iron_subset = subset_data[subset_data[[rm_column]]=="Iron",]
      
      aluminum_subset = subset_data[subset_data[[rm_column]]=="Aluminum",]
      grains_subset = subset_data[subset_data[[rm_column]]=="GrainSize",]
      formula = paste(tm_value,"~",rm_value)
      
      normal_model_iron = summary(lm(formula,data=iron_subset))
      iron_sample_size[i] = nrow(iron_subset)
      iron_r2[i] = round(normal_model_iron$r.squared,3)
      iron_pval[i] = formatC(normal_model_iron$coefficients[2,4],format = "e", digits = 2)
      
      if (i==3) {
        aluminum_sample_size[i] = '-'
        aluminum_r2[i] = '-'
        aluminum_pval[i] = '-'
      }
      else {
        normal_model_aluminum = summary(lm(formula,data=aluminum_subset))
        aluminum_sample_size[i] = nrow(aluminum_subset)
        aluminum_r2[i] = round(normal_model_aluminum$r.squared,3)
        aluminum_pval[i] = formatC(normal_model_aluminum$coefficients[2,4],format = "e", digits = 2)
      }
      
      normal_model_grains = summary(lm(formula,data=grains_subset))
      grains_sample_size[i] = nrow(grains_subset)
      grains_r2[i] = round(normal_model_grains$r.squared,3)
      grains_pval[i] = formatC(normal_model_grains$coefficients[2,4],format = "e", digits = 2)
    }
    
    df = as.data.frame(cbind(names,iron_sample_size,aluminum_sample_size,grains_sample_size,iron_r2,aluminum_r2,grains_r2,iron_pval,aluminum_pval,grains_pval))
    names(df) = c('Trace Metal',rep(c('Iron','Aluminum','Fines'),3))
    return_kable = kable(df,format="html",booktabs=T,escape=F,align="c") %>%
      kable_styling(position = "center") %>%
      add_header_above(c(" ", "Sample Size" = 3, "r-squared" = 3, "p-value" = 3))
    
    
  }
})






with(calibration.env,{
  
  
  
  
  
  
  
  
  
  

  # function calibration_plot generates the necessary information about the
  # calibration data of a trace metal against a contaminant
  calibration_plot = function(reference_metal,trace_metal,clean_sites,dirty_sites,calibration_sites=F){
    ####################################################################################
    # data expected/provided to function
    ####################################################################################
    
    #reference metal = string like "Iron" or "Aluminum"
    #trace metal = string, like "Cadmium", or "Copper"

    #clean_sites = dataframe that looks like this:
    # stationid      lat      long     Stratum TraceMetal  PPM ReferenceMetal   PPH
    # 2 B18-10206 34.39825 -119.8652   Mid Shelf    Arsenic 4.58           Iron 0.761
    # 4 B18-10209 34.28353 -119.3544 Inner Shelf    Arsenic 5.93           Iron 1.050
    # 5 B18-10210 34.24356 -119.3853 Inner Shelf    Arsenic 5.53           Iron 1.360
    # 6 B18-10211 34.22837 -119.3520 Inner Shelf    Arsenic 5.14           Iron 1.430
    # 7 B18-10212 34.19924 -119.2964 Inner Shelf    Arsenic 4.60           Iron 1.290
    # 8 B18-10213 34.17873 -119.3468 Inner Shelf    Arsenic 4.94           Iron 1.490
    
    # dirty_sites = dataframe that looks identical to clean_sites, with known dirty sites
    
    # calibration_sites, boolean that says whether to return plot with calibration sites 
    # overlaid, or with dirty sites overlaid
    
    
    ####################################################################################
    # filter data
    ####################################################################################
    
    # clean_sites and dirty_sites is an entire dataset. only grab the metals we need
    clean_sites = filter(clean_sites,str_detect(TraceMetal, trace_metal) ==TRUE) %>%
      filter(str_detect(ReferenceMetal, reference_metal) ==TRUE)
    dirty_sites = filter(dirty_sites,str_detect(TraceMetal, trace_metal) ==TRUE) %>%
      filter(str_detect(ReferenceMetal, reference_metal) ==TRUE)
    
    
    
    ####################################################################################
    # generate calibration models
    ####################################################################################
    
    formula = paste0("PPM", " ~ ", "PPH")
    model = lm(formula,data=clean_sites)
    
    if(calibration_sites){
      dirty_sites = clean_sites
    }
    
    predicted_PPM = as.data.frame(predict(model,newdata=dirty_sites,interval="prediction"))
    predicted_PPM = cbind(dirty_sites,predicted_PPM)
    predicted_PPM$Actual = dirty_sites$PPM
    predicted_PPM$Residual = predicted_PPM$Actual - predicted_PPM$fit
    palette = brewer.pal(n = 11, "Spectral")
    
    predicted_PPM$Interval = (predicted_PPM$Actual <= predicted_PPM$upr) #& predicted_PPM$Actual >= predicted_PPM$lwr)
    predicted_PPM$distance=predicted_PPM$Actual-predicted_PPM$upr
    predicted_PPM$stationid=dirty_sites$stationid
    predicted_PPM$Human_Addition=ifelse(predicted_PPM$distance<0,0, predicted_PPM$distance)
     
         # Tanya's beautiful plot code
         return_pointsPlot = ggplot(dirty_sites, aes_string(x="PPH", y="PPM"))    +
           geom_point()+
           geom_smooth(method = "lm", fullrange=TRUE, se = TRUE, color = "black",data=clean_sites)+
           theme_bw()+
    geom_line(aes(y = predicted_PPM$lwr), color = "#9E0142", linetype = "dashed")+
      geom_line(aes(y = predicted_PPM$upr), color = "#9E0142", linetype = "dashed")+
      geom_point(aes(colour = predicted_PPM$Interval))+
      scale_colour_manual(name = 'Within Pred. Interval', values = setNames(c("#000000","#3288BD"),c(T, F)))+
           xlab(reference_metal)+
           ylab(trace_metal)+
           ggtitle("Trace Metal-Reference Metal Plots Overlaid With Reference Element Baseline Relationships")
    
         
      
       #####
         
         
         predicted_PPM2 = as.data.frame(predict(model,newdata=clean_sites,interval="prediction"))
         predicted_PPM2 = cbind(clean_sites,predicted_PPM2)
         predicted_PPM2$Actual = clean_sites$PPM
         predicted_PPM2$Residual = predicted_PPM2$Actual - predicted_PPM2$fit
         palette = brewer.pal(n = 11, "Spectral")
         
         predicted_PPM2$Interval = (predicted_PPM2$Actual <= predicted_PPM2$upr & predicted_PPM2$Actual >= predicted_PPM2$lwr)
         
         
         predicted_PPM2$lwrResid = predicted_PPM2$lwr - predicted_PPM2$fit
         predicted_PPM2$uprResid = predicted_PPM2$upr - predicted_PPM2$fit
         
         
         residualsPlot = ggplot(predicted_PPM2,aes(PPH,Residual)) +
           geom_point() +
           geom_smooth(method = lm,se=F)+
           geom_line(aes(y = predicted_PPM2$lwrResid), color = "#9E0142", linetype = "dashed")+
           geom_line(aes(y = predicted_PPM2$uprResid), color = "#9E0142",linetype = "dashed")+
           labs(x = reference_metal, y = "Residuals",title = "Residuals of Calibration")
         
         
         #### Testing Residual Prediction Intervals
         residual_model_lat = lm(Residual~lat,data=predicted_PPM2)
         predicted_residual_lat = as.data.frame(predict(residual_model_lat,newdata=predicted_PPM2,interval="prediction"))    
         
         residual_model_long = lm(Residual~long,data=predicted_PPM2)
         predicted_residual_long = as.data.frame(predict(residual_model_long,newdata=predicted_PPM2,interval="prediction"))      
         
         residual_model_depth = lm(Residual~depth,data=predicted_PPM2)
         predicted_residual_depth = as.data.frame(predict(residual_model_depth,newdata=predicted_PPM2,interval="prediction"))      
         
         residualsByLat =  ggplot(predicted_PPM2,aes(lat,Residual)) +
           geom_point() +
           geom_smooth(method = lm,se=F)+
           geom_line(aes(y = predicted_residual_lat$lwr), color = "#9E0142", linetype = "dashed")+
           geom_line(aes(y = predicted_residual_lat$upr), color = "#9E0142", linetype = "dashed")+
           labs(x = "Lattitude", y = "Residuals",title = "Residuals By Latitude") 
         
         residualsByLong =  ggplot(predicted_PPM2,aes(long,Residual)) +
           geom_point() +
           geom_smooth(method = lm,se=F)+
           geom_line(aes(y = predicted_residual_long$lwr), color = "#9E0142", linetype = "dashed")+
           geom_line(aes(y = predicted_residual_long$upr), color = "#9E0142", linetype = "dashed")+
           labs(x = "Lattitude", y = "Residuals",title = "Residuals By Longitude") 
    
         # New plot for Depth Vs Residuals
    
         residualsByDepth =  ggplot(predicted_PPM2,aes(depth,Residual)) +
         geom_point() +
         geom_smooth(method = lm,se=F)+
           geom_line(aes(y = predicted_residual_depth$lwr), color = "#9E0142", linetype = "dashed")+
           geom_line(aes(y = predicted_residual_depth$upr), color = "#9E0142", linetype = "dashed")+
         labs(x = "Depth", y = "Residuals",title = "Residuals By Depth") 
         
         APByLat =  ggplot(predicted_PPM2,aes(lat,Actual/fit)) +
           geom_point() +
           geom_smooth(method = lm,se=F)+
           labs(x = "Latitude", y = "Correlation",title = "Ratio Actual to Predicted By Latitude")

         APByLong =  ggplot(predicted_PPM2,aes(long,Actual/fit)) +
           geom_point() +
           geom_smooth(method = lm,se=F)+
           labs(x = "Longitude", y = "Correlation",title = "Ratio Actual to Predicted By Longitude")
    
  
        
     StratumPlot = ggplot(predicted_PPM, mapping = aes(x = reorder(Stratum, Human_Addition), Human_Addition,fill=Stratum))+ 
          stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
          #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")+
          geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
          labs(x = "Stratum", y = trace_metal,title = "Mean Anthropogenic Concentration By Stratum")
         
         
         
         
    return_data = list()
    return_data[["clean_sites"]] = clean_sites
    return_data[["dirty_sites"]] = dirty_sites
    return_data[["predicted_PPM"]] = predicted_PPM
    return_data[["pointsPlot"]] = return_pointsPlot
    return_data[["residualsPlot"]] = residualsPlot
    return_data[["residualsByLat"]] = residualsByLat
    return_data[["residualsByLong"]] = residualsByLong
    return_data[["residualsByDepth"]] = residualsByDepth
    return_data[["APByLat"]] = APByLat
    return_data[["APByLong"]] = APByLong
    return_data[["model"]] = model
     return_data[["StratumPlot"]] = StratumPlot
    return(return_data)
    
  }
  
  #test function
  
  
})























#############################################
#                                           #
# setup calibration environment             #
#                                           #  
#############################################




with(calibration.env,{
  dataset = ETL.Load$calibration_data
  dataset$depth[dataset$depth < 0] = 0
  trace_metals = c("Arsenic","Cadmium","Chromium","Copper","Lead","Nickel","Silver","Zinc")
  reference_metals = c("Iron","Aluminum","GrainSize")
  #convert reference metals to pct
  for (reference_metal in reference_metals){
    dataset[[reference_metal]] = dataset[[reference_metal]]/10000
  }
  
  
  # for this analysis, we simply need the trace metals and the reference metals.
  
  dataset = subset(dataset,select=c("stationid","lat","long","depth",reference_metals,trace_metals,"is_contaminated_site","Stratum"))
  # gather the contaminants so we can perform better analysis
  
  dataset = dataset %>% gather("TraceMetal","PPM",trace_metals)
  
  # gather the trace metals
  
  dataset = dataset %>% gather("ReferenceMetal","PPH",reference_metals)
  
  clean_sites = dataset[dataset$is_contaminated_site==0,] %>%
    subset(select=-c(is_contaminated_site)) %>%
    filter(str_detect(Stratum, "Shelf") ==TRUE)
  dirty_sites = dataset[dataset$is_contaminated_site!=0,]  %>%
    subset(select=-c(is_contaminated_site))
  
  
  # [[normalize residuals]]
  sites_normalized = normalize_residuals(clean_sites,"ReferenceMetal","PPH","TraceMetal","PPM")
  clean_sites_normalized=  sites_normalized$clean_sites
  dirty_sites = rbind(dirty_sites, sites_normalized$new_contaminated_sites)
  
})





## Generate a summary table 

with(calibration.env,{
  
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
  
  
  

  B18_station = merge(ETL.Extract$station_completion.trawls,ETL.Extract$chemistry_bight.tbl_chemresults,by="stationid") %>%
    merge(ETL.Extract$station_completion.grabs,by="stationid") %>%
    rename(chemareaweight = `Area Weight All Sites`,stratum = stratum.x) 
  B18_station$result = as.numeric(B18_station$result)
  
  
  
  
  
  
  
  weightedSedimentSummaryStatistics=calculate_AWM_bight18(B18_station)
  
  
  
  
  
})




