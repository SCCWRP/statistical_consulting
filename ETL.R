
# 
# ETL.Extract = new.env()
# ETL.Transform = new.env()
# ETL.Load = new.env()
# calibration.env = new.env()
# 
# 
# 
# 
# 
# 
# with(ETL.Extract,{
#   client_data_prepend = "Data/"
# 
# 
#   grain_size_function = function(client_data_prepend){
# 
#     chemistry_bight_path = paste(client_data_prepend,"CHEMISTRY-BIGHT18.xlsx",sep="")
# 
# 
# 
#     chemical_results <- as.data.frame(read_excel(chemistry_bight_path,sheet='tbl_chemresults',col_types = c("text",
#                                                                                                             "text", "text", "text", "text",
#                                                                                                             "text", "text", "text", "text", "text",
#                                                                                                             "text", "text", "text", "text",
#                                                                                                             "text", "text", "text", "text",
#                                                                                                             "text", "text", "text", "text",
#                                                                                                             "text", "text", "text", "text",
#                                                                                                             "text", "text", "text", "text")))
#     chemical_results = chemical_results[chemical_results$sampletype=="Result",]
#     chemical_results = subset(chemical_results, select=-c(sampletype))
# 
#     chemical_results.grainsize = chemical_results[chemical_results$analyteclass=="GrainSize",]
#     chemical_results.grainsize = subset(chemical_results.grainsize,select=-c(analyteclass))
# 
#     chemical_results.grainsize$saID = paste(chemical_results.grainsize$stationid,chemical_results.grainsize$analytename)
# 
#     chemical_results.grainsize = subset(chemical_results.grainsize,select=-c(bioaccumulationsampleid,
#                                                                              preparationbatchid,
#                                                                              analysisbatchid,
#                                                                              matrix,
#                                                                              analysismethod,
#                                                                              labsampleid,
#                                                                              truevalue,
#                                                                              percentrecovery,sampledate,created_date,qualifier,mdl,rl,lab,analysisdate,globalid,created_user,last_edited_user,last_edited_date,qacode,comments))
# 
#     chemical_results.grainsize = chemical_results.grainsize[order(chemical_results.grainsize$fieldduplicate, decreasing=TRUE),]
#     chemical_results.grainsize = chemical_results.grainsize[!duplicated(chemical_results.grainsize$saID),]
# 
#     chemical_results.grainsize = chemical_results.grainsize[order(chemical_results.grainsize$labreplicate, decreasing=TRUE),]
#     chemical_results.grainsize = chemical_results.grainsize[!duplicated(chemical_results.grainsize$saID),]
# 
#     chemical_results.grainsize = subset(chemical_results.grainsize,select=-c(fieldduplicate,labreplicate,saID,objectid,units))
# 
#     chemical_results.grainsize$result = as.numeric(chemical_results.grainsize$result)
#     chemical_results.grainsize$result[chemical_results.grainsize$result==-88.0] =0
# 
#     chemical_results.grainsize$analytename = as.numeric(gsub(pattern='Phi ',replacement='',x=chemical_results.grainsize$analytename))
#     chemical_results.grainsize = chemical_results.grainsize[chemical_results.grainsize$analytename >= 4.5,]
# 
#     grainsize_totals <- chemical_results.grainsize %>%
#       spread(key = analytename, value = result)
#     grainsize_totals[is.na(grainsize_totals)] <- 0
# 
#     grainsize_totals = data.frame(stationid = grainsize_totals$stationid,GrainSize = rowSums(grainsize_totals[,-1]))
# 
#     return(grainsize_totals)
#   }
# 
# })
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Convert excel sheets to their own CSV files. Could do this automatically, except the data is formatted so it crashes readxl, and readxl is too stupid to understand how to handle the data.
# with(ETL.Extract,{
# 
# 
#   chemistry_bight_path = paste(client_data_prepend,"CHEMISTRY-BIGHT18.xlsx",sep="")
#   ref_phi_conversion_path = paste(client_data_prepend,"Ref - Phi Conversion.xlsx",sep="")
#   station_completion_path = paste(client_data_prepend,"StationCompletionV8.xlsx",sep="")
#   stations_list_path = paste(client_data_prepend,"B18 stations list.csv",sep="")
# 
# 
#   grain_size = grain_size_function(client_data_prepend)
# 
#   chemistry_bight.tbl_chembatch = read_excel(chemistry_bight_path,
#                                              sheet = "tbl_chembatch", col_types = c("text",
#                                                                                     "text", "text", "text", "date", "text",
#                                                                                     "text", "text", "text", "text"))
# 
# 
# 
#   chemistry_bight.tbl_chemresults = read_excel(chemistry_bight_path,
#                                                sheet = "tbl_chemresults", col_types = c("text",
#                                                                                         "text", "text", "text", "text",
#                                                                                         "text", "text", "text", "text", "text",
#                                                                                         "text", "text", "text", "text",
#                                                                                         "text", "text", "text", "text",
#                                                                                         "text", "text", "text", "text",
#                                                                                         "text", "text", "text", "text",
#                                                                                         "text", "text", "text", "text")) # hammer against forehad
# 
#   chemistry_bight.tbl_sample_assignment_table = read_excel(chemistry_bight_path,
#                                                            sheet = "sample_assignment_table")
# 
# 
#   ref_phi_conversion.ref_phi_conversion = read_excel(ref_phi_conversion_path,sheet="Ref___Phi_Conversion")
# 
#   ref_phi_conversion.for_r = read_excel(ref_phi_conversion_path,sheet="for R")
# 
#   ref_phi_conversion.sheet_2 = read_excel(ref_phi_conversion_path,sheet="Sheet2")
# 
#   station_completion.grabs = read_excel(station_completion_path,sheet="Grabs")
#   station_completion.sheet_1 = read_excel(station_completion_path,sheet="Sheet1")
#   station_completion.trawls = read_excel(station_completion_path,sheet="Trawls")
#   station_completion.areas= read_excel(station_completion_path,sheet="areas")
# 
#   b18_stations_list = read_csv(stations_list_path)
# 
# })
# 
# 
# with(ETL.Transform,{
# 
#   # table containing results of different tests performed on samples
#   chemical_tests = ETL.Extract$chemistry_bight.tbl_chemresults
# 
#   # table that helps us convert phi values into a summary statistic for particle size in a
#   # sample
#   phi_conversion_values = ETL.Extract$ref_phi_conversion.for_r
# 
#   # table that contains length/area information about various regions
#   # might need to reformat to get more general information about bays/estuaries/etc.
#   region_information = ETL.Extract$station_completion.areas
# 
#   # information about specific stations where samples come from
#   station_information = ETL.Extract$b18_stations_list
# })
# 
# with(ETL.Transform,{
# 
#   # Only consider test results
#   chemical_results = chemical_tests[chemical_tests$sampletype=="Result",]
#   chemical_results = subset(chemical_results, select=-c(sampletype))
# 
#   # only look at inorganic analyte class
#   chemical_results.inorganics = chemical_results[chemical_results$analyteclass=="Inorganics",]
#   chemical_results.inorganics = subset(chemical_results.inorganics, select=-c(analyteclass))
# 
#   # set a unique id for each station and chemical
#   chemical_results.inorganics$saID = paste(chemical_results.inorganics$stationid,chemical_results.inorganics$analytename)
# 
#   # remove the values that are not needed
#   chemical_results.sans_unnecessary_ids = subset(chemical_results.inorganics,select=-c(bioaccumulationsampleid,
#                                                                                        preparationbatchid,
#                                                                                        analysisbatchid,
#                                                                                        matrix,
#                                                                                        analysismethod,
#                                                                                        labsampleid,
#                                                                                        truevalue,
#                                                                                        percentrecovery))
# 
#   # move the following columns to another table (in case they are needed in the future)
#   chemical_results.ancilliary_data = subset(chemical_results.sans_unnecessary_ids,select=c(objectid,sampledate,created_date,qualifier,mdl,rl,lab,analysisdate,globalid,created_user,last_edited_user,last_edited_date,qacode,comments))
# 
#   chemical_results.without_ancillary_data = subset(chemical_results.sans_unnecessary_ids,select=-c(sampledate,created_date,qualifier,mdl,rl,lab,analysisdate,globalid,created_user,last_edited_user,last_edited_date,qacode,comments))
# 
#   # per e-mail, only keep the most recent fieldduplicate and labreplicate.
#   chemical_results.only_most_recent_fieldduplicate1 = chemical_results.without_ancillary_data[order(chemical_results.without_ancillary_data$fieldduplicate, decreasing=TRUE),]
#   chemical_results.only_most_recent_fieldduplicate2 = chemical_results.only_most_recent_fieldduplicate1[!duplicated(chemical_results.only_most_recent_fieldduplicate1$saID),]
# 
# 
#   chemical_results.only_most_recent_labreplicate1 = chemical_results.only_most_recent_fieldduplicate2[order(chemical_results.only_most_recent_fieldduplicate2$labreplicate, decreasing=TRUE),]
#   chemical_results.only_most_recent_labreplicate2 = chemical_results.only_most_recent_labreplicate1[!duplicated(chemical_results.only_most_recent_labreplicate1$saID),]
# 
# 
#   chemical_results.without_duplicates_ids = subset(chemical_results.only_most_recent_labreplicate2,select=-c(fieldduplicate,labreplicate,saID,objectid,units))
# 
#   # lastly, any results where we get a -88 cannot detect, we'll assume to be 0.
#   chemical_results.without_duplicates_ids$result = as.numeric(chemical_results.without_duplicates_ids$result)
#   chemical_results.without_duplicates_ids$result[chemical_results.without_duplicates_ids$result==-88.0] =0
# 
# 
#   chemical_results.final =chemical_results.without_duplicates_ids %>% spread(analytename,result)
# })
# 
# 
# with(ETL.Transform,{
#   station_information$is_contaminated_site = 0
# 
#   station_information$is_contaminated_site[!is.na(station_information$River_Mout)] = 1
#   station_information$is_contaminated_site[!is.na(station_information$POTW_Statu)] = 1
#   station_information$is_contaminated_site[!is.na(station_information$POTW_Name)] = 1
# 
# 
#   station_information.without_irrelevant_columns = subset(station_information,select=-c(River_Mout,POTW_Statu,POTW_Name,OBJECTID))
# 
#   station_information.final = station_information.without_irrelevant_columns
# 
# })
# 
# 
# with(ETL.Transform,{
# 
#   calibration_data_wo_grain = merge(x=chemical_results.final,y=station_information.final,by="stationid")
#   calibration_data = merge(x=calibration_data_wo_grain,y=ETL.Extract$grain_size,by="stationid")
# 
# })
# 
# 
# save.image("etl.Rdata")

sys.load.image("etl.Rdata",quiet=F)


