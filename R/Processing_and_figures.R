if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  rnaturalearth, rnaturalearthdata,
  osmdata, sf, raster, ggplot2, data.table, magrittr, ggforce, lwgeom,
  ggrepel, gridExtra, tidyr, dplyr, RColorBrewer,
  gganimate, transformr, gifski, viridis, scales,
  Hmisc, parallel, patchwork
)

Sys.setlocale("LC_ALL", "en_US.UTF-8")

path_proj = "" # project path

# (PRE) PROCESSING
process_DATA_FORMAT = F
process_SNOW_DATA = F
process_BREAKAGE_DATA = F

# PROCESSING GAYA 2.0/JLP22 RESULTS
process_GAYA_SCHEDULES = F
process_GAYA_SCHEDULES_CALIBRATE_MODEL_PARAM = F
process_GAYA_SCHEDULES_SUBSTAND_LEVEL = F
process_GAYA_SCHEDULES_SUBSTAND_LEVEL_ADD_NEW = F
process_GAYA_SCHEDULES_STAND_LEVEL = F
process_GAYA_SCHEDULES_PROPERTY_LEVEL = F

# FIGURES
figure_STUDY_AREA = F # FIGURE 1
figure_INITIAL_STATE = F # FIGURE 2
figure_CYCLES = F # FIGURE 3
figure_PRODUCTION_TABLE_TIMESERIES = F # FIGURE 4
figure_SNOW_VARIABLES_TIMESERIES = F # FIGURE 5
figure_CSL_SDP_BY_SI = F # FIGURE 6
figure_CSL_SDP_HISTO_5_10 = F # FIGURE 7
figure_CSL_vs_BREAKP = F # FIGURE 8

figure_DIFF = F # GRAPHICAL ABSTRACT
figure_DIFF_CSL_SDP = F # FIGURE 9

# global variables and functions
interest_rate = 3
period_length = 10
n_per = 15
treatment_offset = 5
period_later = 5
period_last = 10
period_last_value = 15
eur_to_nok = 11.6828
nok_to_eur = 1/eur_to_nok

onecol = 90
twocol = 190
maxHeight = 240

breaking_prob = function(asl, csl, calibration_param, scale_factor){
  result = 1/(1 + exp(-(asl - csl + calibration_param)/scale_factor))
  result[is.nan(result)] = 0
  return(result)  
}

multiyear_prob = function(p, n){
  1 - (1 - p)^n
}

if (process_DATA_FORMAT){
  raw_data = 
    read_excel(paste0(path_root,"\\Data\\fritso_gaya_tilrettelagt_23022024.xlsx")) %>%
    data.table()
  
  calculate_mainSp <- function(ba1, ba2, ba3, n1, n2, n3) {
    if (ba1 > 0 && ba1 > ba2 && ba1 > ba3) {
      return(1)
    } else if (ba2 > 0 && ba2 > ba1 && ba2 > ba3) {
      return(2)
    } else if (ba3 > 0 && ba3 > ba1 && ba3 > ba2) {
      return(3)
    } else if (ba1 == 0 && ba2 == 0 && ba3 == 0) {
      if (n1 >= n2 && n1 >= n3) {
        return(1)
      } else if (n2 >= n1 && n2 >= n3) {
        return(2)
      } else if (n3 >= n1 && n3 >= n2) {
        return(3)
      }
    }
    return(1) 
  }
  
  raw_data[, mainSp := apply(.SD, 1, function(x) calculate_mainSp(x["Ghag"], x["Ghaf"], x["Ghal"], x["Nhag"], x["nhaf"], x["Nhal"])), .SDcols = c("Ghag","Ghaf","Ghal","Nhag","nhaf","Nhal")]
  
  raw_data = raw_data[foresttype==1]
  
  
  formated_data = data.table(
    standID = (raw_data$Standno*100+raw_data$Standpart)*10+raw_data$district,
    area = raw_data$Arealha,
    mainSp = raw_data$mainSp, 
    dev_class = as.numeric(NA),
    `sp(1)` = 1,
    `sp(2)` = 2,
    `sp(3)` = 3,
    `si(1)` = raw_data$SI_G,
    `si(2)` = raw_data$SI_F,
    `si(3)` = raw_data$SI_L,
    `age(1)` = raw_data$yrtot,
    `age(2)` = raw_data$yrtot,
    `age(3)` = raw_data$yrtot,
    `age13(1)` = raw_data$t13g,
    `age13(2)` = raw_data$t13f,
    `age13(3)` = raw_data$t13l,
    `n(1)` = raw_data$Nhag,
    `n(2)` = raw_data$nhaf,
    `n(3)` = raw_data$Nhal,
    `hDom(1)` = raw_data$HO_G,
    `hDom(2)` = raw_data$HO_F,
    `hDom(3)` = raw_data$HO_L,
    `ba(1)` = raw_data$Ghag,
    `ba(2)` = raw_data$Ghaf,
    `ba(3)` = raw_data$Ghal,
    defRegen = 0,
    vegetation_type = raw_data$Vegtype,
    altitude = raw_data$alt,
    slope = raw_data$slope,
    driveDist = raw_data$rdist*100,
    extraCutCostThinning = raw_data$add_tf/100,
    extraCutCostHarvest = raw_data$add_hf/100,
    extraTransportCostThinning = raw_data$add_tt/100,
    extraTransportCostHarvest = raw_data$add_ht/100,
    extPulp = raw_data$add_pulp/100,
    extra_tending_cost = raw_data$add_tending/100
  )
  
  formated_data[`n(1)`+`n(2)`+`n(3)`==0,]$dev_class = 1
  formated_data[`n(1)`+`n(2)`+`n(3)`>0&`ba(1)`+`ba(2)`+`ba(3)`==0,]$dev_class = 2
  formated_data[`n(1)`+`n(2)`+`n(3)`>0&`ba(1)`+`ba(2)`+`ba(3)`>0,]$dev_class = 3
  
  # this is the stand inventory input file for GAYA 2.0
  fwrite(formated_data, paste0(path_root,"\\Data\\marcsman_stands.csv"), sep = ",", col.names = TRUE)
}

if (process_SNOW_DATA){
  path_snow_data = paste0(path_proj, "\\Data\\Snow\\SnC8-NMBU-20240312T165206Z-001\\SnC8-NMBU")
  
  tiff_to_dt = function(tiff_file){
    year_month = substr(tiff_file, nchar(tiff_file) - 9, nchar(tiff_file) - 4)
    snow_brick = brick(tiff_file)
    snow_dt = as.data.frame(snow_brick, xy = T) %>% 
      setDT() %>%
      melt(id.vars = c("x", "y"), variable.name = "band", value.name = "val") %>%
      .[band %like% "^Swet",.(max_sl = max(val)),by=.(x,y)] %>%
      .[,y_m := year_month]
    return(snow_dt)
  }
  
  tiff_files = list.files(path = path_snow_data, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
  
  max_sl_dt = lapply(tiff_files, tiff_to_dt) %>%
    rbindlist() %>%
    .[,.(max_sl = max(max_sl)),by=.(x,y)]
  
  first_snow_brick = brick(tiff_files[1])
  
  coordinates(max_sl_dt) = ~x+y 
  proj4string(max_sl_dt) = crs(first_snow_brick) 
  
  forest_stands_shape = st_read(paste0(path_proj,"\\Data\\Data_Fritsoe_Tron\\Bestand20t.shp")) %>%
    mutate(stand = (BestNr * 100 + DelNr) * 10 + Teig)
  
  max_sl_raster = rasterize(max_sl_dt, first_snow_brick, field = "max_sl", fun = max) %>%
    projectRaster(crs = st_crs(forest_stands_shape)$proj4string)
  
  names(max_sl_raster) = c("max_sl")
  
  #writeRaster(max_sl_raster, filename = paste0(path_proj,"\\Data\\max_sl_raster.tif"), format = "GTiff", overwrite = TRUE)
  
  max_sl_dt = as.data.frame(max_sl_raster, xy = T) %>% 
    setDT()
  
  forest_stands_shape$max_sl = extract(max_sl_raster, st_centroid(forest_stands_shape))
  
  forest_stands_max_sl = as.data.table(forest_stands_shape) %>%
    .[,.(max_sl = max(max_sl)),by=stand]
  
  fwrite(forest_stands_max_sl,paste0(path_proj,"\\Data\\max_sl_data.csv"))
}

if (process_BREAKAGE_DATA){
  h1_threshold = 9
  hdiff_threshold = c(0,1,2,3,4,5)
  tree_data = fread(paste0(path_proj, "\\Data\\tree_match_v1.csv")) %>% .[h1 > h1_threshold]
  break_prob_6yr = sapply(hdiff_threshold, FUN = function(x){nrow(tree_data[h2 < (h1 - x)]) / nrow(tree_data)})
  break_prob_1yr = 1 - ( 1 - break_prob_6yr)^(1/6)
  breaking_probability_data = data.table(probability = break_prob_1yr, hdiff = hdiff_threshold)
  fwrite(breaking_probability_data, paste0(path_proj,"\\Data\\breaking_probability_data.csv"))
}

if (process_GAYA_SCHEDULES){
  csl_min = 0.01
  csl_max = 500
  
  if (process_GAYA_SCHEDULES_SUBSTAND_LEVEL){
    results_standard = fread(paste0(path_proj,"\\Results\\PT_marcsman_standard09.csv")) %>% .[, management := "standard"]
    results_fritzoe = fread(paste0(path_proj,"\\Results\\PT_marcsman_fritzoe09.csv")) %>% .[, management := "fritzoe"]
    
    if (!exists("forest_stands_max_sl")) forest_stands_max_sl = fread(paste0(path_proj,"\\Data\\max_sl_data.csv"))
    
    discount = function(val, interest, years) {
      return (val / ((1 + interest / 100) ^ years))
    }
    
    results_substand_lvl = # since some stands are split by LJP22
      rbind(results_standard,results_fritzoe) %>% 
      .[period%in%1:n_per] %>%
      .[, csl_after := ifelse(csl_after == 0, NaN, ifelse(csl_after < csl_min, csl_min, ifelse(csl_after > csl_max, csl_max, csl_after)))] %>%
      .[, csl_start := ifelse(csl_start == 0, NaN, ifelse(csl_start < csl_min, csl_min, ifelse(csl_start > csl_max, csl_max, csl_start)))] %>%
      .[, treatment := case_when(
        (trt %in% 1:15)  ~ "Tending",
        (management == "standard" & trt %in% c(16:17,20:21)) ~ "Thinning",
        (management == "fritzoe" & trt %in% 18:19) ~ "Thinning",
        (management == "standard" & trt %in% 18:19) | (management == "fritzoe" & trt %in% 16:17) ~ "Harvest",
        TRUE ~ "Let grow"
      )]  %>%
      left_join(forest_stands_max_sl, by = "stand") %>%
      .[,net_dsc := discount(net,interest_rate,(period-1)*period_length+treatment_offset)]
    
    if (process_GAYA_SCHEDULES_CALIBRATE_MODEL_PARAM){
      spread = quantile(results_substand_lvl[period==1,max_sl-csl_start],probs=c(.05,.95),na.rm=T)
      scale_factor_opt = .2*(max(spread)-min(spread))  
      forest_stands_results_P1 = results_substand_lvl[period==1&!is.na(max_sl)&!is.na(csl_start),.(asl=max_sl, csl=csl_start)]
      
      if (!exists("breaking_probability_data")) breaking_probability_data = fread(paste0(path_proj,"\\Data\\breaking_probability_data.csv"))
      target_bp = breaking_probability_data[hdiff==2,probability]
      tolerance_bp = 0.0000001
      
      binary_search_calibration_param = function(calibration_param_min, calibration_param_max){
        calibration_param_middle = 0.5 * (calibration_param_min + calibration_param_max)
        forest_stands_results_P1[,bp:=breaking_prob(asl, csl, calibration_param_middle, scale_factor_opt)]
        mean_bp = mean(forest_stands_results_P1$bp)
        print(paste0("CP = ",calibration_param_middle," BP = ",mean_bp))
        if (abs(mean_bp-target_bp)<tolerance_bp) return(calibration_param_middle)
        if (mean_bp>target_bp) 
          return(binary_search_calibration_param(calibration_param_min, calibration_param_middle))
        else 
          return(binary_search_calibration_param(calibration_param_middle, calibration_param_max))
      }
      
      calibration_param_opt = binary_search_calibration_param(-1000,1000)
      logistic_parameters = data.table(param=c("scale_factor","calibration_param"),
                                       value=c(scale_factor_opt, calibration_param_opt))
      fwrite(logistic_parameters,paste0(path_proj,"\\Data\\logistic_parameters.csv"))
    } 
    
    if (!exists("logistic_parameters")) logistic_parameters = fread(paste0(path_proj,"\\Data\\logistic_parameters.csv"))
    
    results_substand_lvl = 
      results_substand_lvl[,break_p_before:=breaking_prob(asl = max_sl, csl_start, logistic_parameters[param=="calibration_param",value], logistic_parameters[param=="scale_factor",value])] %>%
      .[,break_p_after:=breaking_prob(asl = max_sl, csl_after, logistic_parameters[param=="calibration_param",value], logistic_parameters[param=="scale_factor",value])] 
    
    results_substand_lvl_split = split(results_substand_lvl, by = c("stand", "schedule", "management"))
    
    calculate_probabilities = function(schedule_data){
      per = 1
      npv_at_risk = 0
      
      prob_acc = 0
      prob_after_previous = 0
      
      schedule_data$prob_acc = 0
      
      while (per <= n_per){
        prob_before = multiyear_prob(schedule_data[period == per, break_p_before], 5)
        prob_after = multiyear_prob(schedule_data[period == per, break_p_after], 5)
        
        prob_tot = 1 - (1 - prob_after_previous)*(1 - prob_before)   
        prob_acc = 1 - (1 - prob_acc)*(1 - prob_tot)
        
        schedule_data$prob_acc[per] = prob_acc
        
        if (schedule_data[period == per, treatment] %in% c("Harvest", "Thinning"))
          npv_at_risk = npv_at_risk + prob_acc * schedule_data[period == per, net_dsc]
        
        prob_after_previous = prob_after
        
        if (schedule_data[period == per, treatment] == "Harvest"){
          prob_acc = 0
          prob_after_previous = 0
        } 
        per = per + 1
      }
      schedule_data$npv_at_risk = npv_at_risk
      return(schedule_data)
    }
    
    no_of_threads = 8
    cl = makeCluster(no_of_threads)
    
    clusterEvalQ(cl, library(data.table))
    
    clusterExport(cl, c("multiyear_prob",
                        "n_per"),
                  envir = environment())
    
    system.time({
      results_substand_lvl = parLapply(cl, results_substand_lvl_split, fun = calculate_probabilities) %>%
        rbindlist()
    })
    stopCluster(cl)
    
    results_substand_lvl[,break_p_period:=1-(1-multiyear_prob(break_p_before,5))*(1-multiyear_prob(break_p_after,5))]
    results_substand_lvl[,csl := rowMeans(.SD, na.rm = TRUE), .SDcols = c("csl_start", "csl_after")]  
    results_substand_lvl[,break_p:=rowMeans(.SD, na.rm = TRUE), .SDcols = c("break_p_before", "break_p_after")]  
    
    results_substand_lvl[,value_at_risk:=net*prob_acc]
    fwrite(results_substand_lvl,paste0(path_proj,"\\Results\\results_substand_lvl.csv"))
  }
  
  if (process_GAYA_SCHEDULES_SUBSTAND_LEVEL_ADD_NEW){
    if (!exists("results_substand_lvl")) results_substand_lvl = fread(paste0(path_proj,"\\Results\\results_substand_lvl.csv"))
    results_standard_new = fread("C:/Users/victo/OneDrive/Desktop/Marcsman/Results/PT_marcsman_standard09.csv") 
    results_fritzoe_new = fread("C:/Users/victo/OneDrive/Desktop/Marcsman/Results/PT_marcsman_fritzoe09.csv")
    results_substand_lvl_new = rbind(results_standard_new,results_fritzoe_new) %>% .[period %in% 1:n_per]
    new_cols = setdiff(names(results_substand_lvl_new), names(results_substand_lvl))
    results_substand_lvl = cbind(results_substand_lvl, results_substand_lvl_new[,..new_cols])
    fwrite(results_substand_lvl, paste0(path_proj,"\\Results\\results_substand_lvl.csv"))
  }
  
  if (process_GAYA_SCHEDULES_STAND_LEVEL){
    if (!exists("results_substand_lvl")) results_substand_lvl = fread(paste0(path_proj,"\\Results\\results_substand_lvl.csv"))
    results_stand_lvl = 
      results_substand_lvl[, .(
        n = weighted.mean(ifelse(n > 0, n, NA), area, na.rm = TRUE),
        reg_n = weighted.mean(ifelse(reg_n > 0, n, NA), area, na.rm = TRUE),
        n_tot = weighted.mean(ifelse(n + reg_n > 0, n + reg_n, NA), area, na.rm = TRUE),
        h = weighted.mean(ifelse(h > 0, h, NA), area, na.rm = TRUE),
        d = weighted.mean(ifelse(d > 0, d, NA), area, na.rm = TRUE),
        dh_ratio = weighted.mean(d/h, area, na.rm = TRUE),
        asl = pmax(max_sl, na.rm = TRUE),
        break_p_before = weighted.mean(break_p_before, area, na.rm = TRUE),
        break_p_after = weighted.mean(break_p_after, area, na.rm = TRUE),
        break_p = weighted.mean(break_p, area, na.rm = TRUE),
        break_p_period = weighted.mean(break_p_period, area, na.rm = TRUE),
        npv_at_risk = weighted.mean(npv_at_risk, area, na.rm = TRUE),
        value_at_risk_ha = weighted.mean(value_at_risk, area, na.rm = TRUE),
        value_at_risk_tot = sum(value_at_risk*area, na.rm = TRUE),
        csl_start = weighted.mean(csl_start, area, na.rm = TRUE),
        csl_after = weighted.mean(csl_after, area, na.rm = TRUE),
        csl = weighted.mean(csl, area, na.rm = TRUE),
        net = sum(net*area, na.rm = TRUE),
        cost = sum(cost*area, na.rm = TRUE),
        income = sum(income*area, na.rm = TRUE),
        price_m3 = weighted.mean(price_m3[price_m3 != 0 & !is.na(price_m3)], area[price_m3 != 0 & !is.na(price_m3)]),
        costF_m3 = weighted.mean(costF_m3[costF_m3 != 0 & !is.na(costF_m3)], area[costF_m3 != 0 & !is.na(costF_m3)]),
        costC_m3 = weighted.mean(costC_m3[costC_m3 != 0 & !is.na(costC_m3)], area[costC_m3 != 0 & !is.na(costC_m3)]),
        costT_m3 = weighted.mean(costT_m3[costT_m3 != 0 & !is.na(costT_m3)], area[costT_m3 != 0 & !is.na(costT_m3)]),
        costR = sum(costR*area, na.rm = TRUE),
        vol = sum(vol*area, na.rm = TRUE),
        vol_trt = sum(vol_trt*area, na.rm = TRUE),
        vol_trt_saw = sum(vol_trt_saw*area, na.rm = TRUE),
        vol_trt_pul = sum(vol_trt_pul*area, na.rm = TRUE),
        vol_trt_S = sum(vol_trt_S*area, na.rm = TRUE),
        vol_trt_P = sum(vol_trt_P*area, na.rm = TRUE),
        vol_trt_B = sum(vol_trt_B*area, na.rm = TRUE),
        area = sum(area, na.rm = TRUE)
      ), by = .(stand, period, management)]
    
    # add species
    stand_data = fread(paste0(path_proj,"\\Data\\marcsman_stands.csv")) %>% .[,.(stand=standID,species=mainSp)]
    results_stand_lvl = left_join(results_stand_lvl,stand_data)
    fwrite(results_stand_lvl,paste0(path_proj,"\\Results\\results_stand_lvl.csv"))
  }
  
  if (process_GAYA_SCHEDULES_PROPERTY_LEVEL){
    if (!exists("results_substand_lvl")) results_substand_lvl = fread(paste0(path_proj,"\\Results\\results_substand_lvl.csv"))
    results_property_lvl = 
      results_substand_lvl[, .(
        n = weighted.mean(ifelse(n > 0, n, NA), area, na.rm = TRUE),
        reg_n = weighted.mean(ifelse(reg_n > 0, n, NA), area, na.rm = TRUE),
        n_tot = weighted.mean(ifelse(n + reg_n > 0, n + reg_n, NA), area, na.rm = TRUE),
        h = weighted.mean(ifelse(h > 0, h, NA), area, na.rm = TRUE),
        d = weighted.mean(ifelse(d > 0, d, NA), area, na.rm = TRUE),
        dh_ratio = weighted.mean(d/h, area, na.rm = TRUE),
        
        csl_m2 = weighted.mean(csl, area, na.rm = TRUE),
        break_p_proc = 100*weighted.mean(break_p, area, na.rm = TRUE),
        value_at_risk_ha_yr = weighted.mean(value_at_risk, area, na.rm = TRUE)/10,
        
        vol_harvest_tot = sum(vol_trt * area, na.rm = TRUE),
        vol_harvest_thinning_tot = sum(vol_trt[treatment=="Thinning"] * area[treatment=="Thinning"], na.rm = TRUE),
        vol_harvest_harvest_tot = sum(vol_trt[treatment=="Harvest"] * area[treatment=="Harvest"], na.rm = TRUE),
        vol_harvest_saw_tot = sum(vol_trt_saw * area, na.rm = TRUE),
        vol_harvest_pul_tot = sum(vol_trt_pul * area, na.rm = TRUE),
        vol_harvest_S_tot = sum(vol_trt_S * area, na.rm = TRUE),
        vol_harvest_P_tot = sum(vol_trt_P * area, na.rm = TRUE),
        vol_harvest_B_tot = sum(vol_trt_B * area, na.rm = TRUE),
        vol_standing_tot = sum(vol * area, na.rm = TRUE),
        
        value_harv_tot = sum(price_m3 * vol_trt * V.a1 * area, na.rm = TRUE),
        cost_logging_tot = sum(costC_m3 * vol_trt * V.a1 * area, na.rm = TRUE),
        cost_forward_tot = sum(costF_m3 * vol_trt * V.a1 * area, na.rm = TRUE),
        cost_harvest_tot = sum(costT_m3 * vol_trt * V.a1 * area, na.rm = TRUE),
        income_harvest_tot = sum(price_m3 * vol_trt * V.a1 * area, na.rm = TRUE) -
          sum(costT_m3 * vol_trt * V.a1 * area, na.rm = TRUE),
        cost_planting_tot = sum(ifelse(reg != 0, costR * area, 0), na.rm = TRUE),
        cost_tending_tot = sum(ifelse(reg == 0, costR * area, 0), na.rm = TRUE),  
        net_tot = sum(net*area, na.rm = TRUE),
        
        value_harvest_m3 = weighted.mean(price_m3[price_m3 != 0], area[price_m3 != 0], na.rm = TRUE),
        cost_logging_m3 = weighted.mean(costC_m3[costC_m3 != 0], area[costC_m3 != 0], na.rm = TRUE),
        cost_forward_m3 = weighted.mean(costF_m3[costF_m3 != 0], area[costF_m3 != 0], na.rm = TRUE),
        cost_harvest_m3 = weighted.mean(costT_m3[costT_m3 != 0], area[costT_m3 != 0], na.rm = TRUE),
        income_harvest_m3 = weighted.mean(price_m3[price_m3 != 0], area[price_m3 != 0], na.rm = TRUE) - 
          weighted.mean(costT_m3[costT_m3 != 0], area[costT_m3 != 0], na.rm = TRUE),
        
        cost_planting_ha = weighted.mean(costR[costR != 0 & reg != 0], area[costR != 0 & reg != 0], na.rm = TRUE),
        cost_tending_ha = weighted.mean(costR[costR != 0 & reg == 0], area[costR != 0 & reg == 0], na.rm = TRUE),
        
        area_tending = sum(area[treatment=="Tending"], na.rm = TRUE),
        area_thinning = sum(area[treatment=="Thinning"], na.rm = TRUE),
        area_harvest = sum(area[treatment=="Harvest"], na.rm = TRUE)
      ), 
      by = .(period, management)]
    
    fwrite(results_property_lvl,paste0(path_proj,"\\Results\\results_property_lvl.csv"))
  }
}

if (figure_STUDY_AREA){
  vot_bbox = getbb("Vestfold og Telemark, Norway")
  
  vot_boundaries =  opq(bbox = vot_bbox) %>% 
    add_osm_feature(key = 'boundary', value = 'administrative') %>%
    osmdata_sf() 
  
  vot_boundaries = subset(vot_boundaries$osm_multipolygons, name=="Vestfold og Telemark")
  
  coastline_data = opq(bbox = vot_bbox) %>%
    add_osm_feature(key = 'natural', value = 'coastline') %>%
    osmdata_sf()
  
  coastline_sf = coastline_data$osm_lines  
  coastline_sf = st_transform(coastline_sf, st_crs(vot_boundaries))
  
  blade = coastline_sf %>% st_union %>% st_line_merge
  
  # get land borders of individual municipalities
  
  # SANDEFJORD
  Sandefjord = getbb('Sandefjord, Norway',format_out = "sf_polygon")[[2]]
  Sandefjord = st_split(st_geometry(Sandefjord), st_geometry(blade))
  Sandefjord_land = st_cast(Sandefjord[[1]], 'POLYGON') %>% st_geometry %>% st_sf
  st_crs(Sandefjord_land) = st_crs(Sandefjord)
  Sandefjord_land$name = "Sandefjord"
  
  # LARVIK
  Larvik = getbb('Larvik, Norway',format_out = "sf_polygon")
  Larvik = st_split(st_geometry(Larvik), st_geometry(blade))
  Larvik_land = st_cast(Larvik[[1]], 'POLYGON') %>% st_geometry %>% st_sf
  st_crs(Larvik_land) = st_crs(Larvik)
  Larvik_land$name = "Larvik"
  
  # PORSGRUNN
  Porsgrunn = getbb('Porsgrunn, Norway',format_out = "sf_polygon")
  Porsgrunn_split = st_split(st_geometry(Porsgrunn), st_geometry(blade))
  Porsgrunn_water = st_cast(Porsgrunn_split[[1]], 'POLYGON') %>% st_geometry %>% st_sf
  st_crs(Porsgrunn_water) = st_crs(Porsgrunn)
  Porsgrunn_land = st_difference(Porsgrunn, Porsgrunn_water) %>% st_geometry %>% st_sf
  st_crs(Porsgrunn_land) = st_crs(Porsgrunn)
  Porsgrunn_land$name = "Porsgrunn"
  
  # SILJAN
  Siljan = getbb('Siljan, Norway',format_out = "sf_polygon")
  Siljan = st_split(st_geometry(Siljan), st_geometry(blade))
  Siljan_land = st_cast(Siljan[[1]], 'POLYGON') %>% st_geometry %>% st_sf
  st_crs(Siljan_land) = st_crs(Siljan)
  Siljan_land$name = "Siljan"
  
  # KONGSBERG
  Kongsberg = getbb('Kongsberg, Norway',format_out = "sf_polygon")
  Kongsberg = st_split(st_geometry(Kongsberg), st_geometry(blade))
  Kongsberg_land = st_cast(Kongsberg[[1]], 'POLYGON') %>% st_geometry %>% st_sf
  st_crs(Kongsberg_land) = st_crs(Kongsberg)
  Kongsberg_land$name = "Kongsberg"
  
  all.kommune = rbind(Sandefjord_land,
                      Larvik_land,
                      Porsgrunn_land[1,],
                      Siljan_land,
                      Kongsberg_land)
  
  bbox = st_bbox(all.kommune)
  
  bbox["xmin"] = bbox["xmin"] - 0.1
  bbox["xmax"] = bbox["xmax"] + 0.1
  bbox["ymin"] = bbox["ymin"] - 0.1
  bbox["ymax"] = bbox["ymax"] + 0.1
  
  bbox_polygon = st_as_sfc(bbox)
  split_polygons = st_split(bbox_polygon, blade)
  split_polygons_sf = st_collection_extract(split_polygons, "POLYGON")
  
  water_sf = split_polygons_sf[2]
  land_sf = split_polygons_sf[-2]
  
  # start with the municipality centroids
  all.kommune.centroids = st_centroid(all.kommune)
  labels_dt = data.table(name = all.kommune$name,
                         x = st_coordinates(all.kommune.centroids)[,1],
                         y = st_coordinates(all.kommune.centroids)[,2])
  unit.y = 0.01
  unit.x = 0.01
  
  # manually adjust some of the label positions
  labels_dt[name=="Sandefjord"]$y = labels_dt[name=="Sandefjord"]$y+2*unit.y
  labels_dt[name=="Siljan"]$x = labels_dt[name=="Siljan"]$x-unit.x
  labels_dt[name=="Porsgrunn"]$x = labels_dt[name=="Porsgrunn"]$x-2*unit.x
  labels_dt[name=="Porsgrunn"]$y = labels_dt[name=="Porsgrunn"]$y+unit.y
  labels_dt[name=="Larvik"]$x = labels_dt[name=="Larvik"]$x+4*unit.x
  labels_dt[name=="Larvik"]$y = labels_dt[name=="Larvik"]$y-2*unit.y
  
  # color scheme
  land_fill = "gray80"
  water_fill = "#ADD8E6"
  study_area_fill = "#98FB98"
  plots_fill = "gray80"
  plots_color = "black"
  
  bbox_study_area = st_bbox(all.kommune)
  
  forest_stands = st_read(paste0(path_proj,"\\Data\\Data_Fritsoe_Tron\\Bestand20t.shp"))
  
  
  figure_study_area = 
    ggplot() +
    geom_sf(data = land_sf, fill = land_fill, color = NA) +
    geom_sf(data = forest_stands, color = NA, fill = study_area_fill) +
    geom_sf(data = all.kommune, fill = NA, lwd = 0.5, color="gray40") +
    geom_sf(data = water_sf, fill = water_fill, lwd = 0.3, color = "black") +
    geom_text_repel(data = labels_dt, aes(x = x, y = y, label = name), size = 3, color = "black",force=0) +
    coord_sf(xlim = c(9.5, 10.25), ylim = c(58.9, 59.6)) +
    labs(title="Study area", x ="Longitude", y ="Latitude") +
    theme_bw()
  
  # NORTHERN EUROPE
  northern_europe = c("Norway","Sweden","Denmark","Finland","Estonia","Latvia","Lithuania")
  northern_europe_map = ne_countries(country = northern_europe,
                                     scale = 50, 
                                     returnclass = 'sf')
  
  bbox_northern_europe = st_bbox(c(xmin = 4, ymin = 20, xmax = 50, ymax = 70),st_crs(northern_europe_map))
  
  forest_stands = st_transform(forest_stands, st_crs(northern_europe_map))
  
  northern_europe_map_cropped = st_crop(northern_europe_map, bbox_northern_europe)
  northern_europe_labels = data.table(name = northern_europe,x = c(9.5,15,9,27,25.5,25.8,24.2),y = c(62,63,56.5,63,59,56.85,55.5))
  
  bbox = st_bbox(forest_stands)
  bbox["xmin"] = bbox["xmin"] - 0.5
  bbox["xmax"] = bbox["xmax"] + 0.5
  bbox["ymin"] = bbox["ymin"] - 0.5
  bbox["ymax"] = bbox["ymax"] + 0.5
  
  bbox_polygon = st_as_sfc(bbox)
  
  figure_northern_europe = 
    ggplot() + 
    geom_sf(data = northern_europe_map_cropped, fill = land_fill) + 
    geom_sf(data = forest_stands, fill = "black",color=NA) +
    geom_text(data = northern_europe_labels, aes(x = x,y = y,label = name), size = 4) + 
    geom_sf(data = bbox_polygon, fill=NA, lwd = 0.5, color= "red") + 
    labs(title="Northern Europe", x ="Longitude", y = "Latitude") +
    theme_bw()
  
  figure = grid.arrange(figure_northern_europe,
                        figure_study_area,
                        nrow = 1)
  ggsave(
    filename = paste0(path_proj,"\\Figures\\figure_STUDY_AREA.jpg"), 
    plot = figure, 
    width = 1.2*twocol, height = 1.2*maxHeight*.5, units = "mm",dpi = 1000)
  
} # FIGURE 1

if (figure_CYCLES){
  if (!exists("results_substand_lvl")) results_substand_lvl = fread(paste0(path_proj,"\\Results\\results_substand_lvl.csv"))
  
  results_substand_lvl_cycles = results_substand_lvl[,.(stand,schedule,management,period,program_track,treatment,area)]
  split_list = split(results_substand_lvl_cycles, by = c("stand", "schedule", "management"))
  
  apply_cycle_logic = function(dt) {
    dt[, cycle := 1]
    first_program_track = 0
    first_program_track = dt[program_track == 1, min(period, na.rm = TRUE)]
    
    if (first_program_track %in% 1:n_per) {
      dt[period >= first_program_track, cycle := 2]
      first_harvest = 0 
      first_harvest = dt[period > first_program_track & treatment == "Harvest", min(period, na.rm = TRUE)]
      if (first_harvest %in% 1:n_per) dt[period > first_harvest, cycle := 3]
    }
    return(dt)
  }
  
  split_list = lapply(split_list, apply_cycle_logic)
  results_substand_lvl_cycles = rbindlist(split_list)
  
  percentage_data = results_substand_lvl_cycles[, .(area = sum(area)), by = .(period, cycle, management)] %>%
    merge(results_substand_lvl_cycles[, .(area = sum(area)), by = .(period, management)], by = c("period", "management"), suffixes = c("_cycle", "_total")) %>%
    .[, percentage := (area_cycle / area_total) * 100]
  
  cycle_labels = c("1" = "Current rotation", "2" = "Second rotation", "3" = "Third rotation")
  manag_levels = c("standard", "fritzoe")
  percentage_data$management = factor(percentage_data$management, levels = manag_levels)
  labels = c( 
    "standard" = "Standard management",
    "fritzoe" = "Reduced density management")
  
  figure =
    ggplot(percentage_data[period<=period_last], aes(x = factor(period), y = percentage, fill = factor(cycle))) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~management, labeller = labeller(management = labels), dir='v') + 
    labs(x = "Period", y = "Area (%)", fill = "") +
    theme_minimal() +
    scale_fill_grey(labels = cycle_labels) + 
    theme(legend.position = "bottom")
  
  ggsave(
    filename = paste0(path_proj,"\\Figures\\figure_ROTATIONS.jpg"), 
    plot = figure, 
    width = 1.2*onecol, height = 1.2*maxHeight*.3, units = "mm",dpi = 1000)
} # FIGURE 3

if (figure_INITIAL_STATE | figure_DIFF){
  if (!exists("results_stand_lvl")) results_stand_lvl = fread(paste0(path_proj,"\\Results\\results_stand_lvl.csv"))
  
  pmean = function(x) {
    x = x[is.finite(x)]
    if (length(x) == 0) 
      return(as.numeric(NA))
    else return(mean(x))
  }
  
  results_stand_lvl_x =
    results_stand_lvl[,.(
      area = first(area[period == 1]),
      asl = pmean(asl),
      csl_start = first(csl_start[period == 1]),
      csl_later = mean(csl[period >= period_later & period <= period_last], na.rm = T),
      break_p_proc_start = first(break_p[period == 1])*100,
      break_p_proc_later = mean(break_p[period >= period_later & period <= period_last], na.rm = T)*100, 
      value_at_risk_ha_yr = nok_to_eur*mean(value_at_risk_ha[period >= period_later & period <= period_last_value], na.rm = T)/10
    ), 
    by = .(stand,management)]
  
  results_stand_lvl_later = results_stand_lvl_x[,.(
    area,
    asl, 
    break_p_proc = break_p_proc_later,
    csl_m2 = csl_later,
    value_at_risk_ha_yr,
    stand,
    management)]
  
  results_stand_lvl_start = results_stand_lvl_x[management == "standard",.(
    area,
    asl,
    break_p_proc = break_p_proc_start,
    csl_m2 = csl_start,
    value_at_risk_ha_yr,
    stand,
    management = "initial")]
  
  results_stand_lvl_selected = rbind(results_stand_lvl_start, results_stand_lvl_later)  
  
  results_property_lvl_selected = 
    results_stand_lvl_selected[,.(
      asl = weighted.mean(asl, area, na.rm = T),
      break_p_proc = weighted.mean(break_p_proc, area, na.rm = T),
      csl_m2 = weighted.mean(csl_m2, area, na.rm = T),
      value_at_risk_ha_yr = weighted.mean(value_at_risk_ha_yr[value_at_risk_ha_yr >= 0], area[value_at_risk_ha_yr >= 0], na.rm = T)
    ),
    by = management]
  
  results_stand_lvl_selected_long = 
    melt(results_stand_lvl_selected,
         id.vars = c("stand", "management"), 
         variable.name = "variable", 
         value.name = "val")
  
  fwrite(results_stand_lvl_selected_long, paste0(path_proj,"\\Results\\results_stand_lvl_selected_long.csv"))
  
  forest_stands = 
    st_read(paste0(path_proj,"\\Data\\Data_Fritsoe_Tron\\Bestand20t.shp")) %>%
    mutate(stand = (BestNr*100+DelNr)*10+Teig) %>%
    left_join(results_stand_lvl_selected_long, by = "stand") 
}

if (figure_INITIAL_STATE){
  
  plot_map_hist = function (var_name,scale_direction,plot_title,x_label,mean_val,histo_bins,geom_text_x_add,geom_text_y,precision,x_min,x_max,y_min,y_max,hide_axis.title.y,hide_axis.text.y){
    
    filtered_data = filter(forest_stands, management == "initial" & variable == var_name)
    quantiles = as.numeric(quantile(filtered_data$val, probs = seq(0, 1, by = 0.001), na.rm = TRUE))
    normalized_quantiles = (quantiles - min(quantiles)) / (max(quantiles) - min(quantiles))
    
    figure_map = 
      ggplot(filtered_data) +
      geom_sf(color = NA, aes(fill = val)) +
      theme_bw() +
      theme(legend.position = "none") + 
      labs(title = plot_title, x = "Longitude", y = "Latitude") +
      scale_fill_viridis_c(
        values = normalized_quantiles,
        option = "H",
        direction = scale_direction
      )
    
    figure_hist = 
      ggplot(filtered_data, aes(x=val)) +
      geom_histogram(aes(fill = ..x..), bins = histo_bins) + 
      geom_vline(xintercept = mean_val, linetype = "dashed") +
      geom_text(x = mean_val + geom_text_x_add, y = geom_text_y ,label = sprintf(precision, mean_val), size = 4) +
      theme_bw() +
      theme(legend.position = "none",
            strip.background = element_blank(),  
            strip.text = element_blank()
      ) + 
      labs(x = x_label, y = "Count") +
      scale_fill_viridis_c(
        values = normalized_quantiles,
        option = "H",
        direction = scale_direction
      ) + 
      xlim(x_min, x_max)+
      ylim(y_min, y_max)
    
    if (hide_axis.title.y) {
      figure_map = figure_map + theme(axis.title.y = element_blank())
      figure_hist = figure_hist + theme(axis.title.y = element_blank())
    }
    
    if (hide_axis.text.y) {
      figure_map = figure_map + theme(axis.text.y = element_blank())
      figure_hist = figure_hist + theme(axis.text.y = element_blank())
    }
    
    
    h_ratio = .4
    g_m = ggplotGrob(figure_map)
    g_h = ggplotGrob(figure_hist)
    
    g_h$heights[9] = g_h$heights[9]*h_ratio
    
    figure = rbind(g_m, g_h, size="first")
    figure
  }
  
  init_csl_plot = 
    plot_map_hist(
      var_name = "csl_m2",
      scale_direction = -1,
      plot_title = "Critical snow load",
      x_label = bquote("Critical snow load ("*Kg~m^-2*")"),
      mean_val = results_property_lvl_selected[management == "initial", csl_m2],
      histo_bins = 60,
      geom_text_x_add = 20,
      geom_text_y = 1000,
      precision = "%.2f",
      x_min = 0, x_max = 200,
      y_min = 0, y_max = 1250,
      hide_axis.title.y = F, hide_axis.text.y = F
    )
  
  init_asl_plot = 
    plot_map_hist(
      var_name = "asl",
      scale_direction = 1,
      plot_title = "Expected maximum snow load",
      x_label = bquote("Expected snow load ("*Kg~m^-2*")"),
      mean_val = results_property_lvl_selected[management == "initial", asl],
      histo_bins = 60,
      geom_text_x_add = 10,
      geom_text_y = 1000,
      precision = "%.2f",
      x_min = 40, x_max = 120,
      y_min = 0, y_max = 1250,
      hide_axis.title.y = T, hide_axis.text.y = T
    )
  
  init_break_p_proc_plot = 
    plot_map_hist(
      var_name = "break_p_proc",
      scale_direction = 1,
      plot_title = "Snow damage probability",
      x_label = bquote("Snow damage probability ("*`%`*~yr^-1*")"),
      mean_val = results_property_lvl_selected[management == "initial", break_p_proc],
      histo_bins = 60,
      geom_text_x_add = 0.2,
      geom_text_y = 1000,
      precision = "%.2f",
      x_min = 0, x_max = 2,
      y_min = 0, y_max = 1250,
      hide_axis.title.y = T, hide_axis.text.y = T
    )
  
  figure = cbind(init_csl_plot, init_asl_plot, init_break_p_proc_plot, size="first")
  
  ggsave(
    filename = paste0(path_proj,"\\Figures\\figure_INITIAL_STATE.jpg"), 
    plot = figure, 
    width = 1.2*twocol, height = 1.2*maxHeight*.6, units = "mm",dpi = 1000)
} # FIGURE 2

if (figure_DIFF){
  if (!exists("results_stand_lvl")) results_stand_lvl = fread(paste0(path_proj,"\\Results\\results_stand_lvl.csv"))
  
  results_property_lvl_selected_per = 
    results_stand_lvl[, .(
      n_tot = weighted.mean(n_tot, area, na.rm = TRUE),
      dh_ratio = weighted.mean(dh_ratio, area, na.rm = TRUE),
      csl_m2 = weighted.mean(csl, area, na.rm = TRUE),
      break_p_proc = 100*weighted.mean(break_p, area, na.rm = TRUE),
      value_at_risk_ha_yr = nok_to_eur*weighted.mean(value_at_risk_ha, area, na.rm = TRUE)/10
    ), 
    by = .(period, management)]
  
  results_property_lvl_selected_per$management = factor(results_property_lvl_selected_per$management, 
                                                        levels = c("standard", "fritzoe"))
  
  results_property_lvl_selected_per_long = 
    melt(results_property_lvl_selected_per[period %in% 1:period_last,
                                           .(period, management, csl_m2, break_p_proc, n_tot, dh_ratio)], 
         id.vars = c("period", "management"), 
         variable.name = "variable", 
         value.name = "val")
  
  if (figure_DIFF_CSL_SDP){
    
    red_to_blue_11 =  c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#f0f5f7','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695','#313695')
    
    pseudo_log = function(x){
      if (abs(x)<exp(1)) 
        return(x/exp(1))
      else 
        return(log(abs(x))*(x/abs(x)))
    }
    
    range.val.from.pos = function(pos, min, max){
      range.width = max-min
      val = pos*(max-min)+min
      return(val)
    }
    
    range.pos.from.val = function(val, min, max){
      pos = (val-min)/(max-min)
      return(pos)
    }
    
    get.color.scale = function(val.min, val.max, colors, nbins){
      val.min = val.min 
      val.max = val.max 
      colors.n = length(colors)
      colors.mid.index = round(range.val.from.pos(.5,1,colors.n))
      zero.pos.range = range.pos.from.val(0,val.min,val.max)
      
      colors.min.index = 1
      colors.max.index = colors.n
      
      if (zero.pos.range<0)
        colors.min.index = round(range.val.from.pos(range.pos.from.val(val.min,0,val.max),colors.mid.index,colors.n))
      if (zero.pos.range>1)
        colors.max.index = round(range.val.from.pos(range.pos.from.val(val.max,val.min,0),1,colors.mid.index))
      
      colors.index.subset = seq(from=colors.min.index, to=colors.max.index, by=1)
      colors.subset = colors[colors.index.subset]
      
      if (zero.pos.range>0&zero.pos.range<=1)
        sc.values = c(seq(from = 0, to = zero.pos.range, length.out = colors.mid.index),seq(from = zero.pos.range, to = 1, length.out = colors.mid.index)[2:colors.mid.index])
      
      sc = scale_colour_gradientn(
        colors = colors.subset,
        space = "Lab",
        na.value = "grey50",
        guide = NULL,
        aesthetics = "fill",
        #right=F,
        trans = "pseudo_log",
        name = bquote(""*Mg~ha^-1~yr^-1),
        n.breaks = nbins
        #values = sc.values
      )
      
      original_seq = seq(val.min,val.max,length=nbins)
      trans_seq = apply(as.array(original_seq),FUN=pseudo_log,MARGIN=1)
      #rescaled_seq = apply(as.array(trans_seq), FUN=range.pos.from.val, min = min(trans_seq),max = max(trans_seq), MARGIN=1)
      
      if (zero.pos.range>=0&zero.pos.range<=1){
        zero.index.nbins = round(range.val.from.pos(zero.pos.range,1,nbins))
        rescaled_seq_1 = apply(as.array(trans_seq[1:zero.index.nbins]),FUN=sc$rescaler,to=c(0,0.5),from=c(min(trans_seq),0),MARGIN=1)
        rescaled_seq_2 = apply(as.array(trans_seq[(zero.index.nbins+1):nbins]),FUN=sc$rescaler,to=c(0.5,1),from=c(0,max(trans_seq)),MARGIN=1)
        rescaled_seq = c(rescaled_seq_1,rescaled_seq_2)
      } else {
        rescaled_seq = apply(as.array(trans_seq),FUN=sc$rescaler,to=c(0,1),from=c(min(trans_seq),max(trans_seq)),MARGIN=1)
      }
      
      cols = apply(as.array(rescaled_seq),FUN=sc$palette,MARGIN=1)
      
      sc = scale_colour_gradientn(
        colors = cols,
        space = "Lab",
        na.value = "grey50",
        guide = NULL,
        aesthetics = "fill",
        #right=F,
        #trans = "pseudo_log",
        name = bquote(""*Mg~ha^-1~yr^-1),
        n.breaks = nbins
      )
      
      return (list(cols = cols, sc = sc))
    }
    
    val_by_manag = results_property_lvl_selected[management != "initial"]
    
    ### CSL DIFF ###
    filtered_data = filter(forest_stands, management %in% c("fritzoe", "standard") & variable == "csl_m2")
    
    spread_data = filtered_data %>%
      pivot_wider(names_from = management, values_from = val) 
    
    spread_data = spread_data %>%
      mutate(val_diff = fritzoe - standard) 
    
    multiplier = 5
    nbins = 100
    cs = get.color.scale(min(spread_data$val_diff)*multiplier, max(spread_data$val_diff)*multiplier, red_to_blue_11, nbins)
    
    figure_map_csl = 
      ggplot(spread_data) +
      geom_sf(color = NA, aes(fill = val_diff)) +
      cs$sc + 
      theme_bw() +
      theme(legend.position = "none") + 
      labs(title = "Difference in CSL", x = "Longitude", y = "Latitude", fill = "Difference in CSL")
    
    figure_hist_csl =
      ggplot(spread_data, aes(x = val_diff)) +
      geom_histogram(bins = nbins, fill = cs$cols) +
      theme_bw() +
      theme(legend.position = "none",
            strip.background = element_blank(),  
            strip.text = element_blank()
      ) + 
      labs(x = bquote("Difference in critical snow load ("*Kg~m^-2*")"), y = "Count") + 
      ylim(0, 3500)
    
    figure_timeseries_csl =
      ggplot(results_property_lvl_selected_per_long[variable == "csl_m2"], 
             aes(x = period, 
                 y = val, 
                 color = management)) +
      annotate("rect", xmin = 5, xmax = 10, ymin = -Inf, ymax = Inf, fill = "gray90") +
      geom_line() +
      geom_segment(data = val_by_manag, 
                   aes(x = 5, xend = 10, y = csl_m2, yend = csl_m2, color = management),
                   linetype = "dashed",
                   show.legend = FALSE) + 
      geom_text(data = val_by_manag, 
                aes(x = 5.1, y = csl_m2, label = round(csl_m2, 2), color = management), 
                hjust = 0, 
                vjust = -0.5,
                show.legend = FALSE,
                size = 3) +
      scale_y_continuous(labels = label_number())+
      labs(x = "Period",
           y = bquote("Critical snow load ("*Kg~m^-2*")")) +
      theme_bw()+
      scale_x_continuous(breaks = 1:period_last_value) +
      theme(legend.position = "none",
            legend.title = element_blank(),
            axis.title.x = element_blank()) +
      scale_color_manual(
        values = c("standard" = "black", "fritzoe" = "#00BFFF"),  
        labels = c("standard" = "Standard management", "fritzoe" = "Reduced density management")
      ) 
    
    figure =
      (figure_map_csl + 
         (figure_hist_csl / figure_timeseries_csl)) + 
      plot_layout(widths = c(1, 1))
    
    ggsave(
      filename = paste0(path_proj,"\\Figures\\figure_CSL_DIFF.jpg"), 
      plot = figure, 
      width = 1.2*twocol, height = 1.2*maxHeight*.5, units = "mm",dpi = 1000)
    
    
    ### SDP DIFF ###
    filtered_data = filter(forest_stands, management %in% c("fritzoe", "standard") & 
                             variable == "break_p_proc")
    
    spread_data = filtered_data %>%
      pivot_wider(names_from = management, values_from = val) 
    
    spread_data = spread_data %>%
      mutate(val_diff = fritzoe - standard) 
    
    multiplier = 200
    nbins = 100
    cs = get.color.scale(min(spread_data$val_diff)*multiplier, max(spread_data$val_diff)*multiplier, rev(red_to_blue_11), nbins)
    
    figure_map_sdp = 
      ggplot(spread_data) +
      geom_sf(color = NA, aes(fill = val_diff)) +
      cs$sc + 
      theme_bw() +
      theme(legend.position = "none",
            axis.title.y = element_blank(),
            axis.text.y = element_blank()) + 
      
      labs(title = "Difference in SDP", x = "Longitude", y = "Latitude", fill = "Difference in SDP")
    
    figure_hist_sdp =
      ggplot(spread_data, aes(x=val_diff)) +
      geom_histogram(bins = nbins, fill = cs$cols) +
      theme_bw() +
      theme(legend.position = "none",
            strip.background = element_blank(),  
            strip.text = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank()
      ) + 
      labs(x = bquote("Difference in Snow damage probability ("*`%`*~yr^-1*")"), y = "Count") +
      ylim(0, 3500)
    
    figure_timeseries_sdp =
      ggplot(results_property_lvl_selected_per_long[variable == "break_p_proc"], 
             aes(x = period, 
                 y = val, 
                 color = management)) +
      annotate("rect", xmin = 5, xmax = 10, ymin = -Inf, ymax = Inf, fill = "gray90") +
      geom_line() +
      geom_segment(data = val_by_manag, 
                   aes(x = 5, xend = 10, y = break_p_proc, yend = break_p_proc, color = management),
                   linetype = "dashed",
                   show.legend = FALSE) + 
      geom_text(data = val_by_manag, 
                aes(x = 5.1, y = break_p_proc, label = round(break_p_proc, 3), color = management), 
                hjust = 0, 
                vjust = -0.5,
                show.legend = FALSE,
                size = 3) +    
      scale_y_continuous(labels = label_number()) +
      labs(x = "Period",
           y = bquote("Snow damage probability ("*`%`*~yr^-1*")")) +
      theme_bw() +
      scale_x_continuous(breaks = 1:period_last_value) +
      theme(legend.position = "bottom",
            legend.title = element_blank()) +
      scale_color_manual(
        values = c("standard" = "black", "fritzoe" = "#00BFFF"),  
        labels = c("standard" = "Standard management", "fritzoe" = "Reduced density management")
      ) 
    
    figure =
      (figure_map_sdp + 
         (figure_hist_sdp / figure_timeseries_sdp)) + 
      plot_layout(widths = c(1, 1))
    
    ggsave(
      filename = paste0(path_proj,"\\Figures\\figure_SDP_DIFF.jpg"), 
      plot = figure, 
      width = 1.2*twocol, height = 1.2*maxHeight*.5, units = "mm",dpi = 1000)
    
    ### NEW ARRANGEMENT ###
    
    figure_timeseries =
      (figure_timeseries_csl / figure_timeseries_sdp) + 
      plot_layout(widths = c(1, 1))
    
    h_ratio = .4
    
    g_m = ggplotGrob(figure_map_csl)
    g_h = ggplotGrob(figure_hist_csl)
    g_h$heights[9] = g_h$heights[9]*h_ratio
    figure_csl = rbind(g_m, g_h, size="first")
    
    g_m = ggplotGrob(figure_map_sdp)
    g_h = ggplotGrob(figure_hist_sdp)
    g_h$heights[9] = g_h$heights[9]*h_ratio
    figure_sdp = rbind(g_m, g_h, size="first")
    
    figure_maps = cbind(figure_csl, figure_sdp, size = "first")
    
    ggsave(
      filename = paste0(path_proj,"\\Figures\\figure_TIMESERIES_DIFF.jpg"), 
      plot = figure_timeseries, 
      width = 0.7*twocol, height = maxHeight*.6, units = "mm",dpi = 1000)
    
    ggsave(
      filename = paste0(path_proj,"\\Figures\\figure_MAPS_DIFF.jpg"), 
      plot = figure_maps, 
      width = 1.2*twocol, height = 1.2*maxHeight*0.8, units = "mm",dpi = 1000)
  }
} # GRAPHICAL ABSTRACT

if (figure_SNOW_VARIABLES_TIMESERIES){
  results_property_lvl = fread(paste0(path_proj,"\\Results\\results_property_lvl.csv"))
  
  results_property_lvl_long = 
    melt(results_property_lvl[period %in% 1:period_last,
                              .(period, management, csl_m2, break_p_proc, n_tot, dh_ratio)], 
         id.vars = c("period", "management"), 
         variable.name = "variable", 
         value.name = "val") %>%
    .[, variable := factor(variable, levels = c("n_tot", "dh_ratio", "csl_m2", "break_p_proc"))]
  
  variable_labels = as_labeller(c(
    csl_m2 = "Critical~snow~load~(Kg*~m^{-2})",
    break_p_proc = "Snow~damage~probability~('%'*~yr^{-1})",
    n_tot = "Number~of~trees~(ha^{-1})",
    dh_ratio = "Diameter (cm) / Height (m)~ratio"
  ), label_parsed)
  
  val_by_manag = 
    results_property_lvl_long[period %in% 5:10,
                              .(val_mean = mean(val, na.rm = TRUE)),
                              by = .(variable, management)]%>%
    .[, label := fifelse(variable == "csl_m2",       sprintf("%.2f", val_mean),
                         fifelse(variable == "break_p_proc", sprintf("%.3f", val_mean),
                                 fifelse(variable == "n_tot",        sprintf("%.0f", val_mean),
                                         fifelse(variable == "dh_ratio",     sprintf("%.3f", val_mean), NA_character_))))]
  
  fwrite(val_by_manag, paste0(path_proj,"\\Results\\val_by_manag.csv"))
  
  figure_SNOW_VARIABLES_TIMESERIES = 
    ggplot(results_property_lvl_long, 
           aes(x = period, 
               y = val, 
               color = management)) +
    facet_wrap(~ variable, ncol = 2, scales = "free", labeller = labeller(variable = variable_labels)) + 
    annotate("rect", xmin = 5, xmax = 10, ymin = -Inf, ymax = Inf, fill = "gray90") +
    geom_line() +
    geom_segment(data = val_by_manag, 
                 aes(x = 5, xend = 10, y = val_mean, yend = val_mean, color = management),
                 linetype = "dashed",
                 show.legend = FALSE) +
    geom_text(data = val_by_manag, 
              aes(x = 5.1, y = val_mean, label = label, color = management),
              hjust = 0, vjust = -0.5,
              size = 3,
              show.legend = FALSE) +
    scale_y_continuous(labels = label_number())+
    labs(x = "Period") +
    theme_bw()+
    scale_x_continuous(breaks = 1:period_last_value) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.title.y = element_blank())+
    scale_color_manual(
      values = c("standard" = "black", "fritzoe" = "#00BFFF"),  
      labels = c("standard" = "Standard management", "fritzoe" = "Reduced density management")
    )
  
  ggsave(
    filename = paste0(path_proj,"\\Figures\\figure_SNOW_VARIABLES_TIMESERIES.jpg"), 
    plot = figure_SNOW_VARIABLES_TIMESERIES, 
    width = twocol, height = maxHeight*.6, units = "mm",dpi = 1000)
} # FIGURE 5

if (figure_CSL_SDP_HISTO_5_10){
  results_substand_lvl = fread(paste0(path_proj,"\\Results\\results_substand_lvl.csv"))
  val_by_manag = 
    fread(paste0(path_proj,"\\Results\\val_by_manag.csv")) %>%
    .[, variable := factor(variable, levels = c("n_tot", "dh_ratio", "csl_m2", "break_p_proc"))]
  
  results_substand_lvl_long = 
    melt(results_substand_lvl[period %in% 5:10,
                              .(period, management, 
                                csl_m2 = csl, 
                                break_p_proc = ifelse(break_p>0,100*break_p,NA), 
                                n_tot = ifelse(n + reg_n > 0, n + reg_n, NA), 
                                dh_ratio = d/h)],
         id.vars = c("period", "management"), 
         variable.name = "variable", 
         value.name = "val") %>%
    .[, variable := factor(variable, levels = c("n_tot", "dh_ratio", "csl_m2", "break_p_proc"))] 
  
  variable_labels = as_labeller(c(
    csl_m2 = "Critical~snow~load~(Kg*~m^{-2})",
    break_p_proc = "Snow~damage~probability~('%'*~yr^{-1})",
    n_tot = "Number~of~trees~(ha^{-1})",
    dh_ratio = "Diameter(cm) / Height(m)~ratio"
  ), label_parsed)
  
  figure_CSL_SDP_HISTO_5_10 =
    ggplot(results_substand_lvl_long[period %in% 5:10],
           aes(x = val, fill = management)) +
    geom_histogram(aes(y = after_stat(100*count / sum(count))), 
                   position = "identity", 
                   alpha = 0.5, 
                   bins = 30) +
    facet_wrap(~ variable, scales = "free", labeller = variable_labels) +
    geom_text_repel(data = val_by_manag,
                    aes(x = val_mean, y = Inf, label = label, color = management),
                    nudge_y = -5,  # push down a bit from top
                    direction = "y",
                    hjust = 0,
                    vjust = 0,
                    size = 3,
                    segment.color = NA,  
                    show.legend = FALSE)+
    # Vertical dashed line for mean
    geom_vline(data = val_by_manag,
               aes(xintercept = val_mean, color = management),
               linetype = "dashed",
               show.legend = FALSE) +
    scale_fill_manual(
      values = c("standard" = "black", "fritzoe" = "#00BFFF"),
      labels = c("standard" = "Standard management", "fritzoe" = "Reduced density management")
    ) +
    scale_color_manual(
      values = c("standard" = "black", "fritzoe" = "#00BFFF")
    ) +
    labs(x = NULL, y = "Freqency (%)") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.text = element_text(size = 11)
    )
  
  ggsave(
    filename = paste0(path_proj,"\\Figures\\figure_CSL_SDP_HISTO_5_10.jpg"), 
    plot = figure_CSL_SDP_HISTO_5_10, 
    width = twocol, height = maxHeight*.4, units = "mm",dpi = 1000)
  
  figure_CSL_SDP_ECDF_5_10 =
    ggplot(results_substand_lvl_long[period %in% 5:10],
           aes(x = val, color = management)) +
    stat_ecdf(geom = "step", aes(y = after_stat(..y..) * 100), size = 0.8) +
    facet_wrap(~ variable, scales = "free", labeller = variable_labels) +
    
    geom_vline(data = val_by_manag,
               aes(xintercept = val_mean, color = management),
               linetype = "dashed",
               show.legend = FALSE) +
    
    geom_text_repel(data = val_by_manag,
                    aes(x = val_mean, y = 100, label = label, color = management),
                    direction = "y",
                    nudge_y = -5,
                    hjust = 0,
                    vjust = 0,
                    size = 3,
                    segment.color = NA,
                    show.legend = FALSE) +
    
    scale_color_manual(
      values = c("standard" = "black", "fritzoe" = "#00BFFF"),
      labels = c("standard" = "Standard management", "fritzoe" = "Reduced density management")
    ) +
    labs(x = NULL, y = "Cumulative frequency (%)") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.text = element_text(size = 11)
    )
} # FIGURE 7

if (figure_CSL_vs_BREAKP){
  results_substand_lvl = fread(paste0(path_proj,"\\Results\\results_substand_lvl.csv"))
  
  scatter_data = 
    results_substand_lvl[period %in% 5:10,
                         .(period, management, 
                           csl_m2 = csl, 
                           emsl = max_sl,
                           break_p_proc = 100 * break_p)] %>%
    .[,management:=factor(management,levels = c("standard","fritzoe"))]
  
  # Create the plot
  figure_CSL_vs_BREAKP = 
    ggplot(scatter_data, aes(x = csl_m2, y = break_p_proc)) +
    geom_hex(aes(fill = ..count.. / tapply(..count.., ..PANEL.., sum)[..PANEL..]), bins = 40) +
    facet_wrap(~ management, labeller = as_labeller(c(
      standard = "Standard management",
      fritzoe = "Reduced density management"
    ))) +
    scale_fill_viridis_c(option = "D", name = "Freq. (%)", labels = scales::percent) +
    labs(x = "Critical snow load (Kg m⁻²)",
         y = "Snow damage probability (% yr⁻¹)") +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) 
  
  figure_CSL_vs_EMSL = 
    ggplot(scatter_data, aes(x = csl_m2, y = emsl)) +
    geom_hex(aes(fill = ..count.. / tapply(..count.., ..PANEL.., sum)[..PANEL..]), bins = 40) +
    facet_wrap(~ management, labeller = as_labeller(c(
      standard = "Standard management",
      fritzoe = "Reduced density management"
    ))) +
    scale_fill_viridis_c(option = "D", name = "Freq. (%)", labels = scales::percent) +
    labs(x = "Critical snow load (Kg m⁻²)",
         y = "Expected maximum snow load (Kg m⁻²)") +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = paste0(path_proj,"\\Figures\\figure_CSL_vs_BREAKP.jpg"), 
    plot = figure_CSL_vs_BREAKP, 
    width = twocol, height = maxHeight*.4, units = "mm",dpi = 1000)
  
  ggsave(
    filename = paste0(path_proj,"\\Figures\\figure_CSL_vs_EMSL.jpg"), 
    plot = figure_CSL_vs_EMSL, 
    width = twocol, height = maxHeight*.4, units = "mm",dpi = 1000)
} # FIGURE 8

if (figure_CSL_SDP_BY_SI){
  results_substand_lvl = fread(paste0(path_proj,"\\Results\\results_substand_lvl.csv"))
  
  val_by_manag_si = 
    results_substand_lvl[period %in% 5:10,
                         .(csl = weighted.mean(csl, area, na.rm = TRUE),
                           break_p = 100*weighted.mean(break_p, area, na.rm = TRUE)),
                         by = .(si, management)]
  
  val_by_manag_long = melt(val_by_manag_si,
                           id.vars = c("si", "management"),
                           measure.vars = c("csl", "break_p"),
                           variable.name = "variable",
                           value.name = "val")
  
  figure_CSL_SDP_BY_SI = 
    ggplot(val_by_manag_long, aes(x = si, y = val, color = management)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ variable, scales = "free_y", labeller = as_labeller(c(
      csl = "Critical snow load (Kg m⁻²)",
      break_p = "Snow damage probability (% yr⁻¹)"
    ))) +
    scale_color_manual(
      values = c("standard" = "black", "fritzoe" = "#00BFFF"),  
      labels = c("standard" = "Standard management", "fritzoe" = "Reduced density management")
    ) +
    scale_x_continuous(breaks = unique(val_by_manag_long$si)) +
    labs(x = "Site index (m)", y = NULL, color = NULL) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 11)
    )
  
  ggsave(
    filename = paste0(path_proj,"\\Figures\\figure_CSL_SDP_BY_SI.jpg"), 
    plot = figure_CSL_SDP_BY_SI, 
    width = twocol, height = maxHeight*.4, units = "mm",dpi = 1000)
} # FIGURE 6

if (figure_PRODUCTION_TABLE_TIMESERIES){
  if (!exists("results_property_lvl")) results_property_lvl = fread(paste0(path_proj,"\\Results\\results_property_lvl.csv"))
  
  prod_table = results_property_lvl[period %in% 1:10,.(period, management, 
                                                       area_harvest, area_thinning, area_tending,
                                                       value_harvest_m3 = nok_to_eur*value_harvest_m3, 
                                                       cost_harvest_m3 = nok_to_eur*cost_harvest_m3,
                                                       vol_standing_tot)]
  
  results_table_long = melt(prod_table, id.vars = c("period", "management"), variable.name = "variable", value.name = "val")
  results_table_long$management <- factor(results_table_long$management, 
                                          levels = c("standard", "fritzoe"))
  
  facet_labels = c(
    area_harvest = bquote("Final harvest area ("*ha*")"),
    area_thinning = bquote("Thinning area ("*ha*")"),
    area_tending = bquote("Tending area ("*ha*")"),
    
    value_harvest_m3 = bquote("Timber value harvest ("*EUR~m^-3*")"),
    cost_harvest_m3 = bquote("Cost harvest ("*EUR~m^-3*")"),
    vol_standing_tot = bquote("Standing volume ("*m^3*")")
  )
  
  vlabeller = function (variable, value) {return(facet_labels[value])}
  
  figure = 
    ggplot(results_table_long, aes(x = period, y = val, color = management)) +
    geom_line() +
    scale_y_continuous(labels = label_number()) +
    labs(x = "Period") +
    theme_bw() +
    facet_wrap(~variable, scales = "free_y", labeller = vlabeller, ncol = 2, dir = "v") + 
    scale_x_continuous(breaks = 1:period_last) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.title.y = element_blank()) +
    scale_color_manual(
      values = c("standard" = "black", "fritzoe" = "#00BFFF"),  
      labels = c("standard" = "Standard management", "fritzoe" = "Reduced density management")
    ) 
  
  ggsave(
    filename = paste0(path_proj,"\\Figures\\figure_PRODUCTION_TABLE_TIMESERIES.jpg"), 
    plot = figure, 
    width = twocol, height = maxHeight*.5, units = "mm",dpi = 1000)
} # FIGURE 4