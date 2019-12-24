library(miscTools)

#define constants and variables
flux_grams_HH=1*12.0107/1000000*60*30 #half hourly conversion from umol/m2/s to grams/m2/30mins
flux_grams_HR=1*12.0107/1000000*60*60 #hourly conversion from umol/m2/s to grams/m2/hour
obs_year_HH=48*365
obs_year_HR=24*365
time_HH = 60*30
time_HR = 60*60

#clean data frame, assign timestamp, establish columnames for flux, wind direction and year, remove missing values
clean_df = function (df, QcCol="NEE_VUT_REF_QC",wdCol="WD",TimestampCol="TIMESTAMP_START") {
  #align column names and declare year
  date=as.character(df[[TimestampCol]])
  # df$ws<-df$WS
  df$wd<-df[[wdCol]]
  df$Year<-substr(date, 1,4) 
  df$time_string<-substr(date, 9,12)
  df$time<-as.numeric(df$time_string)
 
  
  #add column that contains the eight wind sectors
  if ("wd" %in% colnames(df)) {
    df$eight_sec <- with (df, ifelse (is.na(wd), NA,
                              ifelse (wd == -9999, NA ,
                              ifelse (wd < 22.5 | wd >337.5,"N", 
                              ifelse (wd >22.5 & wd <67.5,"NE", 
                              ifelse (wd >67.5 & wd <112.5,"E",        
                              ifelse (wd >112.5 & wd <157.5,"SE", 
                              ifelse (wd >157.5 & wd <202.5,"S",  
                              ifelse (wd >202.5 & wd <270-22.5,"SW", 
                              ifelse (wd >270-22.5 & wd<270+22.5,"W", "NW")))))))))) 
   
    df$time_int <- with (df, ifelse (is.na(time), NA,
                              ifelse (time < 300, "0-3", 
                              ifelse (time >= 300 & time <600,"3-6", 
                              ifelse (time >= 600 & time <900,"6-9",        
                              ifelse (time >= 900 & time <1200,"9-12",
                              ifelse (time >= 1200 & time <1500,"12-15",
                              ifelse (time >= 1500 & time <1800,"15-18",
                              ifelse (time >= 1800 & time <2100,"18-21","21-00")))))))))
    

  
    ##remove missing values for wind direction and only select observational data
    df[df == -9999] <- NA
    ##only select observational data based on QC flag = 0
    df = subset(df, !is.na(wd) & df[[QcCol]] == 0)
    return(df) 
  } else {
    print("check if wind direction is present")
  }
  
}

#function to calculate traditional budget
calculate_uncorrected <- function (df, year,targetCol="NEE_VUT_REF",frequency="30min") {
  #subset data frame by year
  subsetted = subset(df, df$Year==year)
  #calculate traditional budget based on simple sum of all observations and convert units 
  if(frequency == "30min") {
    budget<-sum(subsetted[[targetCol]])*time_HH
  } else if (frequency == "60min") {
    budget<-sum(subsetted[[targetCol]])*time_HR
  } else {
    print("only 30min or 60min frequency supported")
  }
  return(budget)
}


#function to calculate standardized budgets based on the average wind pattern for all observation years
calculate_standardized <- function (df, year,targetCol="NEE_VUT_REF",frequency="30min") {
  
  no_years = length(unique(df$Year))
  
  #subset into eight wind sectors
  df_N = subset(df, df$eight_sec == "N")
  df_NE = subset(df, df$eight_sec == "NE")
  df_E = subset(df, df$eight_sec == "E")
  df_SE = subset(df, df$eight_sec == "SE")
  df_S = subset(df, df$eight_sec == "S")
  df_SW = subset(df, df$eight_sec == "SW")
  df_W = subset(df, df$eight_sec == "W")
  df_NW = subset(df, df$eight_sec == "NW")
  
  #subset dataframe by year [i]
  subsetted = subset(df, df$Year==year)
  
  #subset dataframe of year [i] into dataframes for each wind sector
  subsetted_N = subset(subsetted, subsetted$eight_sec == "N")
  subsetted_NE = subset(subsetted, subsetted$eight_sec == "NE")
  subsetted_E = subset(subsetted, subsetted$eight_sec == "E")
  subsetted_SE = subset(subsetted, subsetted$eight_sec == "SE")
  subsetted_S = subset(subsetted, subsetted$eight_sec == "S")
  subsetted_SW = subset(subsetted, subsetted$eight_sec == "SW")
  subsetted_W = subset(subsetted, subsetted$eight_sec == "W")
  subsetted_NW = subset(subsetted, subsetted$eight_sec == "NW")
  
  #carbon budget adjustment as outlined in Griebel et al. 2016:
  #step 1: define the average wind pattern as standardized wind sector contribution to total observations of all years
  N_frac = length(df_N[[targetCol]])/no_years
  NE_frac = length(df_NE[[targetCol]])/no_years
  E_frac = length(df_E[[targetCol]])/no_years
  SE_frac = length(df_SE[[targetCol]])/no_years
  S_frac = length(df_S[[targetCol]])/no_years
  SW_frac = length(df_SW[[targetCol]])/no_years
  W_frac = length(df_W[[targetCol]])/no_years
  NW_frac = length(df_NW[[targetCol]])/no_years

  #step 2: calculate the mean carbon uptake for each sector and multiply with the average contribution of each sector
  N_cor  <- N_frac*mean(subsetted_N[[targetCol]])
  NE_cor <- NE_frac*mean(subsetted_NE[[targetCol]])
  E_cor  <- E_frac*mean(subsetted_E[[targetCol]])
  SE_cor <- SE_frac*mean(subsetted_SE[[targetCol]])
  S_cor  <- S_frac*mean(subsetted_S[[targetCol]])
  SW_cor <- SW_frac*mean(subsetted_SW[[targetCol]])
  W_cor  <- W_frac*mean(subsetted_W[[targetCol]])
  NW_cor <- NW_frac*mean(subsetted_NW[[targetCol]])  
  
  #step 3: integrate across all sectors
   if(frequency == "30min") {
    new_budget_year <- sum(N_cor, NE_cor, E_cor, SE_cor, S_cor, SW_cor, W_cor, NW_cor, na.rm = TRUE)*time_HH
  } else if (frequency == "60min") {
    new_budget_year <- sum(N_cor, NE_cor, E_cor, SE_cor, S_cor, SW_cor, W_cor, NW_cor, na.rm = TRUE)*time_HR
  } else {
    print("only 30min or 60min frequency supported")
  }
  return(new_budget_year)
}  


#function to calculate space equitable budgets where every wind sector contributes equally
calculate_space_equitable <- function (df, year,targetCol="NEE_VUT_REF",frequency="30min",normalize=TRUE) {
  # 
  # #subset into eight wind sectors
  # df_N = subset(df, df$eight_sec == "N")
  # df_NE = subset(df, df$eight_sec == "NE")
  # df_E = subset(df, df$eight_sec == "E")
  # df_SE = subset(df, df$eight_sec == "SE")
  # df_S = subset(df, df$eight_sec == "S")
  # df_SW = subset(df, df$eight_sec == "SW")
  # df_W = subset(df, df$eight_sec == "W")
  # df_NW = subset(df, df$eight_sec == "NW")
  
  #subset dataframe by year [i]
  subsetted = subset(df, df$Year==year)
  
  #step 1: subset dataframe of year [i] into dataframes for each wind sector
  subsetted_N = subset(subsetted, subsetted$eight_sec == "N")
  subsetted_NE = subset(subsetted, subsetted$eight_sec == "NE")
  subsetted_E = subset(subsetted, subsetted$eight_sec == "E")
  subsetted_SE = subset(subsetted, subsetted$eight_sec == "SE")
  subsetted_S = subset(subsetted, subsetted$eight_sec == "S")
  subsetted_SW = subset(subsetted, subsetted$eight_sec == "SW")
  subsetted_W = subset(subsetted, subsetted$eight_sec == "W")
  subsetted_NW = subset(subsetted, subsetted$eight_sec == "NW")
  
  #step 2: calculate mean carbon uptake for each sector and each year 
  if(frequency == "30min") {
    N_cor  <- (48*365*0.125)*mean(subsetted_N[[targetCol]])
    NE_cor <- (48*365*0.125)*mean(subsetted_NE[[targetCol]])
    E_cor  <- (48*365*0.125)*mean(subsetted_E[[targetCol]])
    SE_cor <- (48*365*0.125)*mean(subsetted_SE[[targetCol]])
    S_cor  <- (48*365*0.125)*mean(subsetted_S[[targetCol]])
    SW_cor <- (48*365*0.125)*mean(subsetted_SW[[targetCol]])
    W_cor  <- (48*365*0.125)*mean(subsetted_W[[targetCol]])
    NW_cor <- (48*365*0.125)*mean(subsetted_NW[[targetCol]]) 
    equit_budget_year <- sum(N_cor, NE_cor, E_cor, SE_cor, S_cor, SW_cor, W_cor, NW_cor, na.rm = TRUE)*time_HH
   } else if (frequency == "60min") {
     N_cor  <- (24*365*0.125)*mean(subsetted_N[[targetCol]])
     NE_cor <- (24*365*0.125)*mean(subsetted_NE[[targetCol]])
     E_cor  <- (24*365*0.125)*mean(subsetted_E[[targetCol]])
     SE_cor <- (24*365*0.125)*mean(subsetted_SE[[targetCol]])
     S_cor  <- (24*365*0.125)*mean(subsetted_S[[targetCol]])
     SW_cor <- (24*365*0.125)*mean(subsetted_SW[[targetCol]])
     W_cor  <- (24*365*0.125)*mean(subsetted_W[[targetCol]])
     NW_cor <- (24*365*0.125)*mean(subsetted_NW[[targetCol]]) 
    
     #step 3: sum across all sectors 
    equit_budget_year <- sum(N_cor, NE_cor, E_cor, SE_cor, S_cor, SW_cor, W_cor, NW_cor, na.rm = TRUE)*time_HR
  } else {
    print("only 30min or 60min frequency supported")
  }
  #check if normalization to number of observations is true  
  if (normalize==TRUE){
    if(frequency == "30min"){
    equit_budget_year=equit_budget_year/(obs_year_HH)*dim(subset(df, df$Year==year))[1]
  } else {
    equit_budget_year=equit_budget_year/(obs_year_HR)*dim(subset(df, df$Year==year))[1]
  }}
  return(equit_budget_year)
}  


#function to calculate space-time equitable budgets where every wind sector contributes equally
calculate_space_time_equitable <- function (df, year,targetCol="NEE_VUT_REF",frequency="30min",normalize=TRUE) {
 #subset data frame by year
  subsetted = subset(df, df$Year==year)
  wind_sectors<-list("N","NE","E","SE","S","SW","W","NW")
  time_sectors<-list("0-3","3-6","6-9","9-12","12-15","15-18","18-21","21-00")
  
  #create empty data frame that will be filled with matrix that has columns as sectors and rows as time periods
  spti<-data.frame()
  
  for (w in wind_sectors) {
     df_wind<-subset(subsetted,subsetted$eight_sec==w)
     for (t in time_sectors) {
        df_time<-subset(df_wind,df_wind$time_int==t)
        spti[t,w]=mean(df_time[[targetCol]])
     }
  }
  
  #average by time first (row average) and then by space (column average)
  space_means<-rowMeans(spti,na.rm=TRUE)
  space_time_means<-mean(space_means,na.rm=TRUE)
  if(frequency == "30min"){
    space_time_equit_budget_year <- space_time_means*time_HH*obs_year_HH
  } else {
    space_time_equit_budget_year <- space_time_means*time_HR*obs_year_HR
  }
  
  #check if normalization to number of observations is true 
  if (normalize==TRUE){
    if(frequency == "30min"){
  space_time_equit_budget_year <- space_time_equit_budget_year/(obs_year_HH)*dim(subset(df, df$Year==year))[1]
    } else {
           space_time_equit_budget_year <- space_time_equit_budget_year/(obs_year_HR)*dim(subset(df, df$Year==year))[1]
    }}
  return(space_time_equit_budget_year)
}

  
#this function requires an input file which is the path to the original .csv file
makeFast = function (infile, 
                     windirCol='WD', 
                     QcCol="NEE_VUT_REF_QC",
                     targetCol="NEE_VUT_REF",
                     frequency="30min",
                     normalize=TRUE) {
  
  # create empty DataFrame to store results and specify dst_path
  results = data.frame(matrix(ncol = 6, nrow = 0))
  colnames(results)= c('year', 'traditional_budget', 'standardized_budget', 'space_equitable_budget', 'space_time_equitable_budget','n_obs')
  
  
  # calculate traditional, standardized, space-equitable and space-time-equitable budgets
  
  #read file
  df = read.csv(infile,stringsAsFactors = FALSE)
  
  #clean data frame, assign timestamp, establish columnames for flux, wind direction and year, remove missing values
  df = clean_df(df,QcCol=QcCol)
  
  #identify unique years
  years = unique(df$Year)
  
  # set initial row number
  i=1
  
  # loop through each year
  for (y in years){
    # calculate uncorrected and corrected budgets
    budget = calculate_uncorrected(df, y,targetCol=targetCol,frequency=frequency) 
    budget_corr = calculate_standardized(df, y,targetCol=targetCol,frequency=frequency) 
    budget_equit = calculate_space_equitable(df, y,targetCol=targetCol,frequency=frequency,normalize=normalize)
    budget_space_time_equit = calculate_space_time_equitable(df, y,targetCol=targetCol,frequency=frequency,normalize=normalize)
    # calculate number of observations results are based on
    n_obs = dim(subset(df, df$Year==y))[1]
    
    # generate rows for site and year y 
    annual_result = c(y, budget, budget_corr, budget_equit, budget_space_time_equit, n_obs)
    results <- insertRow(as.matrix(results), i, annual_result)
    
    # increase row number by 1
    i=i+1
  }
  
  return(results)

}
