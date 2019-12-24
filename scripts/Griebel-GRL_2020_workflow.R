rm(list = ls())

#define location of source script that contains functions
source('scripts/eddy4r.fast_functions.R')

#check duration
start.time <- Sys.time()

#call funtion to calculate traditional, standardized, space-equitable and space-time-equitable budgets
results=makeFast('input/AU-Cum_example.csv')

#alternative  example for using latent heat flux instead of carbon fluxes:
#results=makeFast('input/AU-Cum_example.csv',targetCol="LE_F_MDS",QcCol="LE_F_MDS_QC",normalize=FALSE)

#check duration
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

