# high-dimentional time series project of Jesse van der Ende, Amine Ziddi
# and Lena Monsch


# read the datasets

# first option to work with:
loc <- read.csv("locations032017.csv")
air <- read.csv("TaiwanAirBox032017.csv", col.names = c("time", paste0("col",2:517)))
# hourly PM2.5 measurements taken in March 2017 in Taiwan (in micrograms per cubic meters, μg/m2)
# 516 series each with 31 x 24 = 744 observations (i.e. k= 516 and T = 744)

head(loc)
head(air[,1:5])



# second option to work with:
elec <- read.csv("PElectricity1344.csv", col.names = paste0("col",1:1344))
# weekly series of electricity price each hour of each day during 678 weeks in 
# eight regions of New England (7 × 24 × 8 = 1344)
# Note that here k > T

# L: but how are we supposed to know which data belongs to which region
head(elec[,1:5])

