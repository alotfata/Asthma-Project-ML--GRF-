#install.packages("reldist")
#install.packages("gplm")
#install.packages("jtools")
#load libraries
require(caTools)
library(tidycensus)
library(tidyverse)
library(dplyr)
library(stringr)
library(usmap)
library(readxl)
library(sf)
library(reshape2)
library(cdlTools)
library(tigris)
library(reldist)
library(psych)
library(car)
library(spdep)
#library(INLA)
library(tmap)
library(gplm)
library(MASS)
library(raster)
library(jtools)
library(RColorBrewer)

#-----------------------------------------------------------------------------------------
# #preliminaries

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set workspace dir
 
census_api_key("767b7c8f877abf0396b842d5a3f5e92aa5ba05c1") #get census API key!
 
 
# #-----------------------------------------------------------------------------------------
#response variable

#cook county census tracts
cook <- read.csv("../data/cook_tracts.csv")
cook$GEOID <- as.character(cook$GEOID)
cook <- subset(cook, select = c(GEOID, ALAND))  #land area

#asthma data
asthma <- read.csv("../data/PLACES__Local_Data_for_Better_Health__Census_Tract_Data_2020_release.csv") %>%
  subset(., Year == 2018 & StateDesc == "Illinois" & CountyName == "Cook" & Measure == "Current asthma among adults aged >=18 years",
         select=c(LocationName, Data_Value))
asthma$LocationName <- as.character(asthma$LocationName)
write.csv(asthma,"../data/asthma.csv")

# #-----------------------------------------------------------------------------------------
# #tract-level predictor variables 

#tree canopy
canopy <- read.csv("../data/TreeCanopyArea.csv") %>%
  subset(., select = -c(ZONE_CODE))
canopy$FIPS <- as.character(canopy$FIPS)

#demographics
Tot_pop <- subset(get_acs(geography = "tract",
                          variables = c(var = "S0101_C01_001E"),
                          state = "IL",
                          county = "Cook",
                          year = 2018),
                  select = -c(NAME,variable, moe))

Pop_m <- subset(get_acs(geography = "tract",
                        variables = c(var = "B01001_002"),
                        state = "IL",
                        county = "Cook",
                        year = 2018),
                select = -c(NAME,variable, moe))
Pop_m$estimate <- Pop_m$estimate / Tot_pop$estimate
                   
Pop_f <- Pop_m
Pop_f$estimate <- 1 - Pop_f$estimate 
                   

# #race, ethnicity variables
Pop_white <- get_acs(geography = "tract",
                            variables = c(var = "B01001A_001"),
                            state = "IL",
                            county = "Cook",
                            year = 2018)
Pop_white$ratio <- Pop_white$estimate / Tot_pop$estimate
Pop_white <- subset(Pop_white, select=c("GEOID", "ratio"))


# 
Pop_black <- get_acs(geography = "tract",
                            variables = c(var = "B01001B_001"),
                            state = "IL",
                            county = "Cook",
                            year = 2018)
Pop_black$ratio <- Pop_black$estimate / Tot_pop$estimate
Pop_black <- subset(Pop_black, select=c("GEOID", "ratio"))

Pop_asian <- get_acs(geography = "tract",
                           variables = c(var = "B01001D_001"),
                           state = "IL",
                           county = "Cook",
                           year = 2018)
Pop_asian$ratio <- Pop_asian$estimate / Tot_pop$estimate
Pop_asian <- subset(Pop_asian, select=c("GEOID", "ratio"))

Pop_hisp <- get_acs(geography = "tract",
                           variables = c(var = "DP05_0071E"),
                           state = "IL",
                           county = "Cook",
                           year = 2018)
Pop_hisp$ratio <- Pop_hisp$estimate / Tot_pop$estimate
Pop_hisp <- subset(Pop_hisp, select=c("GEOID", "ratio"))


# #hospital density
# h_pts <- st_read("../data/CMAP_MedicalCenter/MedicalCenter_CMAP.shp") %>% #read hospital points
#   st_coordinates() %>%
#   as.data.frame()
# 
# h <- bandwidth.scott(h_pts, kernel = "gaussian", product = TRUE)  #compute bandwidth (using scott's method)
# 
# h_kde <- kde2d(h_pts$X, h_pts$Y, h = h, n = c(1646,2045), lims = c(range(h_pts$X), range(h_pts$Y))) %>% #compute density
#   raster()
# 
# geom <- st_read("../data/cook_tracts.shp") %>%  #read tract geometries
#   st_transform(., crs = 3435)
# 
# hosp <- extract(h_kde, geom, fun = mean, na.rm=TRUE, df=TRUE, cellnumbers=TRUE) %>% #average density for each tract
#   cbind(geom$GEOID, .) %>%
#   subset(., select = c("geom$GEOID","layer"))
# colnames(hosp) <- c("GEOID", "MEDCENTER")
# save(hosp,file="hosp.Rda")

load("hosp.Rda")

#social vulnerability
vulner <- read.csv("../data/CDC SVI Data Illinois.csv") %>%
  subset(., select = c(FIPS, RPL_THEME1, RPL_THEME2, RPL_THEME3, RPL_THEME4))
vulner$FIPS <- as.character(vulner$FIPS)

#pm2.5
# PM25 <- read.csv("../data/Daily_Census_Tract-Level_PM2.5_Concentrations__2016.csv")
# PM25_2 <- subset(PM25, substr(ctfips, 1, 5) == 17031, select=c(ctfips, date, DS_PM_pred, DS_PM_stdd))
# PM25_sp <- aggregate(PM25_2$DS_PM_pred, by=list(Category=PM25_2$ctfips), FUN=median)
# colnames(PM25_sp) <- c("FIPS", "PM25")
# PM25_sp$FIPS <- as.character(PM25_sp$FIPS)
# # plot
# g <- ggplot(data = PM25_2, aes(x=date, y=DS_PM_pred)) + geom_line(aes(colour=ctfips))
# png("PM25.png")
# print(g)
# dev.off()
# save(PM25_sp, file="PM25.Rda")

load("PM25.Rda")

#-----------------------------------------------------------------------------------------
#join

df <- left_join(cook, asthma, by = c("GEOID" = "LocationName")) %>%
  left_join(., PM25_sp, by = c("GEOID" = "FIPS")) %>%
  left_join(., vulner, by = c("GEOID" = "FIPS")) %>%
  left_join(., canopy, by = c("GEOID" = "FIPS")) %>%
  left_join(., Tot_pop, by = "GEOID") %>%
  left_join(., Pop_f, by = "GEOID") %>%
  left_join(., Pop_white, by = "GEOID") %>%
  left_join(., Pop_black, by = "GEOID") %>%
  left_join(., Pop_hisp, by = "GEOID") %>%
  left_join(., Pop_asian, by = "GEOID") %>%
  left_join(., hosp, by = "GEOID")

df$RATIO <- df$AREA/df$ALAND  #canopy ratio
df$estimate.x <- df$estimate.x/df$ALAND   #pop density

colnames(df) <- c("GEOID", "ALAND", "ASTHMA", "PM25", "RPL_THEME1", "RPL_THEME2", "RPL_THEME3", "RPL_THEME4", "AREA", "POPDENS",
                  "POPF", "POPW", "POPB", "POPH", "POPA", "MEDCENTER", "RATIO")

df <- df[complete.cases(df), ]

#df$INCOME <- (df$INCOME - min(df$INCOME))/ (max(df$INCOME) - min(df$INCOME))

df <- subset(df, select = -c(ALAND))



#-----------------------------------------------------------------------------------------
#correlation matrix

df_subset <- df %>% 
  subset(., select = c(ASTHMA,RPL_THEME2,RPL_THEME3,RPL_THEME4,RATIO,MEDCENTER,PM25)) %>%
  rename(., "SVI Theme 2" = "RPL_THEME2") %>%
  rename(., "SVI Theme 3" = "RPL_THEME3") %>%
  rename(., "SVI Theme 4" = "RPL_THEME4") %>%
  rename(., "Tree Cover" = "RATIO") %>%
  rename(., "Medical Center" = "MEDCENTER")


corr <- cor(df_subset, method = c("pearson", "kendall", "spearman"))
round(corr, 2)


pairs.panels(df_subset,
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)


#-----------------------------------------------------------------------------------------
#first model

geom <- st_read("../data/cook_tracts.shp")
df <- left_join(df, geom, by = "GEOID")

#drop columns
drops <- c("ALAND","AWATER","Shape__Are","Shape__Len","GEOID_2","Area") #drop columns
df <- df[ , !(names(df) %in% drops)] #drop columns
df$GEOID <- as.character(df$GEOID)

#scale variables
normalit<-function(m){
  (m - min(m))/(max(m)-min(m))
}

df$RPL_THEME1 <- normalit(df$RPL_THEME1)
df$RPL_THEME2 <- normalit(df$RPL_THEME2)
df$RPL_THEME3 <- normalit(df$RPL_THEME3)
df$RPL_THEME4 <- normalit(df$RPL_THEME4)

df$POPF <- normalit(df$POPF)
df$POPW <- normalit(df$POPW)
df$POPB <- normalit(df$POPB)
df$POPH <- normalit(df$POPH)
df$POPA <- normalit(df$POPA)
df$MEDCENTER <- normalit(df$MEDCENTER)
df$RATIO <- normalit(df$RATIO)
df$PM25 <- normalit(df$PM25)

#OLS
#mod6 <- lm(ASTHMA~RPL_THEME2+RPL_THEME4+RATIO+MEDCENTER+PM25+POPW, data = df)
mod7 <- lm(ASTHMA~RPL_THEME2+RPL_THEME3+RPL_THEME4+RATIO+MEDCENTER+PM25, data = df)
#summ(mod6)  #summary
summ(mod7)
# sort(vif(mod6))  #variance inflation factors
sort(vif(mod7))  #variance inflation factors
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
#plot(mod6)  #residual plots
plot(mod7)  #residual plots


#Spatial Analysis of Residuals of first model
#-------------------------------------------------------------------------------------------------------------

#spatial neighbors
w <- poly2nb(df$geometry, row.names=df$GEOID)
summary(w)

#wm <- as.matrix(nblag(w, maxlag = 3)[[3]])


wm <- nb2mat(w, style='S')

rwm <- mat2listw(wm, style='W')

lm.morantest(mod7, rwm, alternative="two.sided")

#map residuals
df <- df %>%
  cbind(., mod7$residuals) %>%
  cbind(., mod7$fitted.values)

names(df)[names(df) == "mod7$residuals"] <- "modres"
names(df)[names(df) == "mod7$fitted.values"] <- "modfit"


library(RColorBrewer)
pal <- brewer.pal(7, "OrRd") # we select 7 colors from the palette
class(pal)

plot(st_as_sf(df)["modres"],
     main = "Model 1 residuals",
     breaks = "quantile", nbreaks = 7,
     pal = pal)

#Lagrgange Multiplier tests
lm.LMtests(mod7, rwm, test = c("LMerr","LMlag","RLMerr","RLMlag","SARMA"))

#-----------------------------------------------------------------------------------------
#spatial lag model

mod7_lag <- lagsarlm(ASTHMA~RPL_THEME2+RPL_THEME3+RPL_THEME4+RATIO+MEDCENTER+PM25, data = df, rwm)

#mod7_err <- errorsarlm(ASTHMA~RPL_THEME2+RPL_THEME4+RATIO+MEDCENTER+PM25+POPW, data = df, rwm)

summary(mod7_lag)
#summary(mod6_err)

#direct, total and indirect impacts
W <- as(rwm, "CsparseMatrix")
trMC <- trW(W, type="MC")
im<-impacts(mod1_lag, tr=trMC, R=100)
sums<-summary(im,  zstats=T)
#To print the coefficients
data.frame(sums$res)
data.frame(sums$pzmat)

#-----------------------------------------------------------------------------------------
#spatial autocorrelation of ASTHMA


#neighborhood matrix
nb <- poly2nb(df$geometry)
head(nb)

ww <-  nb2listw(nb, style='W')
moran.plot(df$ASTHMA, ww)

#Moran's I
moran(df$ASTHMA, ww, n=length(ww$neighbours), S0=Szero(ww))
moran.mc(df$ASTHMA, ww, nsim=99999)
locm_bm <- localmoran(df$ASTHMA, ww)
names(locm_bm) <- c("Ii", "E.Ii", "Var.Ii", "Z.Ii", "Pr()")

summary(locm_bm)

df <- cbind(df, locm_bm)

localM <- tm_shape(st_as_sf(df)) + 
          tm_fill("Ii", 
          palette = "RdBu",
          style = "pretty") +
          tm_borders(alpha=.4)

tmap_save(localM,filename = "../figures/localM.png")


#scale the variable of interest and save it to a new column
df$ASTHMA_scaled <- scale(df$ASTHMA) %>% as.vector()
#create a spatial lag variable and save it to a new column
df$ASTHMA_lag <- lag.listw(ww, df$ASTHMA_scaled)
summary(df$ASTHMA)
summary(df$ASTHMA_lag)
x <- df$ASTHMA
y <- df$ASTHMA_lag
xx <- data_frame(x,y)
moran.plot(x, ww)
#dataframe$new_variable <- ifelse(dataframe$some_numeric_var < 100, "smaller than 100", "not smaller than 100")
df <- st_as_sf(df) %>% 
  mutate(quad_sig = ifelse(df$ASTHMA_scaled > 0 & 
                             df$ASTHMA_lag > 0 & 
                             locm_bm[,5] <= 0.05, 
                           "high-high",
                           ifelse(df$ASTHMA_scaled <= 0 & 
                                    df$ASTHMA_lag <= 0 & 
                                    locm_bm[,5] <= 0.05, 
                                  "low-low", 
                                  ifelse(df$ASTHMA_scaled > 0 & 
                                           df$ASTHMA_lag <= 0 & 
                                           locm_bm[,5] <= 0.05, 
                                         "high-low",
                                         ifelse(df$ASTHMA_scaled <= 0 & 
                                                  df$ASTHMA_lag > 0 & 
                                                  locm_bm[,5] <= 0.05,
                                                "low-high", 
                                                "non-significant")))))
table(df$quad_sig)
nrow(locm_bm[locm_bm[,5] <= 0.05,])
m <- qtm(df, fill="quad_sig", fill.title="ASTHMA LISA")
tmap_save(m,filename = "../figures/ASTHMA_LISA_v2.png")
st_write(df, "../data/localM.shp", append = FALSE)

#Getis-Ord
local_g <- localG(df$ASTHMA, ww)
local_g <- cbind(df, as.matrix(local_g))

names(local_g)[names(local_g)=="as.matrix.local_g."] <- "gstat"

gstar <-tm_shape(local_g) + 
          tm_fill("gstat", 
          palette = "RdBu",
          style = "pretty") +
          tm_borders(alpha=.4)

tmap_save(gstar,filename = "../figures/GSTAR.png")
st_write(local_g,"../data/localG.shp", append = FALSE)

#-----------------------------------------------------------------------------------------

#satscan
ast <- cbind(df$GEOID,df$ASTHMA) %>%
  as.data.frame()
cents <- st_coordinates(st_centroid(df)$geometry)

cents_geoid <- cbind(df$GEOID, cents) %>%
  as.data.frame()

sat <- left_join(ast, cents_geoid, by = "V1") %>%
  left_join(., Tot_pop, by = c("V1" = "GEOID"))

names(sat) <- c("GEOID", "ASTHMA", "X", "Y", "POP")

write.csv(sat, "../data/satscan.csv")


#-----------------------------------------------------------------------------------------
#variable maps
variable_map <- tm_shape(df) +
  tm_fill(c("RPL_THEME2", "RPL_THEME3", "RPL_THEME4", "RATIO", "MEDCENTER", "PM25"),
          style="fisher", 
          palette=list("Reds", "Oranges", "Greens", "Blues", "Purples", "Greys"),
          title=c("Theme 2", "Theme 3", "Theme 4", "Tree Cover", "Medical Centers", "PM2.5")) +
  tm_borders() +
  tm_layout(legend.position = c("left", "bottom")) 
tmap_save(variable_map,filename = "../figures/variable_map.png")

#-----------------------------------------------------------------------------------------

save.image("ws.RData")
