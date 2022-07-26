library(sf)
library(RColorBrewer)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set workspace dir


load("ws.RData")


pal <- brewer.pal(7, "OrRd") # we select 7 colors from the palette
class(pal)


jpeg(file="../figures/asthma.jpg")
plot(st_as_sf(df)["ASTHMA"],
     main = "Current Asthma Prevalence Among Adults Aged 18 and older (%)",
     breaks = "quantile", nbreaks = 7,
     pal = pal)
dev.off()

jpeg(file="../figures/RPL2.jpg")
plot(st_as_sf(df)["RPL_THEME2"],
     main = "RPL_THEME2",
     breaks = "quantile", nbreaks = 7,
     pal = pal)
dev.off()

jpeg(file="../figures/RPL3.jpg")
plot(st_as_sf(df)["RPL_THEME3"],
     main = "RPL_THEME3",
     breaks = "quantile", nbreaks = 7,
     pal = pal)
dev.off()

jpeg(file="../figures/RPL4.jpg")
plot(st_as_sf(df)["RPL_THEME4"],
     main = "RPL_THEME4",
     breaks = "quantile", nbreaks = 7,
     pal = pal)
dev.off()

jpeg(file="../figures/canopy.jpg")
plot(st_as_sf(df)["RATIO"],
     main = "Canopy Cover (%)",
     breaks = "quantile", nbreaks = 7,
     pal = pal)
dev.off()

jpeg(file="../figures/access.jpg")
plot(st_as_sf(df)["HOSP"],
     main = "Hospital Density",
     breaks = "quantile", nbreaks = 7,
     pal = pal)
dev.off()

