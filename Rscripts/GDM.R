#make folders with csv files ----------------------------------------------------
#Note: I created another folder to extract the species and occ that beneficial dataset had _
#The merged datasets of Isopods n>20 occ, is in Ku_ENM2
setwd('C:\\Users\\admin\\Desktop\\GDM_Species\\Peracarida_shallow') 
x <- read.csv('C:\\Users\\admin\\Desktop\\Crustacea_OBIS\\Peracarida_shallow.csv')
sp <- unique(x$scientific) #str(sp)

n = 20 #number of records 
cols = 2:4 # columns extracted from dataset - species, long, lat
for (i in 1:length(sp)) {
  sptable <- x[x$scientific == sp[i], cols]
  if (dim(sptable)[1] >= n) {
      write.csv(sptable, paste(paste(sp[i], collapse = '_'), ".csv", sep = ""), row.names = FALSE)
      }
  
}



##GDM analysis------------------------------


install.packages("gdm")
library(gdm)
library(raster)

env.files <- list.files("~/Downloads/Marianna/Variables_deep/present/", full.names = TRUE)
env <- stack(env.files)
siteraster <- env[[1]]
siteraster[!is.na(siteraster)] <- seq(1:length(siteraster[!is.na(siteraster)]))

spdata <- data.frame(scientific = "dummy",decimalLon = 0,decimalLat = 0)
for(i in list.files("~/Downloads/Marianna/Copepoda_deep/copepod_csv/", full.names = TRUE)){
  thisdata <- read.csv(i)
  spdata <- rbind(spdata, thisdata)
}
spdata <- spdata[-1,]
head(spdata)
spdata$site <- extract(siteraster, spdata[,2:3])
colnames(spdata) <- c("species", "long", "lat", "site")
spdata <- spdata[,c("species", "site", "lat", "long")]
envtab <- cbind(extract(env, spdata[,4:3]), spdata[,2:4])
head(envtab)



gdmtab <- formatsitepair(bioData = spdata, bioFormat=2, XColumn="long", YColumn="lat",
                         sppColumn="species", siteColumn="site", predData=env)


gdm.copepods <- gdm(gdmtab, geo=T)
gdm.copepods
length(gdm.copepods$predictors)
plot(gdm.copepods, plot.layout = c(2,3))


future.env.files <- list.files("~/Downloads/Marianna/Variables_deep/2050_rcp26/", full.names = TRUE)
future.env <- stack(future.env.files)

future.env[is.na(env)] <- NA
env[is.na(future.env)] <- NA

timepred <- predict(gdm.copepods, env, time=T, predRasts=future.env)
plot(timepred)


# transform climate rasters & plot pattern
rastTrans <- gdm.transform(gdm.copepods, env)

rastDat <- na.omit(getValues(rastTrans))
#rastDat <- sampleRandom(rastTrans, 50000) # can use if rasters are large
pcaSamp <- prcomp(rastDat)
# note the use of the 'index' argument
pcaRast <- predict(rastTrans, pcaSamp, index=1:3)
# scale rasters
pcaRast[[1]] <- (pcaRast[[1]]-pcaRast[[1]]@data@min) /
  (pcaRast[[1]]@data@max-pcaRast[[1]]@data@min)*255
pcaRast[[2]] <- (pcaRast[[2]]-pcaRast[[2]]@data@min) /
  (pcaRast[[2]]@data@max-pcaRast[[2]]@data@min)*255
pcaRast[[3]] <- (pcaRast[[3]]-pcaRast[[3]]@data@min) /
  (pcaRast[[3]]@data@max-pcaRast[[3]]@data@min)*255
plotRGB(pcaRast, r=1, g=2, b=3)
