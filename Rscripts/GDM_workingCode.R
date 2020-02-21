##GDM small project---------------
#Simoes, MVP-----------------
setwd('C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new')
x_total<- read.csv('ck_total.csv')
#dim(x_total)
df2 <- read.csv('ck_deep.csv') #the one i want to exclude
ck_shallow <- x_total[!(x_total$scientific %in% df2$scientific),]
write.csv(ck_shallow, 'ck_shallow2.csv', , row.names = F)
df2 <- read.csv('ck_shallow.csv') #the one i want to exclude
ck_deep <- x_total[!(x_total$scientific %in% df2$scientific),]
write.csv(ck_deep, 'ck_deep2.csv', , row.names = F)
#dim(ck_deep)
write.csv(ck_deep, 'ck_deep2.csv', , row.names = F)
deep <-read.csv("ck_deep2.csv")
shallow <-read.csv("ck_shallow2.csv")
#dim(shallow)
deep$scientific

#Make folders with csv files of species with  >30 records ----------------------------------------------------
#Note> I created another folder to extract the species and occ that beneficial dataset had _
#The merged datasets of Isopods n>20 occ, is in Ku_ENM2
setwd('C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Deep')
x <- read.csv('ck_shallow2.csv') #head(x)
sp <- unique(x$scientific) #str(sp)

n = 20 #number of records 
cols = 1:3 # columns extracted from dataset - species, long, lat
for (i in 1:length(sp)) {
  sptable <- x[x$scientific == sp[i], cols]
  if (dim(sptable)[1] >= n) {
    dir.create(paste(sp[i], collapse = "_"))
    write.csv(sptable, paste(paste(sp[i], collapse = '_'), '\\',
                             paste(sp[i], collapse = '_'), ".csv", sep = ""), row.names = FALSE)
  }
  
}
#Crop rasters------------------------
library(sp)
library(rgeos)
library(raster)
library(rworldmap)
install.packages('sf')
library(sf)
htspt <- st_read("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\Hotspots.shp") #head
env1 <-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\present\\cveloc.asc")
env2<-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\present\\etoporst.asc")
env3<-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\present\\ice.asc")
env4<-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\present\\salinity.asc")
env5<-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\present\\temp.asc")
a <- mask(env1,htspt)
b <- mask(env2,htspt)
c <- mask(env3,htspt)
d <- mask(env4,htspt)
e <- mask(env5,htspt)
writeRaster(a, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\present\\cveloc.asc", format="ascii",overwrite=TRUE)
writeRaster(b, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\present\\etoporst.asc", format="ascii",overwrite=TRUE)
writeRaster(c, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\present\\ice.asc", format="ascii",overwrite=TRUE)
writeRaster(d, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\present\\salinity.asc", format="ascii",overwrite=TRUE)
writeRaster(e, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\present\\temp.asc", format="ascii",overwrite=TRUE)

#Future layers-------------------
#2050_rcp26
env1_f <-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\2050_rcp26\\cveloc.asc")
env2_f<-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\2050_rcp26\\etoporst.asc")
env3_f<-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\2050_rcp26\\ice.asc")
env4_f<-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\2050_rcp26\\salinity.asc")
env5_f<-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\2050_rcp26\\temp.asc")
a <- mask(env1_f,htspt)
b <- mask(env2_f,htspt)
c <- mask(env3_f,htspt)
d <- mask(env4_f,htspt)
e <- mask(env5_f,htspt)
writeRaster(a, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\2050_rcp26\\cveloc.asc", format="ascii",overwrite=TRUE)
writeRaster(b, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\2050_rcp26\\etoporst.asc", format="ascii",overwrite=TRUE)
writeRaster(c, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\2050_rcp26\\ice.asc", format="ascii",overwrite=TRUE)
writeRaster(d, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\2050_rcp26\\salinity.asc", format="ascii",overwrite=TRUE)
writeRaster(e, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\2050_rcp26\\temp.asc", format="ascii",overwrite=TRUE)

#2050_rcp85
env1_f <-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\2050_rcp85\\cveloc.asc")
env2_f<-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\2050_rcp85\\etoporst.asc")
env3_f<-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\2050_rcp85\\ice.asc")
env4_f<-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\2050_rcp85\\salinity.asc")
env5_f<-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\2050_rcp85\\temp.asc")
a <- mask(env1_f,htspt)
b <- mask(env2_f,htspt)
c <- mask(env3_f,htspt)
d <- mask(env4_f,htspt)
e <- mask(env5_f,htspt)
writeRaster(a, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\2050_rcp85\\cveloc.asc", format="ascii",overwrite=TRUE)
writeRaster(b, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\2050_rcp85\\etoporst.asc", format="ascii",overwrite=TRUE)
writeRaster(c, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\2050_rcp85\\ice.asc", format="ascii",overwrite=TRUE)
writeRaster(d, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\2050_rcp85\\salinity.asc", format="ascii",overwrite=TRUE)
writeRaster(e, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\2050_rcp85\\temp.asc", format="ascii",overwrite=TRUE)

#2100_rcp26
env1_f <-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\2100_rcp26\\cveloc.asc")
env2_f<-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\2100_rcp26\\etoporst.asc")
env3_f<-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\2100_rcp26\\ice.asc")
env4_f<-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\2100_rcp26\\salinity.asc")
env5_f<-raster("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\2100_rcp26\\temp.asc")
a <- mask(env1_f,htspt)
b <- mask(env2_f,htspt)
c <- mask(env3_f,htspt)
d <- mask(env4_f,htspt)
e <- mask(env5_f,htspt)
writeRaster(a, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\2100_rcp26\\cveloc.asc", format="ascii",overwrite=TRUE)
writeRaster(b, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\2100_rcp26\\etoporst.asc", format="ascii",overwrite=TRUE)
writeRaster(c, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\2100_rcp26\\ice.asc", format="ascii",overwrite=TRUE)
writeRaster(d, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\2100_rcp26\\salinity.asc", format="ascii",overwrite=TRUE)
writeRaster(e, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\2100_rcp26\\temp.asc", format="ascii",overwrite=TRUE)



#Playing with GDM---------------------------------------------
install.packages("gdm")
library(gdm)
library(raster)
setwd('C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow')

env.files <- list.files("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\present", full.names = TRUE)
env <- stack(env.files)
siteraster <- env[[1]]
siteraster[!is.na(siteraster)] <- seq(1:length(siteraster[!is.na(siteraster)]))
plot(siteraster)
spdata <- data.frame(scientific = "dummy",decimalLon = 0,decimalLat = 0)
for(i in list.files("~/Downloads/Marianna/Copepoda_deep/copepod_csv/", full.names = TRUE)){
  thisdata <- read.csv(i)
  spdata <- rbind(spdata, thisdata)
}
spdata <- read.csv('ck_shallow3.csv') #check order of long and lat!
head(spdata)
spdata <- spdata[-1,]
head(spdata$site)
spdata$site <- extract(siteraster, spdata[,2:3])
colnames(spdata) <- c("species", "long", "lat", "site")
spdata <- spdata[,c("species", "site", "lat", "long")]
spdata<-spdata[complete.cases(spdata), ] 
envtab <- cbind(extract(env, spdata[,4:3]), spdata[,2:4])
head(envtab)

gdmtab <- formatsitepair(bioData = spdata, bioFormat=2, XColumn="long", YColumn="lat",
                         sppColumn="species", siteColumn="site", predData=env)
str(gdmtab)
gdmtab<-gdmtab[complete.cases(gdmtab), ] 
head(gdm.all)

gdm.all <- gdm(gdmtab, geo=T)
str(gdm.all)
plot(gdm.all$geo)
length(gdm.all$predictors)
plot(gdm.all$predicted)
summary(gdm.all)

#future-------------------------
future.env.files <- list.files("C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Variables\\Cropped\\2050_rcp26", full.names = TRUE)
future.env <- stack(future.env.files)

future.env[is.na(env)] <- NA
env[is.na(future.env)] <- NA

timepred <- predict(gdm.all, env, time=T, predRasts=future.env)
plot(timepred)
writeRaster(timepred, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Predict_2050_rcp26.asc", format="ascii",overwrite=TRUE)


# transform climate rasters & plot pattern
rastTrans <- gdm.transform(gdm.all, env)
plot(rastTrans)
rastDat <- na.omit(getValues(rastTrans))
#rastDat <- sampleRandom(rastTrans, 50000) # can use if rasters are large
pcaSamp <- prcomp(rastDat)
# note the use of the 'index' argument
pcaRast <- predict(rastTrans, pcaSamp, index=1:3)
writeRaster(pcaRast, filename="C:\\Users\\admin\\Desktop\\GDM_Hotspots\\new\\ck_Shallow\\Predict_2050_rcp262.asc", format="ascii",overwrite=TRUE)

# scale rasters
pcaRast[[1]] <- (pcaRast[[1]]-pcaRast[[1]]@data@min) /
  (pcaRast[[1]]@data@max-pcaRast[[1]]@data@min)*255
pcaRast[[2]] <- (pcaRast[[2]]-pcaRast[[2]]@data@min) /
  (pcaRast[[2]]@data@max-pcaRast[[2]]@data@min)*255
pcaRast[[3]] <- (pcaRast[[3]]-pcaRast[[3]]@data@min) /
  (pcaRast[[3]]@data@max-pcaRast[[3]]@data@min)*255
plotRGB(pcaRast, r=1, g=2, b=3)