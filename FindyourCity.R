library(curl)
library(jsonlite)



library(rgdal)
library(raster)
library(shapefiles)
library(plyr)
library(maps)
library(devtools)
library(rgeos)
github = "https://raw.githubusercontent.com/kwiter/mckFUNCTIONS/master/mckFUNCTIONS.r"
source_url(github)

northAm = F


distance = function(valM,M){
  dims = length(valM)
  resp = rep(0,dim(M)[1])
  for(i in 1:dims){
    resp = resp + (M[,i]-valM[i])^2
  }
  sqrt(resp) 
}

# Convert degrees to radians
deg2rad <- function(deg) return(deg*pi/180)
# Calculates the geodesic distance between two points specified by radian latitude/longitude using
# Vincenty inverse formula for ellipsoids (vif)
gcd.vif <- function(long1, lat1, long2, lat2) {
  
  # WGS-84 ellipsoid parameters
  a <- 6378137         # length of major axis of the ellipsoid (radius at equator)
  b <- 6356752.314245  # ength of minor axis of the ellipsoid (radius at the poles)
  f <- 1/298.257223563 # flattening of the ellipsoid
  
  L <- long2-long1 # difference in longitude
  U1 <- atan((1-f) * tan(lat1)) # reduced latitude
  U2 <- atan((1-f) * tan(lat2)) # reduced latitude
  sinU1 <- sin(U1)
  cosU1 <- cos(U1)
  sinU2 <- sin(U2)
  cosU2 <- cos(U2)
  
  cosSqAlpha <- NULL
  sinSigma <- NULL
  cosSigma <- NULL
  cos2SigmaM <- NULL
  sigma <- NULL
  
  lambda <- L
  lambdaP <- 0
  iterLimit <- 100
  while (abs(lambda-lambdaP) > 1e-12 & iterLimit>0) {
    sinLambda <- sin(lambda)
    cosLambda <- cos(lambda)
    sinSigma <- sqrt( (cosU2*sinLambda) * (cosU2*sinLambda) +
                        (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda) )
    if (sinSigma==0) return(0)  # Co-incident points
    cosSigma <- sinU1*sinU2 + cosU1*cosU2*cosLambda
    sigma <- atan2(sinSigma, cosSigma)
    sinAlpha <- cosU1 * cosU2 * sinLambda / sinSigma
    cosSqAlpha <- 1 - sinAlpha*sinAlpha
    cos2SigmaM <- cosSigma - 2*sinU1*sinU2/cosSqAlpha
    if (is.na(cos2SigmaM)) cos2SigmaM <- 0  # Equatorial line: cosSqAlpha=0
    C <- f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha))
    lambdaP <- lambda
    lambda <- L + (1-C) * f * sinAlpha *
      (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)))
    iterLimit <- iterLimit - 1
  }
  if (iterLimit==0) return(NA)  # formula failed to converge
  uSq <- cosSqAlpha * (a*a - b*b) / (b*b)
  A <- 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
  B <- uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)))
  deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM^2) -
                                             B/6*cos2SigmaM*(-3+4*sinSigma^2)*(-3+4*cos2SigmaM^2)))
  s <- b*A*(sigma-deltaSigma) / 1000
  
  return(s) # Distance in km
}

# Calculates the geodesic distance between two points specified by degrees (DD) latitude/longitude using
# Haversine formula (hf), Spherical Law of Cosines (slc) and Vincenty inverse formula for ellipsoids (vif)
gcd <- function(long1, lat1, long2, lat2) {
  
  # Convert degrees to radians
  long1 <- deg2rad(long1)
  lat1 <- deg2rad(lat1)
  long2 <- deg2rad(long2)
  lat2 <- deg2rad(lat2)
  
  gcd.vif(long1, lat1, long2, lat2) 
}


path = "C:/Users/mck14/Dropbox"
path ="/home/kwit/Dropbox"
path ="/home/mckwit/Dropbox"

lgm = stack(raster(paste(path,"/bioclime/cclgmbi_10m/cclgmbi1.tif",sep='')),
            raster(paste(path,"/bioclime/cclgmbi_10m/cclgmbi7.tif",sep='')),
            raster(paste(path,"/bioclime/cclgmbi_10m/cclgmbi12.tif",sep='')),
            raster(paste(path,"/bioclime/cclgmbi_10m/cclgmbi15.tif",sep='')))
#plot(lgm)

mid = stack(raster(paste(path,"/bioclime/ccmidbi_10m/ccmidbi1.tif",sep='')),
            raster(paste(path,"/bioclime/ccmidbi_10m/ccmidbi7.tif",sep='')),
            raster(paste(path,"/bioclime/ccmidbi_10m/ccmidbi12.tif",sep='')),
            raster(paste(path,"/bioclime/ccmidbi_10m/ccmidbi15.tif",sep='')))
#plot(mid)
cur = stack(raster(paste(path,"/bioclime/bio_10m_bil/bio1.bil",sep='')),
            raster(paste(path,"/bioclime/bio_10m_bil/bio7.bil",sep='')),
            raster(paste(path,"/bioclime/bio_10m_bil/bio12.bil",sep='')),
            raster(paste(path,"/bioclime/bio_10m_bil/bio15.bil",sep='')))
#plot(cur)

fut50 = stack(raster(paste(path,"/bioclime/cc60bi50/cc60bi501.tif",sep='')),
              raster(paste(path,"/bioclime/cc60bi50/cc60bi507.tif",sep='')),
              raster(paste(path,"/bioclime/cc60bi50/cc60bi5012.tif",sep='')),
              raster(paste(path,"/bioclime/cc60bi50/cc60bi5015.tif",sep='')))
#plot(fut50)

fut70 = stack(raster(paste(path,"/bioclime/cc60bi70/cc60bi701.tif",sep='')),
              raster(paste(path,"/bioclime/cc60bi70/cc60bi707.tif",sep='')),
              raster(paste(path,"/bioclime/cc60bi70/cc60bi7012.tif",sep='')),
              raster(paste(path,"/bioclime/cc60bi70/cc60bi7015.tif",sep='')))
#plot(fut70)

all = stack(lgm,mid,cur,fut50,fut70)

ext = extent(all)
if(northAm == T) {
  ext@xmax = -40
  #ext@xmin = -135 
  #ext@ymax = 52
  ext@ymin = 0  #25
  all = crop(all, ext)
}

cliMat = cbind(coordinates(all),as.matrix(all))
colnames(cliMat)[1] = 'lon'
colnames(cliMat)[2] = 'lat'

international = F
city ='durham'
country = 'USA'
state = 'nc'
if(international == T) url = paste('http://api.wunderground.com/api/5c8ea9ca81cd5e18/geolookup/q/',country,'/',city,'.json',sep='')
if(international == F) url = paste('http://api.wunderground.com/api/5c8ea9ca81cd5e18/geolookup/q/',state,'/',city,'.json',sep='')
wu = readLines(url)
wu = fromJSON(wu)
lon = as.numeric(wu$location$lon)
lat = as.numeric(wu$location$lat)
LL = c(lon,lat)
whr = which.min(distance(LL,cliMat[,c(1,2)])) #one point

NCcliMat = cliMat[whr,]
curP = NCcliMat[c('bio1','bio7','bio12','bio15')]
whrCU = which(rank(distance(curP,cliMat[,c('bio1','bio7','bio12','bio15')]))< 11)
curP = NCcliMat[c('cc60bi501','cc60bi507', 'cc60bi5012','cc60bi5015')]
whr50I = which(rank(distance(curP,cliMat[,c('bio1','bio7','bio12','bio15')]))< 11)
curP = NCcliMat[c('cc60bi701','cc60bi707', 'cc60bi7012','cc60bi7015')]
whr70I = which(rank(distance(curP,cliMat[,c('bio1','bio7','bio12','bio15')]))< 11)

#whrCU = whrCU[which.min(gcd(LL[1], LL[2], cliMat[whrCU,1], cliMat[whrCU,2]))]
#whr50I = whr50I[which.min(gcd(LL[1], LL[2], cliMat[whr50I,1], cliMat[whr50I,2]))]
#whr70I = whr70I[which.min(gcd(LL[1], LL[2], cliMat[whr70I,1], cliMat[whr70I,2]))]

par(mar=c(2,2,2,1), oma=c(0,0,0,0), bg="white", xpd=FALSE, xaxs="r", yaxs="i", mgp=c(2.1,.3,0), las=1, col.axis="#434343", col.main="#343434", tck=0, lend=1)
plot(ext, type="n", bty="n", las=1,
     xlab='', ylab='', family="Helvetica", cex.main=1.5, cex.axis=0.8, cex.lab=0.8, xaxt="n", yaxt="n")
map(,add=T,col=1)
map('state',add=T,col='black')
points(cliMat[whrCU,1],cliMat[whrCU,2],  bg =rgb(.2,.8,.2,.6),   col=rgb(.2,.8,.2,1),pch=21,cex=1)
points(cliMat[whrCU,1],cliMat[whrCU,2], col='black',pch=20,cex=.1)
points(cliMat[whr50I,1],cliMat[whr50I,2],bg =rgb(.99,.84,.36,.6),col=rgb(.99,.84,.36,1),pch=20,cex=1.5)
points(cliMat[whr50I,1],cliMat[whr50I,2], col='black',pch=20,cex=.1)
points(cliMat[whr70I,1],cliMat[whr70I,2],bg =rgb(.8,.2,.2,.6),   col=rgb(.8,.2,.2,1),pch=20,cex=1.5)
points(cliMat[whr70I,1],cliMat[whr70I,2], col='black',pch=20,cex=.1)
#legend(-62,34,c("Current","2050","2070"),fill=c("#00A0B0","#EDC951","#CC333F"),border='white',bty='n') #us ledgend
#title(paste("The climate of", sts, "from today to 2070",line=0))
map('state',add=T,col=rgb(.8,.8,.8,.5))

url50 = paste('http://api.wunderground.com/api/5c8ea9ca81cd5e18/geolookup/q/',cliMat[whr50I[1],2],',',cliMat[whr50I[1],1],'.json',sep='')
url70 = paste('http://api.wunderground.com/api/5c8ea9ca81cd5e18/geolookup/q/',cliMat[whr70I[1],2],',',cliMat[whr70I[1],1],'.json',sep='')

url70con = paste('http://api.wunderground.com/api/5c8ea9ca81cd5e18/conditions/q/',cliMat[whr70I,2],',',cliMat[whr70I,1],'.json',sep='')
url70con = paste('http://api.wunderground.com/api/5c8ea9ca81cd5e18/hourly/q/',cliMat[whr70I,2],',',cliMat[whr70I,1],'.json',sep='')
#con = readLines(url70con)
#con = fromJSON(con)

conditions = c(wu$location$lon,wu$location$lat,wu$location$city,wu$location$state,wu$location$country_name,wu$location$wuiurl)

wu = readLines(url50)
wu = fromJSON(wu)
conditions = rbind(conditions,c(wu$location$lon,wu$location$lat,wu$location$city,wu$location$state,wu$location$country_name,wu$location$wuiurl))

wu = readLines(url70)
wu = fromJSON(wu)
conditions = rbind(conditions,c(wu$location$lon,wu$location$lat,wu$location$city,wu$location$state,wu$location$country_name,wu$location$wuiurl))
conditions

hour(Sys.time())


plot(1,1,type='n', xlim = range(cliMat[,'bio12'],na.rm=T),range(cliMat[,'bio1'],na.rm=T),col=rgb(.8,.8,.8,.1),pch=20,cex=.75)
whr = sample(1:dim(cliMat)[1],30000)
points(cliMat[whr,'bio12'],cliMat[whr,'bio1'],col=rgb(.8,.8,.8,.5),pch=20,cex=.75)
points(cliMat[whrCU,'bio12'],cliMat[whrCU,'bio1'],  bg =rgb(.2,.8,.2,.6),   col=rgb(.2,.8,.2,1),pch=21,cex=1)
points(cliMat[whrCU,'bio12'],cliMat[whrCU,'bio1'], col='black',pch=20,cex=.1)
points(cliMat[whr50I,'cc60bi5012'],cliMat[whr50I,'cc60bi501'],bg =rgb(.99,.84,.36,.6),col=rgb(.99,.84,.36,1),pch=20,cex=1.5)
points(cliMat[whr50I,'cc60bi5012'],cliMat[whr50I,'cc60bi501'], col='black',pch=20,cex=.1)
points(cliMat[whr70I,'cc60bi7012'],cliMat[whr70I,'cc60bi701'],bg =rgb(.8,.2,.2,.6),   col=rgb(.8,.2,.2,1),pch=20,cex=1.5)
points(cliMat[whr70I,'cc60bi7012'],cliMat[whr70I,'cc60bi701'], col='black',pch=20,cex=.1)

