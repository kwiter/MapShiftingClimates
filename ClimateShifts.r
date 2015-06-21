library(rgdal)
library(raster)
library(shapefiles)
library(maps)

path = "C:/Users/mck14/Dropbox"
path ="/home/kwit/Dropbox"

lgm = stack(raster(paste(path,"/bioclime/cclgmbi_10m/cclgmbi1.tif",sep='')),
            raster(paste(path,"/bioclime/cclgmbi_10m/cclgmbi7.tif",sep='')),
            raster(paste(path,"/bioclime/cclgmbi_10m/cclgmbi12.tif",sep='')))
plot(lgm)

mid = stack(raster(paste(path,"/bioclime/ccmidbi_10m/ccmidbi1.tif",sep='')),
            raster(paste(path,"/bioclime/ccmidbi_10m/ccmidbi7.tif",sep='')),
            raster(paste(path,"/bioclime/ccmidbi_10m/ccmidbi12.tif",sep='')))
plot(mid)
cur = stack(raster(paste(path,"/bioclime/bio_10m_bil/bio1.bil",sep='')),
        raster(paste(path,"/bioclime/bio_10m_bil/bio7.bil",sep='')),
        raster(paste(path,"/bioclime/bio_10m_bil/bio12.bil",sep='')))
plot(cur)

fut50 = stack(raster(paste(path,"/bioclime/cc60bi50/cc60bi501.tif",sep='')),
            raster(paste(path,"/bioclime/cc60bi50/cc60bi507.tif",sep='')),
            raster(paste(path,"/bioclime/cc60bi50/cc60bi5012.tif",sep='')))
plot(fut50)

fut70 = stack(raster(paste(path,"/bioclime/cc60bi70/cc60bi701.tif",sep='')),
              raster(paste(path,"/bioclime/cc60bi70/cc60bi707.tif",sep='')),
              raster(paste(path,"/bioclime/cc60bi70/cc60bi7012.tif",sep='')))
plot(fut70)

all = stack(lgm,mid,cur,fut50,fut70)
cliMat = cbind(coordinates(all),as.matrix(all))
NCcliMat = cbind(coordinates(rr),as.matrix(rr))#clipped raster

data.shape<-readOGR(dsn="/home/kwit/Dropbox/Shapes/cb_2014_us_state_500k",layer="cb_2014_us_state_500k")
nc = data.shape[data.shape@data$NAME =='North Carolina',]
nc <- spTransform(nc,crs(all))
rr <- mask(all, nc)
plot(rr)


whr = apply(cliMat[,-c(1,2)],1,function(x) sum(!is.na(x)) != 0)
cliMat = cliMat[whr,]
colnames(cliMat)[1] = 'lon'
colnames(cliMat)[2] = 'lat'

whrNC = apply(NCcliMat[,-c(1,2)],1,function(x) sum(!is.na(x)) != 0)
NCcliMat = NCcliMat[whrNC,]
colnames(NCcliMat)[1] = 'lon'
colnames(NCcliMat)[2] = 'lat'

durham = c(-78.9072,35.9886)
nc = map('state',"florida")
nc = map(,'china',fill=T)
X = nc$x
Y = nc$y

whr = which.min(distance(durham,cliMat[,c(1,2)])) #one point

whr = apply(cbind(X,Y),1,function(x) which.min(distance(x,cliMat[,c(1,2)]))) #n by c(x,y) points
whr = unique(unlist(whr))

curP = NCcliMat[,c('bio1','bio7','bio12')]

#for single point
whrGM = which(rank(distance(curP,cliMat[,c('cclgmbi1','cclgmbi7', 'cclgmbi12')])) < 101)
whrMH = which(rank(distance(curP,cliMat[,c('ccmidbi1','ccmidbi7', 'ccmidbi12')]))< 101)
whr50 = which(rank(distance(curP,cliMat[,c('cc60bi501','cc60bi507', 'cc60bi5012')]))< 101)
whr70 = which(rank(distance(curP,cliMat[,c('cc60bi701','cc60bi707', 'cc60bi7012')]))< 101)
whrCU = which(rank(distance(curP,cliMat[,c('bio1','bio7','bio12')]))< 101)


whrGM = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('cclgmbi1','cclgmbi7', 'cclgmbi12')])) < 51),.progress = "text")
whrMH = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('ccmidbi1','ccmidbi7', 'ccmidbi12')]))< 51),.progress = "text")
whr50 = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('cc60bi501','cc60bi507', 'cc60bi5012')]))< 51),.progress = "text")
whr70 = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('cc60bi701','cc60bi707', 'cc60bi7012')]))< 51),.progress = "text")
whrCU = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('bio1','bio7','bio12')]))< 51),.progress = "text")

#density of points
whrGM = table(unlist(whrGM))
whrMH = table(unlist(whrMH))
whr50 = table(unlist(whr50))
whr70 = table(unlist(whr70))
whrCU = table(unlist(whrCU))

curP = NCcliMat[,c('cclgmbi1','cclgmbi7', 'cclgmbi12')]
whrGMI = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('bio1','bio7','bio12')])) < 51),.progress = "text")
curP = NCcliMat[,c('ccmidbi1','ccmidbi7', 'ccmidbi12')])]
whrMHI = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('bio1','bio7','bio12')]))< 51),.progress = "text")
curP = NCcliMat[,c('cc60bi501','cc60bi507', 'cc60bi5012')]
whr50I = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('bio1','bio7','bio12')]))< 51),.progress = "text")
curP = NCcliMat[,c('cc60bi701','cc60bi707', 'cc60bi7012')]
whr70I = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('bio1','bio7','bio12')]))< 51),.progress = "text")

#density of points
whrGMI = table(unlist(whrGMI))
whrMHI = table(unlist(whrMHI))
whr50I = table(unlist(whr50I))
whr70I = table(unlist(whr70I))
whrCUI = table(unlist(whrCUI))

#all points
#whrGM = unique(unlist(whrGM))
#whrMH = unique(unlist(whrMH))
#whr50 = unique(unlist(whr50))
#whr70 = unique(unlist(whr70))
#whrCU = unique(unlist(whrCU))

#cliMat[c(whrGM,whrMH,whr,whr50,whr70),c(1,2)]

par(mfrow= c(1,1),mar=c(2,2,2,1))
plot(extent(all),xlab='',ylab='',bty='n',xaxt='n',yaxt='n')
map(,add=T,col=1)
map('state',add=T,col=1)
points(cliMat[c(whrGM),1],cliMat[c(whrGM),2],col='blue',pch=20,cex=.75)
points(cliMat[c(whrMH),1],cliMat[c(whrMH),2],col='green',pch=20,cex=.75)
points(cliMat[c(whr),1],cliMat[c(whr),2],col='grey',pch=20,cex=1.5)
points(cliMat[c(whr50),1],cliMat[c(whr50),2],col='orange',pch=20,cex=.75)
points(cliMat[c(whr70),1],cliMat[c(whr70),2],col='red',pch=20,cex=.75)
points(cliMat[c(whrCU),1],cliMat[c(whrCU),2],col='purple',pch=20,cex=.75)

points(NCcliMat[,1],NCcliMat[,2],cex=.1)

points(cliMat[c(a.n(names(whrGM))),1],cliMat[c(a.n(names(whrGM))),2],col=rgb(.2,.6,.6,whrGM/max(whrGM)),pch=20,cex=.4)
points(cliMat[c(a.n(names(whrMH))),1],cliMat[c(a.n(names(whrMH))),2],col=rgb(.6,.6,.6,whrMH/max(whrMH)),pch=20,cex=.4)
points(cliMat[c(a.n(names(whrCU))),1],cliMat[c(a.n(names(whrCU))),2],col=rgb(.2,.2,.6,whrCU/max(whrCU)),pch=20,cex=.4)
points(cliMat[c(a.n(names(whr50))),1],cliMat[c(a.n(names(whr50))),2],col=rgb(.2,.2,.6,whr50/max(whr50)),pch=20,cex=.4)
points(cliMat[c(a.n(names(whr70))),1],cliMat[c(a.n(names(whr70))),2],col=rgb(.2,.2,.6,whr70/max(whr70)),pch=20,cex=.4)



