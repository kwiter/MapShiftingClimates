###
#Summer 2015 MC Kwit
#Data: 
#  Bioclime:Last glacial Max, mid Holocene, current, 2050, and 2070
#           Community Climate System Model (CCSM4.0)
#           Representative Concentration Pathway RCP6.0   1.3 (0.8 to 1.8) 	2.2 (1.4 to 3.1)
#           Mean annual temp, Differance in max and min temps (Temp Seasonality)
#           Total PPT, Coefficient of Variation of PPT
#  Shapefiles: USA, World
#  Algorithm:
#   Select state
#     For each raster point in state
#       Calculate the euclidean distance between its current climate to given climate scenarios
#       Store the 50 top matches per point
#      For each poi(nt on raster count the number of times it matched the current state climate
#        Normalize (0,1)
#   Store:
#     allWHR.rds
#     cliCoor.rds
##


library(rgdal)
library(raster)
library(shapefiles)
library(plyr)
library(maps)
library(devtools)
library(rgeos)
github = "https://raw.githubusercontent.com/kwiter/mckFUNCTIONS/master/mckFUNCTIONS.r"
source_url(github)

path = "C:/Users/mck14/Dropbox"
path ="/home/kwit/Dropbox"
path ="/home/mckwit/Dropbox"

tmp = pathZip(paths(path,"/bioclime/cclgmbi_10m.zip"),c("cclgmbi1.tif","cclgmbi7.tif","cclgmbi12.tif","cclgmbi15.tif"))
lgm = stack(raster(tmp[1]),
            raster(tmp[2]),
            raster(tmp[3]),
            raster(tmp[4]))
#plot(lgm)

tmp = pathZip(paths(path,"/bioclime/ccmidbi_10m.zip"),
              c("ccmidbi1.tif","ccmidbi7.tif","ccmidbi12.tif","ccmidbi15.tif"))
mid = stack(raster(tmp[1]),
           raster(tmp[2]),
           raster(tmp[3]),
           raster(tmp[4]))
#plot(mid)

tmp = pathZip(paths(path,"/bioclime/bio_10m_bil.zip"),
              c("bio1.bil","bio7.bil","bio12.bil","bio15.bil"))
cur = stack(raster(tmp[1]),
            raster(tmp[2]),
            raster(tmp[3]),
            raster(tmp[4]))
#plot(cur)

tmp = pathZip(paths(path,"/bioclime/cc60bi50.zip"),
              c("cc60bi501.tif","cc60bi507.tif","cc60bi5012.tif","cc60bi5015.tif"))
fut50 = stack(raster(tmp[1]),
              raster(tmp[2]),
              raster(tmp[3]),
              raster(tmp[4]))
#plot(fut50)

tmp = pathZip(paths(path,"/bioclime/cc60bi70.zip"),
              c("cc60bi701.tif","cc60bi707.tif","cc60bi7012.tif","cc60bi7015.tif"))
fut70 = stack(raster(tmp[1]),
              raster(tmp[2]),
              raster(tmp[3]),
              raster(tmp[4]))
#plot(fut70)

all = stack(lgm,mid,cur,fut50,fut70)


ext = extent(all)
ext@xmax = -40
#ext@xmin = -135 
#ext@ymax = 52
ext@ymin = 10  #25
noAm = crop(all, ext)
plot(noAm)
all = noAm
cliMat = cbind(coordinates(all),as.matrix(all))
saveRDS(ext,file = paste(path,'/MapShiftingClimates/extent.rds',sep="") )
saveRDS(cliMat[,1:2],file = paste(path,'/MapShiftingClimates/cliCoor_noAm.rds',sep="") )


data.shape<-readOGR(dsn=paste(path,"/Shapes/cb_2014_us_state_500k",sep=""),layer="cb_2014_us_state_500k")
allWHR= data.frame(state = character(),whr = numeric(),abund=numeric(),type=character())
#allWHR <- readRDS(file = paste(path,'/MapShiftingClimates/allWHR.rds',sep="") ) #include if loop stops
for(k in 1:length(data.shape@data$NAME)){  #Change save name depending on input raster

nc = data.shape[data.shape@data$NAME ==data.shape@data$NAME[k],]
nc <- spTransform(nc,crs(all))
rr <- mask(all, nc)
#plot(rr)
NCcliMat = cbind(coordinates(rr),as.matrix(rr))#clipped raster

print(data.shape@data$NAME[k])

whr = apply(cliMat[,-c(1,2)],1,function(x) sum(!is.na(x)) != 0)
cliMat = cliMat[whr,]
colnames(cliMat)[1] = 'lon'
colnames(cliMat)[2] = 'lat'

whrNC = apply(NCcliMat[,-c(1,2)],1,function(x) sum(!is.na(x)) != 0)
NCcliMat = NCcliMat[whrNC,]

if(is.null(dim(NCcliMat))){
  names(NCcliMat)[1] = 'lon'
  names(NCcliMat)[2] = 'lat'
  curP = NCcliMat[c('bio1','bio7','bio12','bio15')]
  single = TRUE
}else{
  colnames(NCcliMat)[1] = 'lon'
  colnames(NCcliMat)[2] = 'lat'
  curP = NCcliMat[,c('bio1','bio7','bio12','bio15')]
  single = FALSE
}

#durham = c(-78.9072,35.9886)

#whr = which.min(distance(durham,cliMat[,c(1,2)])) #one point

#whr = apply(cbind(X,Y),1,function(x) which.min(distance(x,cliMat[,c(1,2)]))) #n by c(x,y) points
#whr = unique(unlist(whr))

#for single point
if(single){
#  whrGM = which(rank(distance(curP,cliMat[,c('cclgmbi1','cclgmbi7', 'cclgmbi12','cclgmbi15')])) < 101)
#  whrMH = which(rank(distance(curP,cliMat[,c('ccmidbi1','ccmidbi7', 'ccmidbi12','ccmidbi15')]))< 101)
#  whr50 = which(rank(distance(curP,cliMat[,c('cc60bi501','cc60bi507','cc60bi5012','cc60bi5015')]))< 101)
#  whr70 = which(rank(distance(curP,cliMat[,c('cc60bi701','cc60bi707','cc60bi7012','cc60bi7015')]))< 101)
  whrCU = which(rank(distance(curP,cliMat[,c('bio1','bio7','bio12','bio15')]))< 101)
  
#  curP = NCcliMat[c('cclgmbi1','cclgmbi7', 'cclgmbi12', 'cclgmbi15')]
#  whrGMI = which(rank(distance(curP,cliMat[,c('bio1','bio7','bio12','bio15')]))< 101)
#  curP = NCcliMat[c('ccmidbi1','ccmidbi7', 'ccmidbi12', 'ccmidbi15')]
#  whrMHI = which(rank(distance(curP,cliMat[,c('bio1','bio7','bio12','bio15')]))< 101)
  curP = NCcliMat[c('cc60bi501','cc60bi507', 'cc60bi5012','cc60bi5015')]
  whr50I = which(rank(distance(curP,cliMat[,c('bio1','bio7','bio12','bio15')]))< 101)
  curP = NCcliMat[c('cc60bi701','cc60bi707', 'cc60bi7012','cc60bi7015')]
  whr70I = which(rank(distance(curP,cliMat[,c('bio1','bio7','bio12','bio15')]))< 101)
}

if(single==FALSE){
 # whrGM = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('cclgmbi1','cclgmbi7', 'cclgmbi12','cclgmbi15')])) < 51),.progress = "text")
 # whrMH = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('ccmidbi1','ccmidbi7', 'ccmidbi12', 'ccmidbi15')]))< 51),.progress = "text")
 # whr50 = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('cc60bi501','cc60bi507', 'cc60bi5012','cc60bi5015')]))< 51),.progress = "text")
 # whr70 = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('cc60bi701','cc60bi707', 'cc60bi7012','cc60bi7015')]))< 51),.progress = "text")
  whrCU = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('bio1','bio7','bio12','bio15')]))< 51),.progress = "text")

#  curP = NCcliMat[,c('cclgmbi1','cclgmbi7', 'cclgmbi12','cclgmbi15')]
#  whrGMI = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('bio1','bio7','bio12','bio15')])) < 51),.progress = "text")
#  curP = NCcliMat[,c('ccmidbi1','ccmidbi7', 'ccmidbi12', 'ccmidbi15')]
#  whrMHI = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('bio1','bio7','bio12','bio15')]))< 51),.progress = "text")
  curP = NCcliMat[,c('cc60bi501','cc60bi507', 'cc60bi5012','cc60bi5015')]
  whr50I = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('bio1','bio7','bio12','bio15')]))< 51),.progress = "text")
  curP = NCcliMat[,c('cc60bi701','cc60bi707', 'cc60bi7012','cc60bi7015')]
  whr70I = alply(curP,1,function(x) which(rank(distance(x,cliMat[,c('bio1','bio7','bio12','bio15')]))< 51),.progress = "text")
}

#density of points
#whrGM = table(unlist(whrGM))
#whrMH = table(unlist(whrMH))
#whr50 = table(unlist(whr50))
#whr70 = table(unlist(whr70))
whrCU = table(unlist(whrCU))



#density of points
#whrGMI = table(unlist(whrGMI))
#whrMHI = table(unlist(whrMHI))
whr50I = table(unlist(whr50I))
whr70I = table(unlist(whr70I))

#typeAll = c(rep('whrGM',len(whrGM)),rep('whrMH',len(whrMH)),rep('whrCU',len(whrCU)),rep('whr50',len(whr50)),rep('whr70',len(whr70)),
#rep('whrGMI',len(whrGMI)),rep('whrMHI',len(whrMHI)),rep('whr50I',len(whr50I)),rep('whr70I',len(whr70I)))
typeAll = c(rep('whrCU',len(whrCU)),rep('whr50I',len(whr50I)),rep('whr70I',len(whr70I)))

#whrAll = c(whrGM,whrMH,whrCU,whr50,whr70,whrGMI,whrMHI,whr50I,whr70I)
whrAll = c(whrCU,whr50I,whr70I)
typeAll = typeAll[whrAll > 1]
whrAll = whrAll[whrAll > 1]
tmp = data.frame( state = rep(data.shape@data$NAME[k],len(whrAll)),whr = names(whrAll) ,abund = whrAll, type = typeAll)
allWHR = rbind(allWHR,tmp) 
saveRDS(allWHR,file = paste(path,'/MapShiftingClimates/allWHR_noAm.rds',sep="") )
print(data.shape@data$NAME[k])
}

#clean 3 duplicates from allWHR
#Nebraska
#New Hampshire
#New Mexico
whr = which(allWHR$state == sort(unique(allWHR$state))[28])
whr = whr[1:(len(whr)/2)]
allWHR = allWHR[-whr,]
whr = which(allWHR$state == sort(unique(allWHR$state))[30])
whr = whr[1:(len(whr)/2)]
allWHR = allWHR[-whr,]
whr = which(allWHR$state == sort(unique(allWHR$state))[32])
whr = whr[1:(len(whr)/2)]
allWHR = allWHR[-whr,]
saveRDS(allWHR,file = paste(path,'/MapShiftingClimates/allWHR.rds',sep="") )
#all points
#whrGM = unique(unlist(whrGM))
#whrMH = unique(unlist(whrMH))
#whr50 = unique(unlist(whr50))
#whr70 = unique(unlist(whr70))
#whrCU = unique(unlist(whrCU))

#cliMat[c(whrGM,whrMH,whr,whr50,whr70),c(1,2)]
allWHR = readRDS(file = paste(path,'/MapShiftingClimates/allWHR.rds',sep=""))
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

allWHR = readRDS(file = paste(path,'/MapShiftingClimates/allWHR.rds',sep=""))
whrST = which(allWHR$state == 'North Carolina' & (allWHR$type == 'whrCU' |allWHR$type == 'whr50I'|allWHR$type == 'whr70I'))
whr = a.n(a.c(allWHR[whrST,'whr']))
colors = c("#EDC951","#CC333F","#00A0B0","#6A4A3C","#EB6841")
colours = rep(NA,len(whrST))
for(i in 1:len(whrST)){
  colours[i] = adjustcolor(colors[a.n(as.factor(allWHR[whrST[i],'type']))], alpha.f = allWHR[whrST[i],'abund'])
  #if(colours[i] < .1){colours[i] = adjustcolor(colors[a.n(as.factor(allWHR[whrST[i],'type']))], alpha.f = 0)}
}

ext = extent(all)
ext@xmax = -40
ext@xmin = -135 
ext@ymax = 52
ext@ymin = 25  #25
par(mar=c(2,2,2,1), oma=c(0,0,0,0), bg="white", xpd=FALSE, xaxs="r", yaxs="i", mgp=c(2.1,.3,0), las=1, col.axis="#434343", col.main="#343434", tck=0, lend=1)
plot(ext, type="n", bty="n", las=1,
     xlab='', ylab='', family="Helvetica", cex.main=1.5, cex.axis=0.8, cex.lab=0.8, xaxt="n", yaxt="n")
map(,add=T,col=1)
map('state',add=T,col='black')
points(cliCoor[whr,1],cliCoor[whr,2],col=colours,pch=20,cex=.75)
tapply(palette()[a.n(as.factor(allWHR[whrST,'type']))],allWHR[whrST,'type'],unique)
tapply(palette()[a.n(as.factor(allWHR[whrST,'type']))],a.n(as.factor(allWHR[whrST,'type'])),unique)
legend(-62,34,c("Current","2050","2070"),fill=c("#00A0B0","#EDC951","#CC333F"),border='white',bty='n') #us ledgend
title("The climate of North Carolina from today to 2070",line=0.2)
title( "or hello Texarkana",line=-.8)
#legend(47,-30,c("Current","2050","2070"),fill=c("#00A0B0","#EDC951","#CC333F"))



###
specnames = c("acerrubr","lirituli","queralba","querrubr")
s.name = c("ACRU","LITU","QUAL","QURU")
f.name = c("red maple","tulip poplar","white oak","red oak")
dims = sqrt(length(s.name))
par(mfrow = c(floor(dims),ceiling(dims)),
    mar=c(0,0,0,0),xpd=F,pty='s')
for(i in 1:length(specnames)){
  map <- readOGR(paste("little/",s.name[i],sep=""),specnames[i],verbose=F)
  map('world',xlim=c(-98,-65),ylim=c(25,52),col='#E8D9C5',fill=T,mar=c(0,0,0,0), bg="lightblue",border='grey',myborder=0)
  plot(map,add=T,col=rgb(.1,.4,.2,.4),border=NA)
  map('state',xlim=c(-98,-65),ylim=c(25,52),add=T,col='grey')
  points(cliCoor[whr,1],cliCoor[whr,2],col=colours,pch=20,cex=.6)
}  

###

points(cliMat[c(a.n(names(whrGM))),1],cliMat[c(a.n(names(whrGM))),2],col=rgb(.6,.6,.6,whrGM/max(whrGM)),pch=20,cex=.4)
points(cliMat[c(a.n(names(whrMH))),1],cliMat[c(a.n(names(whrMH))),2],col=rgb(.6,.6,.6,whrMH/max(whrMH)),pch=20,cex=.4)
points(cliMat[c(a.n(names(whrCU))),1],cliMat[c(a.n(names(whrCU))),2],col=groups,pch=20,cex=.4)
points(cliMat[c(a.n(names(whr50))),1],cliMat[c(a.n(names(whr50))),2],col=rgb(.2,.2,.6,whr50/max(whr50)),pch=20,cex=.4)
points(cliMat[c(a.n(names(whr70))),1],cliMat[c(a.n(names(whr70))),2],col=rgb(.2,.2,.6,whr70/max(whr70)),pch=20,cex=.4)

points(cliMat[c(a.n(names(whrGMI))),1],cliMat[c(a.n(names(whrGMI))),2],col=rgb(.6,.6,.6,whrGMI/max(whrGMI)),pch=20,cex=.4)
points(cliMat[c(a.n(names(whrMHI))),1],cliMat[c(a.n(names(whrMHI))),2],col=rgb(.6,.6,.6,whrMHI/max(whrMHI)),pch=20,cex=.4)
points(cliMat[c(a.n(names(whrCU))),1],cliMat[c(a.n(names(whrCU))),2],col=rgb(.6,.6,.6,whrCU/max(whrCU)),pch=20,cex=.4)
points(cliMat[c(a.n(names(whr50I))),1],cliMat[c(a.n(names(whr50I))),2],col=rgb(.2,.2,.6,whr50I/max(whr50I)),pch=20,cex=.4)
points(cliMat[c(a.n(names(whr70I))),1],cliMat[c(a.n(names(whr70I))),2],col=rgb(.2,.2,.6,whr70I/max(whr70I)),pch=20,cex=.4)


mydata = cliMat[a.n(names(whr70I)),c(1,2)]
d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="single")
ng = 5
groups <- cutree(fit, k=ng) # cut tree into k clusters

par(mfrow= c(1,1),mar=c(2,2,2,1))
plot(extent(all),xlab='',ylab='',bty='n',xaxt='n',yaxt='n')
map(,add=T,col=1)
map('state',add=T,col=1)
for(j in 1:ng){
if(length(a.n(whr70I[groups==j]))<10)next
tmp = spline.poly(spatialQuantAbund(mydata[groups==j,],a.n(whr70I[groups==j]),nbrks = 10,quants = c(.2,.8)),50)
polygon(tmp,lwd=3,lty=2,border=rgb(.4,.4,.8,1),col=rgb(.4,.4,.8,.5))
}

spatialQuantAbund = function(xyMat,abund,nbrks = 10,quants = c(.025,.975)){
  stdX = sd(xyMat[,1],na.rm=T)
  stdY = sd(xyMat[,2],na.rm=T)
  
  xyMat[,1] = (xyMat[,1]) /stdX
  xyMat[,2] = (xyMat[,2]) /stdY
  
  mx = median(xyMat[,1])
  my = median(xyMat[,2])
  
  tmp = angleTo(c(mx,my),xyMat)
  angs = tmp$ang
  tx = tmp$tx
  ty = tmp$ty
  quad = rep(0,len(tx))
  quad[tx ==  1 & ty ==  1] = 1
  quad[tx ==  1 & ty == -1] = 1
  quad[tx == -1 & ty == -1] = -1
  quad[tx == -1 & ty ==  1] = -1
  
  dis = distance(c(mx,my),xyMat)*abund
  
  
  breaks <- seq(-90,90,length = nbrks)
  bins <- cut(angs,breaks = breaks)
  bins[which(is.na(bins))] = sort(unique(bins))[1]
  iBin = sort(unique(cut(seq(-89,89,by = 1),breaks = breaks)))
  tmp = matrix(NA,len(iBin),3)
  tmp[,1] = breaks[-1] - diff(breaks)[1]/2
  for(i in 1:len(iBin)){
    whr = which(bins == iBin[i])
    if(length(whr) == 0){tmp[i,c(2,3)] = c(0,0)
    }else{
    tmp[i,c(2,3)] = quantile(dis[whr] * quad[whr] / abund[whr],quants)
    }
  }
  xs = (cos(c(tmp[,1],tmp[,1])*2*pi/360)*c(tmp[,2],tmp[,3]) + mx)*stdX 
  ys = (sin(c(tmp[,1],tmp[,1])*2*pi/360)*c(tmp[,2],tmp[,3]) + my)*stdY 
  cbind(xs,ys)
}
