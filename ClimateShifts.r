library(rgdal)
library(raster)
lgm = stack(raster("C:/Users/mck14/Dropbox/bioclime/cclgmbi_10m/cclgmbi1.tif"),
            raster("C:/Users/mck14/Dropbox/bioclime/cclgmbi_10m/cclgmbi7.tif"),
            raster("C:/Users/mck14/Dropbox/bioclime/cclgmbi_10m/cclgmbi12.tif"))
plot(lgm)

mid = stack(raster("C:/Users/mck14/Dropbox/bioclime/ccmidbi_10m/ccmidbi1.tif"),
            raster("C:/Users/mck14/Dropbox/bioclime/ccmidbi_10m/ccmidbi7.tif"),
            raster("C:/Users/mck14/Dropbox/bioclime/ccmidbi_10m/ccmidbi12.tif"))
plot(mid)
cur = stack(raster("C:/Users/mck14/Dropbox/bioclime/bio_10m_bil/bio1.bil"),
        raster("C:/Users/mck14/Dropbox/bioclime/bio_10m_bil/bio7.bil"),
        raster("C:/Users/mck14/Dropbox/bioclime/bio_10m_bil/bio12.bil"))
plot(cur)

fut50 = stack(raster("C:/Users/mck14/Dropbox/bioclime/cc60bi50/cc60bi501.tif"),
            raster("C:/Users/mck14/Dropbox/bioclime/cc60bi50/cc60bi507.tif"),
            raster("C:/Users/mck14/Dropbox/bioclime/cc60bi50/cc60bi5012.tif"))
plot(fut50)

fut70 = stack(raster("C:/Users/mck14/Dropbox/bioclime/cc60bi70/cc60bi701.tif"),
              raster("C:/Users/mck14/Dropbox/bioclime/cc60bi70/cc60bi707.tif"),
              raster("C:/Users/mck14/Dropbox/bioclime/cc60bi70/cc60bi7012.tif"))
plot(fut70)

?extract
extract(fut70,c(-45,45))
head(as.matrix(raster("C:/Users/mck14/Dropbox/bioclime/cclgmbi_10m/cclgmbi1.tif")))
as.vector(extent(fut70))
all = stack(lgm,mid,cur,fut50,fut70)
cliMat = cbind(coordinates(all),as.matrix(all))
head(cliMat)
whr = apply(cliMat[,-c(1,2)],1,function(x) sum(!is.na(x)) != 0)
cliMat = cliMat[whr,]
colnames(cliMat)[1] = 'lon'
colnames(cliMat)[2] = 'lat'


whr = which.min(distance(c(-78.9072,35.9886),cliMat[,c(1,2)]))

curP = cliMat[whr,c('bio1','bio7','bio12')]
whrGM = which(rank(distance(curP,cliMat[,c('cclgmbi1','cclgmbi7', 'cclgmbi12')])) < 21)
whrMH = which(rank(distance(curP,cliMat[,c('ccmidbi1','ccmidbi7', 'ccmidbi12')]))< 41)
whr50 = which(rank(distance(curP,cliMat[,c('cc60bi501','cc60bi507', 'cc60bi5012')]))< 21)
whr70 = which(rank(distance(curP,cliMat[,c('cc60bi701','cc60bi707', 'cc60bi7012')]))< 21)

cliMat[c(whrGM,whrMH,whr,whr50,whr70),c(1,2)]

library(maps)
par(mfrow= c(1,1),mar=c(2,2,2,1))
plot(extent(all),xlab='',ylab='',bty='n',xaxt='n',yaxt='n')
map(,add=T,col=1)
map('state',add=T,col=1)
points(cliMat[c(whrGM),1],cliMat[c(whrGM),2],col='blue',pch=20,cex=.75)
points(cliMat[c(whrMH),1],cliMat[c(whrMH),2],col='green',pch=20,cex=.75)
points(cliMat[c(whr),1],cliMat[c(whr),2],col='grey',pch=20,cex=1.5)
points(cliMat[c(whr50),1],cliMat[c(whr50),2],col='orange',pch=20,cex=.75)
points(cliMat[c(whr70),1],cliMat[c(whr70),2],col='red',pch=20,cex=.75)






