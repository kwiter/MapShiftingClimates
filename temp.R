x1 = rnorm(1000,20,5)
y1 = rnorm(1000,20,5)

x2 = rnorm(1000,30,5)
y2 = rnorm(1000,30,5)

plot(0,0,type='n',xlim = c(0,50),ylim = c(0,50))
points(x1,y1)
points(x2,y2,pch=20)

library(spatstat)

 

x1c= cut(x1,breaks = 0:50,include.lowest = TRUE,labels=seq(.5,49.5,by=1))
y1c= cut(y1,breaks = 0:50,include.lowest = TRUE,labels=seq(.5,49.5,by=1))
x2c= cut(x2,breaks = 0:50,include.lowest = TRUE,labels=seq(.5,49.5,by=1))
y2c= cut(y2,breaks = 0:50,include.lowest = TRUE,labels=seq(.5,49.5,by=1))
one = paste(x1c,y1c,sep="|") #
one = one[grep('NA',one,invert=T)]
two = paste(x2c,y2c,sep="|")
two = two[grep('NA',two,invert=T)]

one = table(one)
oneCoor = TtoC(names(one),symbol = "\\|")
two = table(two)
twoCoor = TtoC(names(two),symbol = "\\|")

plot(0,0,type='n',xlim = c(0,50),ylim = c(0,50))
points(a.n(oneCoor[,1]),a.n(oneCoor[,2]),col = rgb(.2,.2,.6,one/max(one)),pch=15,cex=.6)
points(a.n(twoCoor[,1]),a.n(twoCoor[,2]),col = rgb(.2,.6,.2,two/max(two)),pch=15,cex=.6)
points(a.n(oneCoor[,1]),a.n(oneCoor[,2]),col = rgb(.2,.2,.6),pch=15,cex=.6)

Pgrid = expand.grid(seq(.5,49.5,by=1),seq(.5,49.5,by=1))
Pgrid[,3] = rep(0,nrow(Pgrid))
Pgrid[,4] = rep(0,nrow(Pgrid))
Pgrid[match(apply(oneCoor,1,paste,collapse=""),apply(Pgrid[,c(1,2)],1,paste,collapse="")),3] = one/max(one)
Pgrid[match(apply(twoCoor,1,paste,collapse=""),apply(Pgrid[,c(1,2)],1,paste,collapse="")),4] = two/max(two)
points(Pgrid[,1],Pgrid[,2],col = rgb(.2,.2,.6,Pgrid[,3]),pch=15,cex=.6)
points(Pgrid[,1],Pgrid[,2],col = rgb(.2,.6,.2,Pgrid[,4]),pch=15,cex=.6)
resp = rep(0,nrow(Pgrid))

for(i in 1:nrow(Pgrid)){
  bb = Pgrid[i,c(1,2)]
  aa = rbind(bb + c(-1,1),bb + c(0,1),bb + c(1,1),
             bb + c(-1,0),bb + c(0,0),bb + c(1,0), 
             bb + c(-1,-1), bb + c(0,-1),bb + c(1,-1))
  vals = Pgrid[match(apply(aa,1,paste,collapse=""),apply(Pgrid[,c(1,2)],1,paste,collapse="")),-c(1,2)]

  if(sum(apply(vals,2,mean,na.rm=T)) == 0) next
  if(all(apply(vals,2,function(x) sum(x>0,na.rm=T)/sum(!is.na(x))) < .5)) resp[i] = 0
  resp[i] = which.max(apply(vals,2,mean,na.rm=T)) + 2
  progress(i,nrow(Pgrid),100)
  if(i%%100 ==0)points(Pgrid[,1],Pgrid[,2],col = resp,pch=15,cex=.4)
}
points(Pgrid[,1],Pgrid[,2],col = resp,pch=15,cex=.6)

steps = seq(.5,49.5, by = 1)
den1 = den2 = numeric()
for(i in 1:50){
  slice = Pgrid[Pgrid[,1] == steps[i],]
  den1 = rbind(den1,density(slice[,2],weights = slice[,4],n = 50,from = .5, to = 49.5)$y)
  slice = Pgrid[Pgrid[,2] == steps[i],]
  den2 = rbind(den2,density(slice[,1],weights = slice[,4],n = 50,from = .5, to = 49.5)$y)
}
den = (den1+den2)/2
den[den <.1] = 0
image(den)
r <-raster(
  den,
  xmn=.5, xmx=49.5,
  ymn=.5, ymx=49.5, 
  crs=NA
)
plot(r)

spg = data.frame(x = c(Pgrid[,1]),y=c(Pgrid[,2]),z=c(Pgrid[,3]))
coordinates(spg) = ~x+y
gridded(spg) <- TRUE
spg2 = data.frame(x = c(Pgrid[,1]),y=c(Pgrid[,2]),z=c(Pgrid[,4]))
coordinates(spg2) = ~x+y
gridded(spg2) <- TRUE
r = stack(raster(spg),raster(spg2))
plot(r)
nlayers(r)


steps = seq(xmin(r),xmax(r)-1, by = 1) + .5
den1 = den2 = den11 = den22 = numeric()
for(i in 1:50){
  slice = values(r)[coordinates(r)[,1] == steps[i],]
  slCo = coordinates(r)[coordinates(r)[,1] == steps[i],]
  den1  = rbind(den1,density(slCo[,2],weights = slice[,1],n = 50,from = .5, to = 49.5)$y)
  den11 = rbind(den11,density(slCo[,2],weights = slice[,2],n = 50,from = .5, to = 49.5)$y)
  slice = values(r)[coordinates(r)[,2] == steps[i],]
  slCo = coordinates(r)[coordinates(r)[,2] == steps[i],]
  den2  = rbind(den2,density(slCo[,1],weights = slice[,1],n = 50,from = .5, to = 49.5)$y)
  den22 = rbind(den22,density(slCo[,1],weights = slice[,2],n = 50,from = .5, to = 49.5)$y)
}
den = (den1+den2)/2
den12 = (den11+den22)/2
den[den < .1] = 0
den12[den12 < .1] = 0

rD <-raster(
  den,
  xmn=.5, xmx=49.5,
  ymn=.5, ymx=49.5, 
  crs=NA
)
rD2 <-raster(
  den12,
  xmn=.5, xmx=49.5,
  ymn=.5, ymx=49.5, 
  crs=NA
)
r = stack(rD,rD2)

p = raster(spg)
values(p) = apply(getValues(r),1,which.max) - (apply(getValues(r),1,sum) == 0)
plot(p)

library(raster)

path ="/home/mckwit/Dropbox"
allWHR = readRDS(file = paste(path,'/MapShiftingClimates/allWHR_noAm.rds',sep="")) #pick your file
cliCoor = readRDS(file = paste(path,'/MapShiftingClimates/cliCoor_noAm.rds',sep=""))



##This normalizes state counts to (0,1)
sts = sort(unique(allWHR[,'state']))
typ = sort(unique(allWHR[,'type']))
for(i in 1:len(typ)){
  for(j in 1:len(sts)){
    whrST = which(allWHR[,'type'] == typ[i] & allWHR[,'state'] == sts[j] )
    if(len(whrST)==0) next
    allWHR[whrST,'abund'] = allWHR[whrST,'abund']/max(allWHR[whrST,'abund'],na.rm=T)
  }
  print(i)
}
saveRDS(allWHR,file = paste(path,'/MapShiftingClimates/allWHR2_noAM.rds',sep=""))
allWHR = readRDS(file = paste(path,'/MapShiftingClimates/allWHR2_noAM.rds',sep=""))
cliCoor = readRDS(file = paste(path,'/MapShiftingClimates/cliCoor_noAm.rds',sep=""))
basRast = raster(paste(path,"/bioclime/bio_10m_bil/bio1.bil",sep=""))
cur = allWHR[allWHR[,'type'] == 'whr70I',]
sts = sort(unique(cur[,'state']))
whr = a.n(a.c(cur[,'whr']))
wts = c(4,1,1,1,1,1.414,1.414,1.414,1.414,2,2,2,2)

ext = extent(basRast)
ext@xmax = -40
#ext@xmin = -135 
#ext@ymax = 52
ext@ymin = 10  #25
basRast = crop(basRast, ext)

crsRast = crs(basRast)
coorRast = coordinates(basRast)

tmp = reshape(cur[,-4],direction = 'wide',timevar= 'state',idvar='whr')
colnames(tmp) = gsub("abund.","",colnames(tmp))

subCoor = cliCoor[a.n(a.c(tmp[,'whr'])),]

inRast = match(apply(subCoor,1,function(x) paste(x,collapse="|")),apply(coorRast,1,function(x) paste(x,collapse="|")))

whrMat = cbind(match(coorRast[,1],sort(unique(coorRast[,1]))),match(coorRast[,2],sort(unique(coorRast[,2]))))
whrMatC = apply(whrMat,1,function(x) paste(x,collapse="|")) 

whrMat11 = match(paste(whrMat[,1]+1,whrMat[,2],sep="|"),whrMatC)
whrMat12 = match(paste(whrMat[,1]-1,whrMat[,2],sep="|"),whrMatC)
whrMat13 = match(paste(whrMat[,1]  ,whrMat[,2]+1,sep="|"),whrMatC)
whrMat14 = match(paste(whrMat[,1]  ,whrMat[,2]-1,sep="|"),whrMatC)

whrMat21 = match(paste(whrMat[,1]+1,whrMat[,2]+1,sep="|"),whrMatC)
whrMat22 = match(paste(whrMat[,1]+1,whrMat[,2]-1,sep="|"),whrMatC)
whrMat23 = match(paste(whrMat[,1]-1,whrMat[,2]+1,sep="|"),whrMatC)
whrMat24 = match(paste(whrMat[,1]-1,whrMat[,2]-1,sep="|"),whrMatC)

whrMat31 = match(paste(whrMat[,1]+2,whrMat[,2],sep="|"),whrMatC)
whrMat32 = match(paste(whrMat[,1]-2,whrMat[,2],sep="|"),whrMatC)
whrMat33 = match(paste(whrMat[,1]  ,whrMat[,2]+2,sep="|"),whrMatC)
whrMat34 = match(paste(whrMat[,1]  ,whrMat[,2]-2,sep="|"),whrMatC)

#whrVal = which(!is.na(inRast))
whrVal = inRast

tmpVals = matrix(NA,nrow(coorRast),len(sts))
for(i in 1:len(whrVal)){
  tmpVals[whrVal[i],] = a.n(tmp[i,-1])
  if(i%%1000 == 0)print(i)
}

tmpVals[is.na(tmpVals)] = 0
for(i in 1:3){
t11 = tmpVals[whrMat11,] ; t12 = tmpVals[whrMat12,] ; t13 = tmpVals[whrMat13,] ; t14 = tmpVals[whrMat14,]
t21 = tmpVals[whrMat21,] ; t22 = tmpVals[whrMat22,] ; t23 = tmpVals[whrMat23,] ; t24 = tmpVals[whrMat24,]
t31 = tmpVals[whrMat31,] ; t32 = tmpVals[whrMat32,] ; t33 = tmpVals[whrMat33,] ; t34 = tmpVals[whrMat34,]

t11[is.na(t11)] = 0 ; t12[is.na(t12)] = 0 ; t13[is.na(t13)] = 0 ; t14[is.na(t14)] = 0
t21[is.na(t21)] = 0 ; t22[is.na(t22)] = 0 ; t23[is.na(t23)] = 0 ; t24[is.na(t24)] = 0
t31[is.na(t31)] = 0 ; t32[is.na(t32)] = 0 ; t33[is.na(t33)] = 0 ; t34[is.na(t34)] = 0

tmpVals = tmpVals * 4 + t11 + t12 + t13 + t14 + 
                     (t21/1.414) + (t22/1.414) + (t23/1.414) + (t24/1.414) + 
                     (t31/2) + (t32/2) + (t33/2) + (t34/2) 
}
tmpVals[tmpVals == 0 ] = NA
tmpVals =  apply(tmpVals,1,which.max)
tmpVals = a.n(tmpVals)

#world.shape<-readOGR(dsn=paste(path,"/Shapes/TM_WORLD_BORDERS_SIMPL-0.3",sep=""),layer="TM_WORLD_BORDERS_SIMPL-0.3")
#worldCRS = crs(world.shape)

zerRast = basRast
whr = which(is.na(values(basRast)))
values(zerRast) = tmpVals

values(zerRast)[whr] = rep(NA, length(whr))
#values(zerRast)[whr] = rep(0, length(whr))
plot(zerRast)
tmpShp = rasterToPolygons(zerRast, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=TRUE)
tmpShp@data[[1]] <-  colnames(tmp)[-1][tmpShp@data[[1]]]



factpal <- colorFactor(topo.colors(52), as.factor(tmpShp@data[[1]]))
library(leaflet)
  leaflet(tmpShp) %>%
  addTiles() %>%
  addPolygons(
    stroke =  ~factpal(as.factor(tmpShp@data[[1]])), fillOpacity = 0.2, smoothFactor = 0.5,
    color =  ~factpal(as.factor(tmpShp@data[[1]])),
    popup = ~tmpShp@data[[1]]
  )
  
palette(c("#a6cee3", "#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6"))
plot(tmpShp,col= tmpShp@data[[1]],border = rev(tmpShp@data[[1]]))









allWHR = readRDS(file = paste(path,'/MapShiftingClimates/allWHR2_noAM.rds',sep=""))
whrCliCoor = match(apply(cliCoor,1,function(x) paste(x,collapse = "|")),apply(coorRast,1,function(x) paste(x,collapse = "|")))


#non Vectorized cleanup

cur = allWHR[allWHR[,'type'] == 'whrCU',]
sts = sort(unique(cur[,'state']))
whr = a.n(a.c(cur[,'whr']))
wts = c(4,1,1,1,1,1.414,1.414,1.414,1.414,2,2,2,2)/4
#need to scale abund to state [0,1]

for(k in 1:length(whr)){
  whr2 = which(coorRast[,1] == cliCoor[whr[k],1] & coorRast[,2] == cliCoor[whr[k],2] )
  basis = whrMat[whr2,]
  adj = rbind(basis,c(basis[1]+1,basis[2]),c(basis[1]-1,basis[2]),c(basis[1],basis[2]+1),c(basis[1],basis[2]-1),
              basis - 1,basis + 1,c(basis[1]+1,basis[2]-1),c(basis[1]-1,basis[2]+1),
              c(basis[1]+2,basis[2]),c(basis[1]-2,basis[2]),c(basis[1],basis[2]+2),c(basis[1],basis[2]-2))
  
  whrCC = match(match(apply(adj,1,function(x)  paste(x,collapse = "|")),whrMatC),whrCliCoor)
  whrSt = match(whrCC, a.n(a.c(cur[,'whr'])) )
  #counts = a.n(table(cur[whrSt,'state']))
  weights = tapply(cur[whrSt,'abund']* wts,cur[whrSt,'state'],mean,na.rm=T)
  #weights[is.na(weights)] = 0
  if(FALSE){
    counts = weights = rep(0,len(sts))
  for(i in 1:nrow(adj)){
    #adj[i,] = coorRast[which(whrMat[,1] == adj[i,1] & whrMat[,2] == adj[i,2]),] #from xy to lat lon
    whrCC = which(cliCoor[,1] == adj[i,1] & cliCoor[,2] == adj[i,2] ) #from lat lon raster to lat lon data
    whrSt = which(a.n(a.c(cur[,'whr'])) == whrCC) #where in allWHR
    whrCO = match(cur[whrSt,'state'],sts) #which states
    counts[a.n(names(table(whrCO)))]  = counts[a.n(names(table(whrCO)))] + a.n(table(whrCO))
    if(len(whrSt) == 0)next
    weights[a.n(names(table(whrCO)))] = weights[a.n(names(table(whrCO)))] + cur[whrSt,'abund'] * wts[i]
  }
  }
  values(zerRast)[whr2] = which.max(weights)
  progress(k,length(whr),1000)
}



cur = allWHR[allWHR[,'type'] == 'whrCU',]

r = raster()
for(i in 1: length(sts)){
  den = cur[cur[,'state'] == sts[i],]
  whr = a.n(a.c(den[,'whr']))
  tmp =data.frame(x = cliCoor[whr,1],y = cliCoor[whr,2], z = den[,'abund']/max(den[,'abund']))
  values(zerRast)[match(paste(tmp$x,tmp$y),paste(coorRast[,1],coorRast[,2]))] = tmp$z  
  coordinates(tmp) = ~x+y
  gridded(tmp) <- TRUE
  r2 <-raster(tmp)
  r = stack(r,r2)
}

library(maps)
par(mfrow= c(1,1),mar=c(2,2,2,1))
plot(extent(r2),xlab='',ylab='',bty='n',xaxt='n',yaxt='n')
map(,add=T,col=1)
map('state',add=T,col=1)
plot(r2, add =T)

spatDen = function(r,maskRast){
rY = rX = r
stepy = sort(unique(coordinates(r)[,1]))
stepx = sort(unique(coordinates(r)[,2]))
#for(i in 1:length(stepy)){
#  print(i)
#  whr = coordinates(r)[,1] == stepy[i]
#  slice = values(r)[whr]
#  if(sum(slice) == 0) next
#  coord = coordinates(r)[whr,]
#  values(rY)[whr]  = density(coord[,2],weights = slice,n = len(whr),from = min(coord), to = max(coord))$y
#}
for(i in 1:length(stepx)){
  print(i)
  whr = coordinates(r)[,2] == stepx[i]
  slice = values(r)[whr]
  if(sum(slice) == 0) next
  coord = coordinates(r)[whr,]
  values(rX)[whr]  = density(coord[,1],weights = slice,n = len(whr),from = min(coord), to = max(coord))$y
}

r = (rX+rY)/2
values(r)[values(r) < .1] = 0

}


out = kde2d.weighted(coordinates(r)[1:100000,1],coordinates(r)[1:100000,2],w = values(r)[1:100000],n = len(values(r)[1:100000]))

library(MASS)
kde2d.weighted <- function (x, y, w, h, n = 25, lims = c(range(x), range(y))) {
  nx <- length(x)
  if (length(y) != nx) 
    stop("data vectors must be the same length")
  if (length(w) != nx & length(w) != 1)
    stop("weight vectors must be 1 or length of data")
  gx <- seq(lims[1], lims[2], length = n) # gridpoints x
  gy <- seq(lims[3], lims[4], length = n) # gridpoints y
  if (missing(h)) 
    h <- c(bandwidth.nrd(x), bandwidth.nrd(y));
  if (missing(w)) 
    w <- numeric(nx)+1;
  h <- h/4
  ax <- outer(gx, x, "-")/h[1] # distance of each point to each grid point in x-direction
  ay <- outer(gy, y, "-")/h[2] # distance of each point to each grid point in y-direction
  z <- (matrix(rep(w,n), nrow=n, ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n, nx))/(sum(w) * h[1] * h[2]) # z is the density
  return(list(x = gx, y = gy, z = z))
}