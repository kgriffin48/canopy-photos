# coded by Francesco CHIANUCCI, 16th Feb 2018 (fchianucci@gmail.com)

#the code take a fullframe fisheye image, threshold it, and calculate  
#effective leaf area (Le), leaf area index (L), foliage clumping (LX)
#for image processed using 9 zenith rings and 8 azimuth segments, the
#Miller (1967) theorem, and the gap averaging method of Lang and Xiang (1986) 

rm(list=ls())

library(raster)
library(rtiff)
library(jpeg)

deg2rad <- function(x) x*pi/180
dist2Dorig <- function(x0, x, y0, y) sqrt((y-y0)^2 + (x-x0)^2)

#set the directory in which you have fullframe images:
setwd('/Users/francescochianucci/Dropbox')
dir()

files<-list.files()
#set the common extension (defaul .pngs)
dbf.files <- files[grep(".JPG", files, fixed=T)]
# k=1
bo=NULL
tdf=NULL

#thresholding image
for (k in 1:length(dbf.files)){
  img <- raster(dbf.files[k],band=3)
  threshold_all <- autoThreshold(values(img))
  th <- threshold_all[3]
  mask<-which(values(img)>=th)
  values(img)<-0 
  values(img)[mask]<-1# sky pixels
  writeJPEG(as.matrix((img)),paste0('thd_',dbf.files[k]))  
  
  
  dim(img)
  pxdf <- as.data.frame(coordinates(img))  
  pxdf$DN <- 0 # trees
  pxdf$DN[mask] <- 1 # sky
  
  Xmax <- dim(img)[2]
  Ymax <- dim(img)[1]
  
  Xc<-Xmax/2
  Yc<-Ymax/2
  Rad<-trunc(sqrt(Xc^2+Yc^2))
  
  pxdf$dist <- dist2Dorig(Xc,pxdf$x, Yc,pxdf$y)
  
  #coordinates traslated to origin(Xc,Yc):
  pxdf$dx <- pxdf$x-Xc
  pxdf$dy <- pxdf$y-Yc
  
  #8 azimuth segments
  pxdf$theta <- NA
  pxdf$theta[pxdf$dx>0 & pxdf$dy>=0] <- atan(pxdf$dy[pxdf$dx>0 & pxdf$dy>=0]/pxdf$dx[pxdf$dx>0 & pxdf$dy>=0])
  pxdf$theta[pxdf$dx>0 & pxdf$dy<0] <- atan(pxdf$dy[pxdf$dx>0 & pxdf$dy<0]/pxdf$dx[pxdf$dx>0 & pxdf$dy<0])+2*pi
  pxdf$theta[pxdf$dx<0] <- atan(pxdf$dy[pxdf$dx<0]/pxdf$dx[pxdf$dx<0])+pi
  pxdf$theta[pxdf$dx==0 & pxdf$dy>0] <- pi/2
  pxdf$theta[pxdf$dx==0 & pxdf$dy<0] <- pi*3/2
  
  #9 zenith rings
  r_min=seq(0,trunc(Rad-(Rad/9)),trunc(Rad/9))
  r_max=seq(trunc(Rad/9),Rad,trunc(Rad/9))
  
  rdf <- data.frame(rmin=rep(r_min,each=8), rmax=rep(r_max,each=8),
                    alpha.from=rep(seq(0,7/4*pi,pi/4),length(r_min)), alpha.to=rep(seq(pi/4,2*pi,pi/4),length(r_min)))
  rdf$idsect <- paste('sect',10000+(1:nrow(rdf)),sep='')
  rdf$ring<-rep(seq(10,90,10),each=8)
  rdf$IM<-rep(dbf.files[k],nrow(rdf))
  
  pxdf$settore <- 'tbd'
  head(pxdf)
  j <- 1
  for(j in 1:nrow(rdf)){
    pxdf$settore[pxdf$dist>=rdf$rmin[j] & pxdf$dist<rdf$rmax[j] & 
                   pxdf$theta>=rdf$alpha.from[j] & pxdf$theta<rdf$alpha.to[j]] <- rdf$idsect[j]
    if ((j-1) == (trunc((j-1)/3)*3)){
      print(paste("seg ",j," ","image ", k ," ",date()))
    }
  }
  
  rdf$GF<-round(tapply(pxdf$DN, pxdf$settore, mean),2)[1:nrow(rdf)]
  bo<-data.frame(rbind(bo,rdf))
}



bo$conc <- paste(bo$IM, bo$ring,sep='##')
bo$alpha.from<-bo$alpha.from*180/pi
bo$alpha.to<-bo$alpha.to*180/pi

bo2<-reshape(bo,  v.names="GF",idvar="conc", timevar="alpha.to",direction="wide")
bo2<-bo2[,-c(3,4)]
for(z in 1:nrow(bo2)){bo2$NULL_SEG[z]<-length(bo2[z,6:13][bo2[z,6:13]==0])}
bo2$GF<-apply(bo2[,6:13],1,mean,na.rm=TRUE)

#Plot GF distribution
plot(bty="n",bo2$ring,bo2$GF,main="Gap fraction distribution",xlab="Zenith angle",ylab=expression(P(theta)),xlim=c(0,95),ylim=c(0,1.2),xaxs="i",yaxs="i",xaxt="n",yaxt="n",col=as.factor(bo2$IM))
axis(2, las=2)
axis(1,at=seq(0,90,15))

#replace empty gap segments
for (x in 1:length(bo2[,1]))
bo2[x,6:13][bo2[x,6:13]==0|is.na(bo2[x,6:13])]<-apply(bo2[x,6:13],1,mean,na.rm=TRUE)

#calculate clumping index
LX=NULL
for (j in 1:length(bo2[,1]))
  LX<-c(LX,log(mean(as.numeric(bo2[j,6:13])))/mean(log(as.numeric(bo2[j,6:13]))))
LX[which(is.na(LX))]<-1

#plot Clumping distribution (LX)
plot(bty="n",bo2$ring,LX,main="Foliage Clumping distribution",xlab="Zenith angle",ylab=expression(Omega(theta)),xlim=c(0,95),ylim=c(0.4,1.2),xaxs="i",yaxs="i",xaxt="n",yaxt="n",col=as.factor(bo2$IM))
axis(2, las=2)
axis(1,at=seq(0,90,15))

#calculate effective leaf area Le
Le=NULL
for (w in 1:length(bo2[,1]))
  Le<-c(Le,-log(mean(as.numeric(bo2[w,6:13])))*cos(bo2[w,3]/2*pi/180)*2)

#calculate leaf area L sensu Lang and Xiang (1986)
L=NULL
for (y in 1:length(bo2[,1]))
L<-c(L,mean(-log(as.numeric(bo2[y,6:13])))*cos(bo2[y,3]/2*pi/180)*2)

# normalization weighting for view zenith angles (VZA)
bo2$Le<-Le
bo2$L<-L
bo2$LX<-Le/L
bo2$sinVZA<-sin(bo2$ring/2*pi/180)
bo2$w<-(bo2$sinVZA)/(sum(bo2$sinVZA[1:9]))


#Calculate weighted Le, L, LX
bo2$LXw<-round(LX,2)*bo2$w
bo2$Lew<-round(Le,2)*bo2$w
bo2$Lw<-round(L,2)*bo2$w
bo2$GFw<-bo2$GF*bo2$w

#Calculate Image-level Le, L, LX
(tdf<-data.frame(rbind(tdf,cbind(LXw=tapply(bo2$LXw,bo2$IM,sum),Le=tapply(bo2$Lew,bo2$IM,sum),L=tapply(bo2$Lw,bo2$IM,sum),DIFN=tapply(bo2$GFw,bo2$IM,sum)))))




# write.csv(bo2,"bo2.csv")
  