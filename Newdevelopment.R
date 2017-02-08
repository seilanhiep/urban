library(raster)
#*******************************************************#
###############################################################
###########################************Focal statistic*************####################################

Urbanmanifest= function(currenturbanfile,previousurbanfile,outdrive,threshold1,threshold2){
if(missing(urbandrive)){cat("Please Enter Image Resolution:", "\n")}
urbextt2		=	raster(currenturbanfile)
urbextt1		=	raster(previousurbanfile)

#urbext		=	raster(urbext,layer=numclassurban)
#filter=matrix(c(0,1,0,1,-4,1,0,1,0), nrow=3)		#Laplacian filter
#filter=matrix(c(1,2,1,0,0,0,-1,-2,-1) / 4, nrow=3)	#Sobel filter

####### Manifestation Metrics for measuring urban sprawl #############
neighbor	=	33
neighbor2	=	3
urbnesst1	=	urbextt1
nfocal	=	matrix(rep(c(1:1),times=(neighbor^2)),neighbor,neighbor)
ft1		=	focal(urbnesst1, w=nfocal, pad=TRUE, padValue=0)
urb1kmt1	=	ft1*urbnesst1

### create directory ###
ddri		=	unlist(strsplit(currenturbanfile, split='', fixed=TRUE))[1]
dirlsat	=	substr(currenturbanfile, grep(ddri,currenturbanfile), (nchar(currenturbanfile)- nchar(gsub('^.*/','',currenturbanfile))) -1)
dir.create(paste0(dirlsat,"/",outdrive))
dir = paste0(dirlsat,"/",outdrive)

## statistic cell computation ###
urbmaxt1		=	cellStats(urb1kmt1, stat='max', na.rm=TRUE, asSample=TRUE)
urbcomnesst1	=	(urb1kmt1/urbmaxt1)*100
murbcoret1		=	(urbcomnesst1 >= threshold2)
fmaint1		=	focal(murbcoret1, w=nfocal, fun=max, pad=TRUE, padValue=0) #mean, modal, min or max
surbcoret1		=	(fmaint1*(urbcomnesst1 != 0 & urbcomnesst1 < threshold2))- murbcoret1
surbcoret1		=	(surbcoret1 >0)
urbfringet1		=	((urbcomnesst1 >= threshold1 & urbcomnesst1 < threshold2)-surbcoret1)
urbfringet1		=	(urbfringet1 >0)
nfocal2t1		=	matrix(rep(c(1:1),times=(neighbor2^2)),neighbor2,neighbor2)
f2t1			=	focal((urbcomnesst1 != 0 & urbcomnesst1 < 30), w=nfocal2t1, pad=TRUE, padValue=0)
ribbondevt1		=	(((f2t1>=7)*(urbcomnesst1 != 0 & urbcomnesst1 < threshold1))-urbfringet1-surbcoret1)
ribbondevt1		=	(ribbondevt1 >0)
scatterdevt1	=	(((urbcomnesst1 != 0 & urbcomnesst1 < threshold1))-ribbondevt1-surbcoret1)
scatterdevt1	=	(scatterdevt1 >0)
### compute urban manifestration ###
urbmanifestt1	=	murbcoret1+(surbcoret1*2)+(urbfringet1*3)+(ribbondevt1*4)+(scatterdevt1*5)
	
####### Metrics for measuring urban extent #############
### impervious surface ###
builtupt1		=	urbextt1
builtup1kmt1	=	ft1/cellStats(ft1, stat='max')*100
urbost1		=	(((urbnesst1 ==0)*(builtup1kmt1>=50))- builtupt1)
urbost1		=	(urbost1>0)
urbareat1		=	builtupt1 + urbost1*2
builtup100mt1	=	focal(builtupt1, w=nfocal2t1, pad=TRUE, padValue=0)
perost1		=	(((urbnesst1 ==0)*(builtup100mt1>0))- builtupt1 - urbost1)
perost1		=	(perost1>0)	
urbfootprintt1	=	builtupt1 + urbost1*2 + perost1*3
ost1			=	urbost1 + (perost1*2)

###### New development metrics ###
newdev		=	(urbextt2==1 & urbextt1==0)
newdev		=	(newdev>0)
infill		=	(newdev)*(urbost1)
leapfrog		=	(newdev)*(urbfootprintt1==0)
extdev		=	(newdev)-(infill)-(leapfrog)
newdevelopment	=	infill+extdev*2+leapfrog*3

ndf			=	data.frame(freq(newdevelopment))
ndf$count		=	(round(newdevfreq[,2]*900/1000000,digit=2))
ndf$percet		=	round(ndf[1:nrow(ndf),2]/sum(ndf[2:nrow(ndf),2])*100,digit=2)
ndf

### plot urbanform ###
plot(urbextt1,col=colorRampPalette(1:3)(255),main="Previous Urban Pixel")
plot(urbcomnesst1,col=colorRampPalette(1:8)(255), main="urban compactness")
plot(urbmanifestt1,col=colorRampPalette(1:5)(255), main="urban manifestration")
plot(urbfootprintt1,col=colorRampPalette(1:3)(255), main="urban urbfootprint")
plot(ost1,col=colorRampPalette(1:3)(255), main="urban os")
plot(newdevelopment,col=colorRampPalette(1:3)(255), main="urban new development type")
### write image to file ###
writeRaster(urb1kmt1,filename=paste0(dir,"/","urbancompactness.img"),format="HFA",overwrite=TRUE)
writeRaster(urbmanifestt1,filename=paste0(dir,"/","urbanmanifestration.img"),format="HFA",overwrite=TRUE)
writeRaster(urbfootprintt1,filename=paste0(dir,"/","urbanfootprint.img"),format="HFA",overwrite=TRUE)
writeRaster(ost1,filename=paste0(dir,"/","Openspace.img"),format="HFA",overwrite=TRUE)
writeRaster(newdevelopment,filename=paste0(dir,"/","newdevelopmenttype.img"),format="HFA",overwrite=TRUE)
}

system.time(Urbanmanifest( "D:/tcR/urban2016.img",
				"D:/tcR/urban2000.img",
				"newdevelopment",
				30,
				50))



