library(raster)
#*******************************************************#
#############################################################################################
### in case you want to write urban extent from land use ###
#urbext		=	raster("D:/tcR/LC81260522016050LGN00/class/superclass2016.img")
#urbext		=	raster("D:/tcR/LE71260522002003BKT00/class/superclass2016.img")
#plot(urbext)
#urbext		=	(urbext==2 | urbext==4)
#writeRaster(urbext,filename=paste0("D:/tcR/urbanfile","/","urban2002.img"),format="HFA",overwrite=TRUE)
#filter=matrix(c(0,1,0,1,-4,1,0,1,0), nrow=3)		#Laplacian filter
#filter=matrix(c(1,2,1,0,0,0,-1,-2,-1) / 4, nrow=3)	#Sobel filter
###############################################################
###########################************Urban manifestration current year*************####################################

Urbanmanifest= function(classedimage,outdrive,neighbor,neighbor2,threshold1,threshold2,outputurbcompact,outputurbform){
if(missing(classedimage)){cat("Please Enter Image Resolution:", "\n")}
urbext		=	raster(classedimage)
#filter=matrix(c(0,1,0,1,-4,1,0,1,0), nrow=3)		#Laplacian filter
#filter=matrix(c(1,2,1,0,0,0,-1,-2,-1) / 4, nrow=3)	#Sobel filter

####### Manifestation Metrics for measuring urban sprawl #############
urbness	=	urbext
nfocal	=	matrix(rep(c(1:1),times=(neighbor^2)),neighbor,neighbor)
f		=	focal(urbness, w=nfocal, pad=TRUE, padValue=0)
urb1km	=	f*urbness

### create directory ###
ddri		=	unlist(strsplit(classedimage, split='', fixed=TRUE))[1]
dirlsat	=	substr(classedimage, grep(ddri, classedimage), (nchar(classedimage)- nchar(gsub('^.*/','',classedimage))) -1)
dir.create(paste0(dirlsat,"/",outdrive))
dir = paste0(dirlsat,"/",outdrive)
## statistic cell computation ###
urbmax	=	cellStats(urb1km, stat='max', na.rm=TRUE, asSample=TRUE)
urbcomness	=	(urb1km/urbmax)*100
murbcore	=	(urbcomness >= threshold2)
fmain		=	focal(murbcore, w=nfocal, fun=max, pad=TRUE, padValue=0) #mean, modal, min or max
surbcore	=	(fmain*(urbcomness != 0 & urbcomness < threshold2))- murbcore
surbcore	=	(surbcore >0)
urbfringe	=	((urbcomness >= threshold1 & urbcomness < threshold2)-surbcore)
urbfringe	=	(urbfringe >0)
nfocal2	=	matrix(rep(c(1:1),times=(neighbor2^2)),neighbor2,neighbor2)
f2		=	focal((urbcomness != 0 & urbcomness < 30), w=nfocal2, pad=TRUE, padValue=0)
ribbondev	=	(((f2>=7)*(urbcomness != 0 & urbcomness < threshold1))-urbfringe-surbcore)
ribbondev	=	(ribbondev >0)
scatterdev	=	(((urbcomness != 0 & urbcomness < threshold1))-ribbondev-surbcore)
scatterdev	=	(scatterdev >0)
### compute urban manifestration ###
urbmanifest	=	murbcore+(surbcore*2)+(urbfringe*3)+(ribbondev*4)+(scatterdev*5)
	
####### Metrics for measuring urban extent #############
### impervious surface ###
builtup			=	urbext
builtup1km			=	f/cellStats(f, stat='max')*100
urbos				=	(((urbness ==0)*(builtup1km>=50))- builtup)
urbos				=	(urbos>0)
urbarea			=	builtup + urbos*2
builtup100m			=	focal(builtup, w=nfocal2, pad=TRUE, padValue=0)
peros				=	(((urbness ==0)*(builtup100m>0))- builtup - urbos)
peros				=	(peros>0)	
urbfootprint		=	builtup + urbos*2 + peros*3
os				=	urbos + (peros*2)

### indices ###
urbmat			=	data.frame(freq(urbmanifest))
names(urbmat)		=	c("urbantype","pixel")
urbmat			=	rbind(urbmat,sum(urbmat[2:nrow(urbmat),2]))
urbmat$sizekm2		=	(round(urbmat[,2]*900/1000000,digit=2))
urbmat$percent		=	round(urbmat[1:nrow(urbmat),2]/sum(urbmat[2:(nrow(urbmat)-1),2])*100,digit=2)
urbmat[1:nrow(urbmat),1]=	c("non urban","main urban core","secondary urban core","urban fringe","ribbon development","scatter development","builtup area")
print(urbmat)
#write(t(as.matrix(urbmat)),file=paste0(dir,"/","urban manifest.csv"),append=TRUE,sep="\"")
### plot urbanform ###
plot(urbext,col=colorRampPalette(1:3)(255),main="Urban Pixel")
plot(urbcomness,col=colorRampPalette(1:8)(255), main="urban compactness")
plot(urbmanifest,col=colorRampPalette(1:5)(255), main="urban manifestration")
plot(urbfootprint)
plot(os)
### write image to file ###
writeRaster(urb1km,filename=paste0(dir,"/","urbancompactness.img"),format="HFA",overwrite=TRUE)
writeRaster(urbmanifest,filename=paste0(dir,"/","urbanmanifestration.img"),format="HFA",overwrite=TRUE)
writeRaster(urbfootprint,filename=paste0(dir,"/","urbanfootprint.img"),format="HFA",overwrite=TRUE)
writeRaster(os,filename=paste0(dir,"/","Openspace.img"),format="HFA",overwrite=TRUE)
}

system.time(Urbanmanifest("D:/tcR/urbanfile/urban2016.img",
				"urbanform",
				33,
				3,
				30,
				50))

###########################************Urban manifestration previous year*************####################################
###########################************New development between two years*************####################################

Urbanmanifest= function(currenturbanfile,previousurbanfile,outdrive,threshold1,threshold2){
if(missing(urbandrive)){cat("Please Enter Image Resolution:", "\n")}
library(raster)
urbextt2		=	raster(currenturbanfile)
urbextt1		=	raster(previousurbanfile)
comp			=	compareRaster(urbextt2, urbextt1, stopiffalse=FALSE)
if( comp == TRUE){
urbextt2		=	urbextt2
urbextt1		=	urbextt1
}else{
ext			=	c(xmin(urbextt1),xmin(urbextt2),xmax(urbextt1),xmax(urbextt2),ymin(urbextt1),ymin(urbextt2),ymax(urbextt1),ymax(urbextt2))
ext			=	sort(ext, decreasing = FALSE)
cropext		=	c(ext[2],ext[3],ext[6],ext[7])
urbextt2		=	crop(urbextt2,cropext)
urbextt1		=	crop(urbextt1,cropext)
}
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

### indices ###
urbmat			=	rbind(data.frame(freq(urbmanifestt1)),data.frame(freq(urbfootprintt1)))
names(urbmat)		=	c("urbantype","pixel")
urbmat			=	rbind(urbmat,sum(urbmat[8:nrow(urbmat),2]),sum(urbmat[9:nrow(urbmat),2]))
urbmat$sizekm2		=	(round(urbmat[,2]*900/1000000,digit=2))
urbmat$percent		=	0
urbmat[1:6,4]		=	round(urbmat[1:6,2]/sum(urbmat[2:6,2])*100,digit=2)
urbmat[7:nrow(urbmat),4]=	round(urbmat[7:nrow(urbmat),2]/urbmat[8,2]*100,digit=2)
urbmat[1:nrow(urbmat),1]=	c("non urban","main urban core","secondary urban core","urban fringe","ribbon development","scatter development",
					"non urban excluded OS","builtup area","urban OS","periperial OS","urban footprint","open space")
print(urbmat)
#write(t(as.matrix(urbmat)),file=paste0(dir,"/","urban matrices.csv"),append= FALSE,sep="\"")

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

system.time(Urbanmanifest( "D:/tcR/urbanfile/urban2016.img",
				"D:/tcR/urbanfile/urban2002.img",
				"newdevelopment",
				30,
				50))

