library(maps)
library(mapdata)
library(mapplots)
par(mar = c(0,0,0,0)) ## set plot window margins
par(pin = c(6,6)) ## set plot window size
mycol = rainbow(12, alpha = 0.6) # set colours

setwd("C:/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/")


cords <- read.csv("C:/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/R stats/Adegenet IBD/Mantelcoordinates.csv", header = T)  ## load my microsat coordinates file. 

#geocent <- read.csv("C:/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/R stats/mapplots/geo cluster centre.csv")  ##  Loads my geographic centre data

points(geocent$lat, geocent$lon, pch=19, col=mycol, cex=3)  ## Plots the geocent data on the map.

mtcords <- read.csv("C:/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/R stats/Adegenet IBD/MtDNAcoordinates.csv", header = T)  ## Loads my mtDNA sample coordinates

points(cords$lon, cords$lat, pch=19, col="black", cex=0.5)  ## Plots my microsat sample locations on the map 

points(mtcords$lat, mtcords$lon, pch=19, col="blue", cex=0.5)  ## Plots my mtDNA Genbank sample locations on the map 



pies <- read.csv("/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Complete dataset outputs/Whole EU 22 clusters/WHOLE_EU_mean_clustermemberships.csv", header=T) ## have read this in again and specified that there are headers as I was having difficulty assigning headers to the object. This allows me to call populattions using the $ operator as below.

pies[13,2:50] <- 1

## Plot map ##

map("worldHires", xlim=c(-10, 45), ylim=c(43,70), col="gray90", fill=TRUE)##plots the area of the map that I am interested in

## Add sample locations ## 

points(cords$lon, cords$lat, pch=19, col="black", cex=0.5)  ## Plots my sample locations on the map (just for fun). But dont know how to lable, yet!
points(mtcords$lat, mtcords$lon, pch=19, col="blue", cex=0.5)  ## Plots my mtDNA Genbank sample locations on the map (just for fun). But dont know how to lable, yet!


## Now plot the pies, note that these pies have been offset from their actual population locations so as to allow for minimum overlap ##




mycol1 = adjustcolor("red", alpha.f = 0.6)
mycol2 = adjustcolor("blue", alpha.f = 0.6)
mycol3 = adjustcolor("green", alpha.f = 0.6)
mycol4 = adjustcolor("purple", alpha.f = 0.6)
mycol5 = adjustcolor("orange", alpha.f = 0.6)
mycol6 = adjustcolor("black", alpha.f = 0.6)
mycol7 = adjustcolor("yellow", alpha.f = 0.6)
mycol8 = adjustcolor("pink", alpha.f = 0.6)
mycol9 = adjustcolor("brown", alpha.f = 0.6)
mycol10 = adjustcolor("white", alpha.f = 0.3)



## Finish making these colours and assign the right pops to the right pools!

## UK Pies ##

add.pie(pies$HOLT[13],x=0.4,y=54.5,labels="",radius=0.847,edges=200,clockwise=T, col = mycol1)

add.pie(pies$CCS[13],x=0.1279768,y=50.5 ,labels="",radius=0.664541667,edges=200,clockwise=T, col = mycol1)

add.pie(pies$CAKE[13],x= -1.9,y=53.8,labels="",radius=0.628833334,edges=200,clockwise=T, col = mycol1)

add.pie(pies$BF[13],x=-2,y=52.3,labels="",radius=0.715833334,edges=200,clockwise=T, col = mycol1)

add.pie(pies$RM[13],x=2.3,y=51.,labels="",radius=0.776166667,edges=200,clockwise=T, col = mycol1)

add.pie(pies$OTOM[13],x=3.8,y=52.1 ,labels="",radius=0.635333334,edges=200,clockwise=T, col = mycol1)

add.pie(pies$RAIL[13],x=3.1,y=54.5,labels="",radius=0.77325,edges=200,clockwise=T, col = mycol1)

add.pie(pies$MOAT[13],x=4.9,y=53.3 ,labels="",radius=0.719875,edges=200,clockwise=T, col = mycol1)

add.pie(pies$FFF[13],x=-2.5,y=50.5,labels="",radius=0.711208334,edges=200,clockwise=T, col = mycol1)

## Belgian Pies ##

add.pie(pies$BOK[13],x=7.5,y=50,labels="",radius=0.711208334,edges=200,clockwise=T, col = mycol1)

add.pie(pies$MVW[13],x=2.8,y=49.2,labels="",radius=0.739416667,edges=200,clockwise=T, col = mycol1)
add.pie(pies$MVWZ[13],x=5.2,y=49.6,labels="",radius=0.737166667,edges=200,clockwise=T, col = mycol1)

## N. German pie ##

add.pie(pies	$	FFG[13]	,	x	=	7.558594	,	y	=	51.890053	,	labels	=	""	,	radius	=	1.186083334	,	edges	=	200,	clockwise	=	T, col = mycol1)																																			


####### Pool 2 #######

### Danish Pies ###

add.pie(pies  $  GAM[13]  ,	x	=	9.5	,	y	=	56.5	,	labels	=	""	,	radius	=	1.05	,	edges	=	200,	clockwise	=	T, col = mycol2)																																			
add.pie(pies  $  PED[13]	,	x	=	12.34	,	y	=	55.73	,	labels	=	""	,	radius	=	0.9	,	edges	=	200,	clockwise	=	T, col = mycol2)																																			
add.pie(pies  $  COP[13]	,	x	=	9.5	,	y	=	54.8	,	labels	=	""	,	radius	=	1.05	,	edges	=	200,	clockwise	=	T, col = mycol2)																																			

add.pie(pies  $  GAM[1:11]  ,  x	=	9.5	,	y	=	56.5	,	labels	=	""	,	radius	=	1.05	,	edges	=	200,	clockwise	=	T, col = mycol)																																			
add.pie(pies  $  COP[1:11]  ,	x	=	9.5	,	y	=	54.8	,	labels	=	""	,	radius	=	1.05	,	edges	=	200,	clockwise	=	T, col = mycol)																																			



## S.Swedish pie ###

add.pie(pies  $	SK[13]	,	x	=	13.152523	,	y	=	55.550972	,	labels	=	""	,	radius	=	0.959416667	,	edges	=	200,	clockwise	=	T, col = mycol2)																																			


######## Pool 3 ###########

## Polish Pies ##

add.pie(pies$GR2[13],x=15.357225,y=52.9312583,labels="",radius=1.265833334,edges=200,clockwise=T, col = mycol3)

add.pie(pies  $	GR1[13]	,	x	=	19.5	,	y	=	54.897537	,	labels	=	""	,	radius	=	1.34925	,	edges	=	200,	clockwise	=	T, col = mycol3)																																			

add.pie(pies	$	TU[13]	,	x	=	20.5	,	y	=	50.5	,	labels	=	""	,	radius	=	1.477666667	,	edges	=	200,	clockwise	=	T, col = mycol3)																																			

add.pie(pies$POLEN[13],x=25.022095,y=53,labels="",radius=1.134958334,edges=200,clockwise=T, col = mycol3)


##### Pool 4 #######
##### Central Sweden ######


add.pie(pies$STYV[13],x=14.271862,y=57.561081,labels="",radius=0.870625,edges=200,clockwise=T, col = mycol4)
add.pie(pies$EK[13],x=13,y=61,labels="",radius=1.012708334,edges=200,clockwise=T, col = mycol4)
add.pie(pies$SD[13],x=12.5,y=63,labels="",radius=0.903333334,edges=200,clockwise=T, col = mycol4)
add.pie(pies$GD[13],x=15.5,y=64.5,labels="",radius=1.195166667,edges=200,clockwise=T, col = mycol4)
add.pie(pies$OST[13],x=18.380814,y=61.8,labels="",radius=1.171291667,edges=200,clockwise=T, col = mycol4)
add.pie(pies$LMCA[13],x=9.5,y=59.453506,labels="",radius=1.28,edges=200,clockwise=T, col = mycol4)


add.pie(pies  $	KAP[13]	,	x	=	18.785334	,	y	=	57.849045	,	labels	=	""	,	radius	=	0.7	,	edges	=	200,	clockwise	=	T, col = mycol4)																																			
add.pie(pies  $  WEN[13]  ,	x	=	18.95	,	y	=	59.66 ,	labels	=	""	,	radius	=	1.0	,	edges	=	200,	clockwise	=	T, col = mycol4)																																			
add.pie(pies  $  OBY[13]	,	x	=	17.79	,	y	=	60.21	,	labels	=	""	,	radius	=	0.758	,	edges	=	200,	clockwise	=	T, col = mycol4)																																			
add.pie(pies  $  AL[13]	,	x	=	19.852339	,	y	=	60.359329	,	labels	=	""	,	radius	=	1.239541667	,	edges	=	200,	clockwise	=	T, col = mycol4)																																			


##### Northern Sweden  ######
add.pie(pies  $  UMCA[13]	,	x	=	20.405216	,	y	=	63.712364	,	labels	=	""	,	radius	=	0.945041667	,	edges	=	200,	clockwise	=	T, col = mycol5)



##### N. Finnish pies #####

add.pie(pies  $	CALK[13]	,	x	=	25.758348	,	y	=	62.262291	,	labels	=	""	,	radius	=	0.713958334	,	edges	=	200,	clockwise	=	T, col = mycol6)																																			
add.pie(pies  $	NLP[13]	,	x	=	29.676218	,	y	=	62.680271	,	labels	=	""	,	radius	=	0.630083334	,	edges	=	200,	clockwise	=	T, col = mycol6)																																			
add.pie(pies	$	OU[13]	,	x	=	25.472832	,	y	=	65.012375	,	labels	=	""	,	radius	=	1.046208334	,	edges	=	200,	clockwise	=	T, col = mycol6)																																			

#### Tromso ######

add.pie(pies  $  TROM[13]	,	x	=	18.95	,	y	=	69.65	,	labels	=	""	,	radius	=	0.59	,	edges	=	200,	clockwise	=	T, col = mycol7)																																			


#### N.E. EU #### 

add.pie(pies	$	SA[13]	,	x	=	23.104763	,	y	=	60.371616	,	labels	=	""	,	radius	=	1.102416667	,	edges	=	200,	clockwise	=	T, col = mycol8)																																			
add.pie(pies$VIIKCA[13],x=29,y=60.363937,labels="",radius=1.102666667,edges=200,clockwise=T, col = mycol8)
add.pie(pies	$	MY20[13]	,	x	=	30.248108	,	y	=	55.903034	,	labels	=	""	,	radius	=	1.1505	,	edges	=	200,	clockwise	=	T, col = mycol8)																																			
add.pie(pies  $	EST[13]	,	x	=	24.3	,	y	=	57.8	,	labels	=	""	,	radius	=	0.7	,	edges	=	200,	clockwise	=	T, col = mycol8)																																			
add.pie(pies$EST2[13],x=29.5,y=58,labels="",radius=0.7,edges=200,clockwise=T, col = mycol8)
add.pie(pies  $	UKR[13]	,	x	=	30.52002	,	y	=	52.469398	,	labels	=	""	,	radius	=	1.448208334	,	edges	=	200,	clockwise	=	T, col = mycol8)																																			



## Danubian Pies ##

add.pie(pies$GEW3[13],x=8.502993,y=46.5,labels="",radius=1.172458334,edges=200,clockwise=T, col = mycol9)

add.pie(pies$GEW8[13],x=13.02993,y=46.5,labels="",radius=1.069291667,edges=200,clockwise=T, col = mycol9)

add.pie(pies$CA.CR[13],x=17.5,y=47.8764583,labels="",radius=1.205708334,edges=200,clockwise=T, col = mycol9)


## Lower Europe Pies ##

add.pie(pies	$	PRO[13]	,	x	=	40.46814	,	y	=	47.457809	,	labels	=	""	,	radius	=	1.279916667	,	edges	=	200,	clockwise	=	T, col = mycol10)																																			

## New pies ## NEED TO DO RADIUSES!!




## Finally, need to add a Key ###

keypie1= c(0,1)

add.pie(keypie1,x=-7.5,y=45,labels="",radius=1.5, edges=200,clockwise=T, col = "black")
add.pie(keypie1,x=2,y=45,labels="",radius=.5, edges=200,clockwise=T, col = "black")
add.pie(keypie1,x=-2,y=45,labels="",radius=1., edges=200,clockwise=T, col = "black")
