
source("Prey.models.R")
source("jaw.model.fun.R")
source("~/Documents/R.scripts/ggplot/bw.theme.r")

#config=1 is  no sling Aw mass added to A2 mass, 2 is no sling and Aw mass not added, 3 is sling
run.jaw.prey.2 <- function(spec.n=3,dat=NULL,config=1,Loose=F,min.ML=0.6,OutPutVerb=T,prey=F,press=100,prey.per=NULL,prey.strike.ang=NULL,prey.pos="flat",MaxIts=1000,progress=F,print.every=10,inset=T,out="pdf")	{
#*********** File Control ***************

#remove pdfs in graphics dir


#where graphics will be stored
  
graph.dir <- file.path(getwd(),"graphics",paste(config,spec.n,sep="."))
dir.create(graph.dir,showWarnings=F)
	
#Load specimen data and parameter values
MuscleData<-read.csv(dat)
ParameterValues<- read.csv(file=paste(getwd(),"/ParValPrey.csv",sep=""), head=TRUE)
#print(MuscleData)
#specimen number to index with

SpecimenNumber=spec.n;
#Set morphological variables
SL=MuscleData[SpecimenNumber,"SL"];
JL=MuscleData[SpecimenNumber,"JL"];

#A2 geometry and mass
A2OrIns=MuscleData[SpecimenNumber,"A2OrIns"];
A2OrJoint=MuscleData[SpecimenNumber,"A2OrJoint"]; #change w/ SpecAnge
A2Tendon=MuscleData[SpecimenNumber,"A2Tendon"]
A2Li=MuscleData[SpecimenNumber,"A2Li"];
A2Lo=MuscleData[SpecimenNumber,"A2Lo"];
A2AvgFL=MuscleData[SpecimenNumber,"A2AvgFL"];#change w/ SpecAnge
A2AvgPenn=(pi/180)*MuscleData[SpecimenNumber,"A2AvgPenn"]; #change w/ SpecAng
A2Mass=ifelse(config==1,MuscleData[SpecimenNumber,"A2Mass"]+MuscleData[SpecimenNumber,"AwMass"],MuscleData[SpecimenNumber,"A2Mass"]);

#A3 geometry and mass
A3OrIns =MuscleData[SpecimenNumber,"A3OrIns"]
A3OrJoint = MuscleData[SpecimenNumber,"A3OrJoint"]
A3Tendon= MuscleData[SpecimenNumber,"A3Tendon"]
A3Li =MuscleData[SpecimenNumber,"A3Li"]
A3Lo <- MuscleData[SpecimenNumber,"A3Lo"]
A3AvgFL <- MuscleData[SpecimenNumber,"A3AvgFL"]
A3AvgPenn <- (pi/180)*MuscleData[SpecimenNumber,"A3AvgPenn"]
A3Mass <- MuscleData[SpecimenNumber,"A3Mass"]

#Aw geometry and mass
AwOrIns <- MuscleData[SpecimenNumber,"AwOrIns"]
AwOrJoint <- MuscleData[SpecimenNumber,"AwOrJoint"]
AwTendon <- MuscleData[SpecimenNumber,"AwTendon"]
AwLi <- MuscleData[SpecimenNumber,"AwLi"]
AwLo <- MuscleData[SpecimenNumber,"AwLo"]
AwAvgFL <- MuscleData[SpecimenNumber,"AwAvgFL"]
AwAvgPenn <- (pi/180)*MuscleData[SpecimenNumber,"AwAvgPenn"]
AwMass <- ifelse(config==1|config==2,0,MuscleData[SpecimenNumber,"AwMass"])

#jaw geometry, gape
JawWidth <- MuscleData[SpecimenNumber,"Sspine"]
MandWidth <- MuscleData[SpecimenNumber,"MandibleWidth"]
MandDepth <- MuscleData[SpecimenNumber,"MandibleDepth"]
SymphLength <- MuscleData[SpecimenNumber,"SymphysisLength"]
Gape <- (pi/180)*MuscleData[SpecimenNumber,"Gape"]
SpecAng<- (pi/180)*MuscleData[SpecimenNumber,"SpecAng"]
Pressure <-press #MuscleData[SpecimenNumber,"Pressure"]
print(press)

#species, passed to outfile name
Species <-MuscleData[SpecimenNumber,"Species"]

#Where results will be printed (either verbose or tidy, prey or no prey)
OutFileName <- paste(getwd(),"/",MuscleData$Species[spec.n],prey.pos,"CloseOutFile.csv",sep="")

#Set model parameter input variables
TimeStep <- ParameterValues$Value[1]
TimeStepSecs <-TimeStep/1000
MuscleDen <- ParameterValues$Value[2]
MaxIso <- ParameterValues$Value[3]
MaxFpar <- ParameterValues$Value[4]
Vmax <- ParameterValues$Value[5]
#shape of Hill curve
G <- ParameterValues$Value[6]
#Activation tise time
ActRiseTime <- ParameterValues$Value[7] 
MLOpt <- ParameterValues$Value[8]/100;

#fiber length-tension relationship from Porro et al. (2011)
p.fl <- c(0.508, 0.68, 0.92, 1.0, 1.56) #fl/flo
p.f <- c(0, 0.84, 1, 1, 0) #%force
fl.poly <- lm(p.f ~ poly(p.fl, 4, raw=TRUE))



ThetaClosed <- (pi/180)*ParameterValues$Value[9];
Scale <- ParameterValues$Value[10]
CdPlate <- ParameterValues$Value[11]
Start <- TimeStep*-1;

#prey parameters
PreyProp <- prey.per;
PreyStrikeDist <- 0.85*A2Lo#MuscleData[SpecimenNumber,"strike"];
PreyStrikeAng <-prey.strike.ang*(pi/180);

#Set type of Simulation
Loosejaw <- Loose;
ThinFlat <- prey.pos;

#Calculate Fmax of each muscle w/ out cos(pennationangle)
A3FmaxNoPenn <-(A3Mass)/(MuscleDen*A3AvgFL)*2*MaxIso/10;
A2FmaxNoPenn <-(A2Mass)/(MuscleDen*A2AvgFL)*2*MaxIso/10;
AwFmaxNoPenn <-(AwMass)/(MuscleDen*AwAvgFL)*2*MaxIso/10;

#Calculate PCSA closed
A3PCSAclosed <-(A3Mass*cos(A3AvgPenn))/(MuscleDen*A3AvgFL)*2
A2PCSAclosed <-(A2Mass*cos(A2AvgPenn))/(MuscleDen*A2AvgFL)*2
AwPCSAclosed <-(AwMass*cos(AwAvgPenn))/(MuscleDen*AwAvgFL)*2

#Calculate MA
A3MA <- A3Li/A3Lo
A2MA <- A2Li/A2Lo
AwMA <- AwLi/AwLo
# 
# #Calculate Fmax (N) in closed position
A3FmaxClosed <-MaxIso/10*A3PCSAclosed
A2FmaxClosed <-MaxIso/10*A2PCSAclosed
AwFmaxClosed <-MaxIso/10*AwPCSAclosed

#*************** Muscle geometry *******************

#Angle of line between origin-joint and lower jaw when closed
A3ThetaJointClosed <- acos((A3Li^2+A3OrJoint^2-A3OrIns^2)/(2*A3Li*A3OrJoint))
A2ThetaJointClosed <- acos((A2Li^2+A2OrJoint^2-A2OrIns^2)/(2*A2Li*A2OrJoint))
AwThetaJointClosed <- acos((AwLi^2+AwOrJoint^2-AwOrIns^2)/(2*AwLi*AwOrJoint))

#Angle of line between origin-joint and lower jaw when open
A3ThetaJointOpen <- A3ThetaJointClosed+Gape-SpecAng
A2ThetaJointOpen <- A2ThetaJointClosed+Gape-SpecAng
AwThetaJointOpen <- AwThetaJointClosed+Gape-SpecAng

#Angle of insertion on line of actions when jaw is closed
A3MuscleThetaClosed <-acos((A3Li^2+A3OrIns^2-A3OrJoint^2)/(2*A3Li*A3OrIns))
A2MuscleThetaClosed <-acos((A2Li^2+A2OrIns^2-A2OrJoint^2)/(2*A2Li*A2OrIns))
AwMuscleThetaClosed <-acos((AwLi^2+AwOrIns^2-AwOrJoint^2)/(2*AwLi*AwOrIns))


########## cartesian coordinates closed #######
joint.p <- c(0,0)
start <- rad(-80)
inlever.p <- point.ang.r(joint.p,A2Li,-start)
JawTip <- point.ang.r(joint.p,JL,-start)
A2Or <- point.ang.r(joint.p,A2OrJoint,A2ThetaJointClosed-start)
A3Or <- point.ang.r(joint.p,A3OrJoint,A3ThetaJointClosed-start)
#figure our Aw angles
AwOr <- point.ang.r(joint.p,AwLi,-start)
#insertion point on tendon based on muscle theta and orig-ins length (so distal tendon)
TendonLengthClosed <- sin(AwMuscleThetaClosed)*AwOrIns

#force Aw to insert on middle of tendon according to muscle theta
AwIns<- point.ang.r(AwOr,AwLi-A2Li,-rad(180)-AwMuscleThetaClosed-start)

#geom.df <- t(data.frame(joint.p=joint.p,A2Or=A2Or,A3Or=A3Or,JawTip=JawTip,inlever.p=inlever.p,AwIns=AwOr,nexus=AwIns));colnames(geom.df) <- c("x","y")
#geom.df.melt <- as.data.frame(geom.df)
#p <- ggplot(data=geom.df.melt,aes(x=x,y=y))+geom_point(aes(colour=x),size=6)+bw.theme(begin=-2.5,stop=2.5,x.lab="x",y.lab="y",axis.p=c(0.1,0.8),leg.text=15)+ylim(-2.5,2.5)+line.2pt("joint.p","JawTip",geom.df.melt)+line.2pt("nexus","AwIns",geom.df.melt,col="darkgreen")+line.2pt("A3Or","nexus",geom.df.melt,col="brown")+line.2pt("A2Or","nexus",geom.df.melt,col="red")+line.2pt("inlever.p","nexus",geom.df.melt,col="gray")

#png(filename=paste(graph.dir,"/geom_closed.j.png",sep=""),width=700,height=700)
#print(p)
#dev.off()

################################################################

#Length of Muscle (origin to insertion on tendon + proximal tendon) when closed, based on cartesian coords
A3MLClosed <- dist.2d(AwIns,A3Or)
A2MLClosed<- dist.2d(AwIns,A2Or)
AwMLClosed<- dist.2d(AwIns,AwOr)


#Length along line of action when open

A3OrInsOpen <-(A3Li^2+A3OrJoint^2-cos(A3ThetaJointClosed+Gape-SpecAng)*A3OrJoint*A3Li*2)^(0.5)
A2OrInsOpen <-(A2Li^2+A2OrJoint^2-cos(A2ThetaJointClosed+Gape-SpecAng)*A2OrJoint*A2Li*2)^(0.5)
AwOrInsOpen <-(AwLi^2+AwOrJoint^2-cos(AwThetaJointClosed+Gape-SpecAng)*AwOrJoint*AwLi*2)^(0.5)

#Angle of along line of action relative to lowerjaw  when jaw is open
A3MuscleThetaOpen <-acos((A3Li^2+A3OrInsOpen^2-A3OrJoint^2)/(2*A3Li*A3OrInsOpen))
A2MuscleThetaOpen <-acos((A2Li^2+A2OrInsOpen^2-A2OrJoint^2)/(2*A2Li*A2OrInsOpen))

########## cartesian coordinates open #######
open.ang <- -1*(start+Gape-SpecAng)
inlever.p.open <- point.ang.r(joint.p,A2Li,open.ang)
JawTip.open <- point.ang.r(joint.p,JL,open.ang)
#figure out Aw angles
AwOr.open <- point.ang.r(joint.p,AwLi,open.ang)

#force Aw to insert on middle of tendon according to muscle theta
#mean of two muscle insertions scaled around 180
TendonThetaOpen <-mean(c(A2MuscleThetaOpen,A3MuscleThetaOpen))
AwIns.open<- point.ang.r(inlever.p.open,TendonLengthClosed,(rad(180)-(Gape-SpecAng))-TendonThetaOpen-start) #really the origin and the nexus

nexus <- AwIns.open

#where A3,A2, Aw meet at aponeurosis

A3MLOpen<- dist.2d(AwIns.open,A3Or)
A2MLOpen<- dist.2d(AwIns.open,A2Or)
AwMLOpen<- dist.2d(AwIns.open,AwOr)

A3TendonThetaOpen <- TendonThetaOpen-A3MuscleThetaOpen
A2TendonThetaOpen <- A2MuscleThetaOpen-TendonThetaOpen

#angle between MLs and Aw MLs
A3AwThetaOpen <- cos.ang(AwMLOpen,A3MLOpen,dist.2d(AwOr.open,A3Or))
A2AwThetaOpen <- cos.ang(A2MLOpen,AwMLOpen,dist.2d(AwOr.open,A2Or))


	
#get signs right
#geom.df.open <- t(data.frame(joint.p,A2Or,A3Or,JawTip.open,inlever.p.open,AwIns=AwOr.open,nexus=AwIns.open));colnames(geom.df.open) <- c("x","y")
#geom.df.open.melt <- as.data.frame(geom.df.open)
#p.open <- ggplot(data=geom.df.open.melt,aes(x=x,y=y))+geom_point(aes(colour=x),size=6)+bw.theme(begin=-2.5,stop=2.5,x.lab="x",y.lab="y",axis.p=c(0.1,0.8),leg.text=25)+ylim(-2.5,2.5)+line.2pt("joint.p","JawTip.open",geom.df.open.melt)+theme(legend.position="none")+line.2pt("nexus","AwIns",geom.df.open.melt,col="darkgreen")+line.2pt("A3Or","nexus",geom.df.open.melt,col="brown")+line.2pt("A2Or","nexus",geom.df.open.melt,col="red")+line.2pt("inlever.p.open","nexus",geom.df.open.melt,col="gray")


#png(filename=paste(graph.dir,"/geom_open.j.png",sep=""),width=700,height=700)
#print(p)
#dev.off()
##############


#Aw insertion angle calculated from open geometry
AwMuscleThetaOpen <- asin(sin(pi-TendonThetaOpen)/AwMLOpen*TendonLengthClosed)
#Aw Origin Angle calculated from ^ and subtraction
AwTendonThetaOpen <-pi-(pi-TendonThetaOpen+AwMuscleThetaOpen)

#Total muscle contraction
A3Con <-A3MLOpen-A3MLClosed
A2Con <-A2MLOpen-A2MLClosed
AwCon <-AwMLOpen-AwMLClosed

#Fiber length when jaw is open
A3FLOpen <-sqrt((sin(A3AvgPenn)*A3AvgFL)^2+((cos(A3AvgPenn)*A3AvgFL)+A3Con)^2)
A2FLOpen <-sqrt((sin(A2AvgPenn)*A2AvgFL)^2+((cos(A2AvgPenn)*A2AvgFL)+A2Con)^2)
AwFLOpen <-sqrt((sin(AwAvgPenn)*AwAvgFL)^2+((cos(AwAvgPenn)*AwAvgFL)+AwCon)^2)

#Pennation angle open
A3PennOpen <-asin((sin(A3AvgPenn)*A3AvgFL)/A3FLOpen)
A2PennOpen <-asin((sin(A2AvgPenn)*A2AvgFL)/A2FLOpen)
AwPennOpen <-asin((sin(AwAvgPenn)*AwAvgFL)/AwFLOpen)

#Optimal muscle length (check if these makes sense 21Aug14)
A3MLOpt <-A3MLOpen-(A3MLOpen-A3MLClosed)*MLOpt
A2MLOpt <- A2MLOpen-(A2MLOpen-A2MLClosed)*MLOpt
AwMLOpt <- AwMLOpen-(AwMLOpen-AwMLClosed)*MLOpt
# A3MLOpt <-A3MLOpen*MLOpt
# A2MLOpt <- A2MLOpen*MLOpt
# AwMLOpt <- AwMLOpen*MLOpt

#Optimal fiber length
A3FLOpt <- sqrt((sin(A3AvgPenn)*A3FLOpen)^2+(cos(A3AvgPenn)*A3FLOpen-(A3MLOpen-A3MLOpt))^2)
A2FLOpt <-  sqrt((sin(A2AvgPenn)*A2FLOpen)^2+(cos(A2AvgPenn)*A2FLOpen-(A2MLOpen-A2MLOpt))^2)
AwFLOpt <-  sqrt((sin(AwAvgPenn)*AwFLOpen)^2+(cos(AwAvgPenn)*AwFLOpen-(AwMLOpen-AwMLOpt))^2)

#set initial Fmaxs to 0
A3Fmax=0;A2Fmax=0;AwFmax=0
#**************** Lower-jaw dimensions for loosejaws****************

a <- sqrt(A3Lo^2-(JawWidth/2)^2)
bprime <-SymphLength/a*(JawWidth/2)
s <-sqrt(SymphLength^2+bprime^2)
l <-a-SymphLength
MandRadius <-(MandWidth+MandDepth)/4

#******************** Mass components *******************

#Mass of rami minus symphyseal membrane
MandMass <-(MandRadius/Scale)^2*pi*((A3Lo/Scale)-(s/Scale))*1000

#Cd of rami ellipsoids after Blevins (1994)
CdRami <-0.966*(MandDepth/MandWidth)^(-0.746)

#Moment of inertia of half ellipse
Inormal <-(2/15)*pi*1000*(a/Scale)^3*((JawWidth/2)/Scale)^2

#Moment of interia of partial ellipse (loosejaw)
Ipartellip <-(((JawWidth/2)/Scale)^2*sqrt((a/Scale)^2-(l/Scale)^2)*(2*(a/Scale)^4+(a/Scale)^2*(l/Scale)^2-3*(l/Scale)^4)*pi*1000)/(15*(a/Scale)^2)

#Moment of inertial of rami (loosejaw)
Irami <-(2*(0.25*MandMass*(MandRadius/Scale)^2+(1/3)*MandMass*((A3Lo-s)/Scale)^2))

#Total moment of inertia (loosejaw)
Iloose <-Ipartellip+Irami

#Moment of inertia of distal tendon, distal to muscle insertions (I=(m*L^2)/3), a thin rigid wire, scale to mass of A3 which is about as thick
Itendon <- ((TendonLengthClosed/A3MLClosed)*A3Mass/1000*(TendonLengthClosed/Scale)^2)/3


#Mass and dimensional components of prey
PreySL <- PreyProp*SL
PreyDatFile <- paste(getwd(),"/","PreyDat.csv",sep="")
PreyDat <- read.csv(PreyDatFile)
#models are in mm SL,returning mass in mg, width in cm
PreyMass <- DiaphusMass.fun(PreySL*10,mass=PreyDat$mass,SL=PreyDat$SL)
PreyWidth<- DiaphusWidth.fun(PreySL,width=PreyDat$width,SL=PreyDat$SL)
PreyDepth<- DiaphusDepth.fun(PreySL,depth=PreyDat$depth,SL=PreyDat$SL)

Prey <- prey


#mass moment inertia of prey; /scale removed from PreyStrikeDist 25 June
Iprey <- (PreyMass/1000)*0.25*((PreyWidth/Scale)^2+(PreyDepth/Scale)^2)+((PreyMass/1000)*PreyStrikeDist/Scale^2)


#Drag of prey components

PreyCoM <- PreyStrikeDist #Prey center of Mass
CdThin <- (0.966*(PreyDepth/PreyWidth)^(-0.746)) #Cd in thin position
CdFlat  <- (0.966*(PreyWidth/PreyDepth)^(-0.746)) #Cd in flat position
#review drag term integrals, perhaps they make more sense
# FlatDragTerm <- (PreyDepth/Scale*PreySL/Scale)*CdFlat
# ThinDragTerm <- (PreyWidth/Scale*PreySL/Scale)*CdThin
FlatDistD <- PreyCoM+PreyDepth/2 #Distance (radius) to distal margin of prey in flat position
FlatDistP <- PreyCoM-PreyDepth/2 #Distance (radius) to proximal margin of prey in flat position
ThinDistD <- PreyCoM+PreyWidth/2 #Distance (radius) to distal margin of prey in thin position
ThinDistP <- PreyCoM-PreyWidth/2 #Distance (radius) to proxima margin of prey in thin position ;
FlatDragTerm <- ((FlatDistD/Scale)^4-(FlatDistP/Scale)^4)*CdFlat
ThinDragTerm <- ((ThinDistD/Scale)^4-(ThinDistP/Scale)^4)*CdThin
if(ThinFlat=="flat"){DragTerm <- FlatDragTerm}else{if(ThinFlat=="thin"){DragTerm <- ThinDragTerm}else{stop("ThinOrFlat input unrecognized. Only 'thin' or 'flat' accepted.")}}


#Which jaw mass to use?
if(Loosejaw==T){Ijaw <- Iloose}else{if(Loosejaw==F){Ijaw <- Inormal}else{stop("Loose or Normal simulation input unrecognized. Only 'T' or 'F' accepted.")}}

if(Prey==T){Itot <- Ijaw+Iprey}else{if(Prey==F){Itot <- Ijaw}else{stop("Prey input unrcognized. Only 'T' or 'F' accepted.")}}

#plot start



					#********* MODEL CALCULATIONS *************

#****************** ITERATIONS ************

	
	j.var <- data.frame(Iteration=1,source("var.df.R",local=T)$value)
	j <- data.frame(Iteration=c(2,MaxIts),matrix(0,ncol=ncol(j.var),nrow=MaxIts))
	j[1,] <- j.var
	colnames(j) <- colnames(j.var)
	
	geom.out <- list()

	pb = txtProgressBar(min = 2, max = MaxIts, initial = 2,style=3)
	
	for(n in 2:MaxIts){ifelse(j[n-1,"JawAng"]>ThetaClosed,source("jaw.torque.MS4.R",local=T),n=MaxIts)	
	  setTxtProgressBar(pb,n)
		}


geom.dat <-data.frame(A3FmaxNoPenn,A2FmaxNoPenn,AwFmaxNoPenn,A3PCSAclosed,A2PCSAclosed,AwPCSAclosed,A3FmaxClosed,A2FmaxClosed,AwFmaxClosed,A3ThetaJointClosed,A2ThetaJointClosed,AwThetaJointClosed,A3ThetaJointOpen,A2ThetaJointOpen,AwThetaJointOpen,A3MuscleThetaClosed,A2MuscleThetaClosed,AwMuscleThetaClosed,A3MLOpen,A2MLOpen,AwMLOpen,A3MuscleThetaOpen,A2MuscleThetaOpen,AwMuscleThetaOpen,A3Con,A2Con,AwCon,A3FLOpen,A2FLOpen,AwFLOpen,A3PennOpen,A2PennOpen,AwPennOpen,A3MLOpt,A2MLOpt,AwMLOpt,A3FLOpt,A2FLOpt,AwFLOpt)

res.list <- list(out.dat=j[2:which.max(j$Time),],time=max(j$Time),geom.dat=geom.dat,stat.bite=max(j$StaticBite),max.angvel=max(j$AngVel),max.jawvel.tip=max(j$JawVelTip),geom.out=geom.out)

return(res.list)


}




