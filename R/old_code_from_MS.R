
# from git and MS. Produce results from paper. Feb 23

if(F){
start_old <- function(x=NULL,spec.n=3,config=1,pars=NULL)	{
  #*********** File Control ***************
  if(!is.null(x)){ w <- tryCatch(review_data(x),warning = function(w) w)

  if("warning" %in% class(w)) stop("Input data passed to `x` are broken. They may have missing variables, non-numeric data, etc. Please use `reveiew_data()`.")
  }
  if(!is.null(x) & !is.data.frame(x)) stop("'x' and data passed to `start_sim()` must be a data.frame")

  #Load specimen data and parameter values
  if(!is.null(x)) MuscleData<-x


  if(is.null(x)) MuscleData<-system.file("extdata", "csloani.csv", package = "JawCloseR") %>% read.csv() %>% rownames_to_column(var="spec.n")

  ParameterValues <-system.file("extdata", "par_val.csv", package = "JawCloseR") %>% read.csv()



  if(!is.null(pars)){  if (any(!names(pars) %in% ParameterValues$Parameter))
    stop(
      '`pars` must be either NULL or a named list with names in c("TimeStepMS","Density","MaxIso","MaxFpar","Vmax","Gcon","ActivTime","MLOptimum","GapeClosed","Scale","CdPlate","MLmin","Pressure")'
    )
  }

  if (!is.null(pars))
    for (p in names(pars))
      ParameterValues[ParameterValues$Parameter == p, ]$Value <- pars[[which(names(pars) == p)]]


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
  JawWidth <- MuscleData[SpecimenNumber,"GapeWidth"]
  MandWidth <- MuscleData[SpecimenNumber,"MandWidth"]
  MandDepth <- MuscleData[SpecimenNumber,"MandDepth"]
  SymphLength <- MuscleData[SpecimenNumber,"SymphLength"]
  Gape <- (pi/180)*MuscleData[SpecimenNumber,"Gape"]
  SpecAng<- (pi/180)*MuscleData[SpecimenNumber,"SpecAng"]


  #species, passed to outfile name
  Species <-MuscleData[SpecimenNumber,"Species"]


  get_var <- function(x,v){
    x %>% filter(Parameter==v) %>% pull(Value)
  }


  #Set model parameter input variables
  TimeStep <- ParameterValues %>% get_var(v="TimeStepMS")
  TimeStepSecs <- TimeStep / 1000
  MuscleDen <- ParameterValues %>% get_var(v="Density")
  MaxIso <- ParameterValues %>% get_var(v="MaxIso")
  MaxFpar <- ParameterValues %>% get_var(v="MaxFpar")
  Vmax <- ParameterValues %>% get_var(v="Vmax")
  #shape of Hill curve
  G <- ParameterValues %>% get_var(v="Gcon")
  #Activation tise time
  ActRiseTime <- ParameterValues %>% get_var(v="ActivTime")
  MLOpt <- ParameterValues %>% get_var(v="MLOptimum")/100

  Pressure <-ParameterValues %>% get_var(v="Pressure")


  ThetaClosed <-ParameterValues %>% get_var(v="GapeClosed") %>% rad
  Scale <-  ParameterValues %>% get_var(v="Scale")
  CdPlate <-  ParameterValues %>% get_var(v="CdPlate")
  Start <- TimeStep * -1


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
  inlever.p <- point_ang_r(joint.p,A2Li,-start)
  JawTip <- point_ang_r(joint.p,JL,-start)
  A2Or <- point_ang_r(joint.p,A2OrJoint,A2ThetaJointClosed-start)
  A3Or <- point_ang_r(joint.p,A3OrJoint,A3ThetaJointClosed-start)
  #figure our Aw angles
  AwOr <- point_ang_r(joint.p,AwLi,-start)
  #insertion point on tendon based on muscle theta and orig-ins length (so distal tendon)
  TendonLengthClosed <- sin(AwMuscleThetaClosed)*AwOrIns

  #force Aw to insert on middle of tendon according to muscle theta
  AwIns<- point_ang_r(AwOr,AwLi-A2Li,-rad(180)-AwMuscleThetaClosed-start)

  #geom.df <- t(data.frame(joint.p=joint.p,A2Or=A2Or,A3Or=A3Or,JawTip=JawTip,inlever.p=inlever.p,AwIns=AwOr,nexus=AwIns));colnames(geom.df) <- c("x","y")
  #geom.df.melt <- as.data.frame(geom.df)
  #p <- ggplot(data=geom.df.melt,aes(x=x,y=y))+geom_point(aes(colour=x),size=6)+bw.theme(begin=-2.5,stop=2.5,x.lab="x",y.lab="y",axis.p=c(0.1,0.8),leg.text=15)+ylim(-2.5,2.5)+line.2pt("joint.p","JawTip",geom.df.melt)+line.2pt("nexus","AwIns",geom.df.melt,col="darkgreen")+line.2pt("A3Or","nexus",geom.df.melt,col="brown")+line.2pt("A2Or","nexus",geom.df.melt,col="red")+line.2pt("inlever.p","nexus",geom.df.melt,col="gray")

  #png(filename=paste(graph.dir,"/geom_closed.j.png",sep=""),width=700,height=700)
  #print(p)
  #dev.off()

  ################################################################

  #Length of Muscle (origin to insertion on tendon + proximal tendon) when closed, based on cartesian coords
  A3MLClosed <- dist_2d(AwIns,A3Or)
  A2MLClosed<- dist_2d(AwIns,A2Or)
  AwMLClosed<- dist_2d(AwIns,AwOr)


  #Length along line of action when open

  A3OrInsOpen <-(A3Li^2+A3OrJoint^2-cos(A3ThetaJointClosed+Gape-SpecAng)*A3OrJoint*A3Li*2)^(0.5)
  A2OrInsOpen <-(A2Li^2+A2OrJoint^2-cos(A2ThetaJointClosed+Gape-SpecAng)*A2OrJoint*A2Li*2)^(0.5)
  AwOrInsOpen <-(AwLi^2+AwOrJoint^2-cos(AwThetaJointClosed+Gape-SpecAng)*AwOrJoint*AwLi*2)^(0.5)

  #Angle of along line of action relative to lowerjaw  when jaw is open
  A3MuscleThetaOpen <-acos((A3Li^2+A3OrInsOpen^2-A3OrJoint^2)/(2*A3Li*A3OrInsOpen))
  A2MuscleThetaOpen <-acos((A2Li^2+A2OrInsOpen^2-A2OrJoint^2)/(2*A2Li*A2OrInsOpen))

  ########## cartesian coordinates open #######
  open.ang <- -1*(start+Gape-SpecAng)
  inlever.p.open <- point_ang_r(joint.p,A2Li,open.ang)
  JawTip.open <- point_ang_r(joint.p,JL,open.ang)
  #figure out Aw angles
  AwOr.open <- point_ang_r(joint.p,AwLi,open.ang)

  #force Aw to insert on middle of tendon according to muscle theta
  #mean of two muscle insertions scaled around 180
  TendonThetaOpen <-mean(c(A2MuscleThetaOpen,A3MuscleThetaOpen))
  AwIns.open<- point_ang_r(inlever.p.open,TendonLengthClosed,(rad(180)-(Gape-SpecAng))-TendonThetaOpen-start) #really the origin and the nexus

  nexus <- AwIns.open

  #where A3,A2, Aw meet at aponeurosis

  A3MLOpen<- dist_2d(AwIns.open,A3Or)
  A2MLOpen<- dist_2d(AwIns.open,A2Or)
  AwMLOpen<- dist_2d(AwIns.open,AwOr)

  A3TendonThetaOpen <- TendonThetaOpen-A3MuscleThetaOpen
  A2TendonThetaOpen <- A2MuscleThetaOpen-TendonThetaOpen

  #angle between MLs and Aw MLs
  A3AwThetaOpen <- cos_ang(AwMLOpen,A3MLOpen,dist_2d(AwOr.open,A3Or))
  A2AwThetaOpen <- cos_ang(A2MLOpen,AwMLOpen,dist_2d(AwOr.open,A2Or))



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

  #Optimal muscle length
  A3MLOpt <-A3MLOpen-(A3MLOpen-A3MLClosed)*MLOpt
  A2MLOpt <- A2MLOpen-(A2MLOpen-A2MLClosed)*MLOpt
  AwMLOpt <- AwMLOpen-(AwMLOpen-AwMLClosed)*MLOpt


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

  if(F)
  {  #Mass and dimensional components of prey
    # PreySL <- PreyProp*SL
    # PreyDatFile <- paste(getwd(),"/","PreyDat.csv",sep="")
    # PreyDat <- read.csv(PreyDatFile)
    # #models are in mm SL,returning mass in mg, width in cm
    # PreyMass <- DiaphusMass.fun(PreySL*10,mass=PreyDat$mass,SL=PreyDat$SL)
    # PreyWidth<- DiaphusWidth.fun(PreySL,width=PreyDat$width,SL=PreyDat$SL)
    # PreyDepth<- DiaphusDepth.fun(PreySL,depth=PreyDat$depth,SL=PreyDat$SL)
    #
    Prey <- prey <- 0


    #mass moment inertia of prey; /scale removed from PreyStrikeDist 25 June
    # Iprey <- (PreyMass/1000)*0.25*((PreyWidth/Scale)^2+(PreyDepth/Scale)^2)+((PreyMass/1000)*PreyStrikeDist/Scale^2)
    #

    #Drag of prey components

    # PreyCoM <- PreyStrikeDist #Prey center of Mass
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

  }

  Itot <- Inormal
  #plot start

  start.l <-
    list(
      a=a,
      A2AvgFL = A2AvgFL,
      A2AvgPenn = A2AvgPenn,
      A2AwThetaOpen = A2AwThetaOpen,
      A2Con = A2Con,
      A2FLOpen = A2FLOpen,
      A2FLOpt = A2FLOpt,
      A2Fmax= A2Fmax,
      A2FmaxClosed = A2FmaxClosed,
      A2FmaxNoPenn = A2FmaxNoPenn,
      A2Li = A2Li,
      A2Lo = A2Lo,
      A2MA = A2MA,
      A2Mass = A2Mass,
      A2MLClosed = A2MLClosed,
      A2MLOpen = A2MLOpen,
      A2MLOpt = A2MLOpt,
      A2MuscleThetaClosed = A2MuscleThetaClosed,
      A2MuscleThetaOpen = A2MuscleThetaOpen,
      A2Or = A2Or,
      A2OrIns = A2OrIns,
      A2OrInsOpen = A2OrInsOpen,
      A2OrJoint = A2OrJoint,
      A2PCSAclosed = A2PCSAclosed,
      A2PennOpen = A2PennOpen,
      A2Tendon = A2Tendon,
      A2TendonThetaOpen = A2TendonThetaOpen,
      A2ThetaJointClosed = A2ThetaJointClosed,
      A2ThetaJointOpen = A2ThetaJointOpen,
      A3AvgFL = A3AvgFL,
      A3AvgPenn = A3AvgPenn,
      A3AwThetaOpen = A3AwThetaOpen,
      A3Con = A3Con,
      A3FLOpen = A3FLOpen,
      A3FLOpt = A3FLOpt,
      A3Fmax = A3Fmax,
      A3FmaxClosed = A3FmaxClosed,
      A3FmaxNoPenn = A3FmaxNoPenn,
      A3Li = A3Li,
      A3Lo = A3Lo,
      A3MA = A3MA,
      A3Mass = A3Mass,
      A3MLClosed = A3MLClosed,
      A3MLOpen = A3MLOpen,
      A3MLOpt = A3MLOpt,
      A3MuscleThetaClosed = A3MuscleThetaClosed,
      A3MuscleThetaOpen = A3MuscleThetaOpen,
      A3Or = A3Or,
      A3OrIns = A3OrIns,
      A3OrJoint = A3OrJoint,
      A3PCSAclosed = A3PCSAclosed,
      A3PennOpen = A3PennOpen,
      A3Tendon = A3Tendon,
      A3TendonThetaOpen = A3TendonThetaOpen,
      A3ThetaJointClosed = A3ThetaJointClosed,
      A3ThetaJointOpen = A3ThetaJointOpen,
      ActRiseTime = ActRiseTime,
      AwAvgFL = AwAvgFL,
      AwAvgPenn = AwAvgPenn,
      AwCon = AwCon,
      AwFLOpen = AwFLOpen,
      AwFLOpt = AwFLOpt,
      AwFmaxClosed = AwFmaxClosed,
      AwFmaxNoPenn = AwFmaxNoPenn,
      AwFmax= AwFmax,
      AwIns.open = AwIns.open,
      AwIns = AwIns,
      AwLi = AwLi,
      AwLo = AwLo,
      AwMA = AwMA,
      AwMass = AwMass,
      AwMLClosed = AwMLClosed,
      AwMLOpen = AwMLOpen,
      AwMLOpt = AwMLOpt,
      AwMuscleThetaClosed = AwMuscleThetaClosed,
      AwMuscleThetaOpen = AwMuscleThetaOpen,
      AwOr.open = AwOr.open,
      AwOr = AwOr,
      AwOrIns = AwOrIns,
      AwOrInsOpen = AwOrInsOpen,
      AwOrJoint = AwOrJoint,
      AwPCSAclosed = AwPCSAclosed,
      AwPennOpen = AwPennOpen,
      AwTendon = AwTendon,
      AwTendonThetaOpen = AwTendonThetaOpen,
      AwThetaJointClosed = AwThetaJointClosed,
      AwThetaJointOpen = AwThetaJointOpen,
      CdPlate = CdPlate,
      config=config,
      G = G,
      Gape = Gape,
      JawWidth=JawWidth,
      A2inlever.open = inlever.p.open,
      A2inlever = inlever.p,
      Inormal = Inormal,
      Itot = Itot,
      JawTip.open = JawTip.open,
      JawTip = JawTip,
      JL = JL,
      joint.p = joint.p,
      MaxFpar = MaxFpar,
      MaxIso = MaxIso,
      MLOpt = MLOpt,
      MuscleDen = MuscleDen,
      open.ang = open.ang,
      Pressure = Pressure,
      Scale = Scale,
      SL = SL,
      SpecAng = SpecAng,
      start = start,
      Start = Start,
      TendonLengthClosed = TendonLengthClosed,
      TendonThetaOpen=TendonThetaOpen,
      ThetaClosed = ThetaClosed,
      TimeStep = TimeStep,
      TimeStepSecs = TimeStepSecs,
      Vmax = Vmax
    )
  return(start.l)
}

#' @export

close_old <- function(x,MaxIts=1000,tail.n=20){

  config <- x$config
  j.var <- var_df_old(x=x)
  j <- data.frame(matrix(0,ncol=ncol(j.var),nrow=MaxIts))
  j[1,] <- j.var
  colnames(j) <- colnames(j.var)
  j$Iteration <- 1:MaxIts

  pb = txtProgressBar(min = 0, max =1, initial = 2,style=3)

  restart <- TRUE
  geom.l <- list()

  tendon.c <- NULL
  A2ffv.c <- A3ffv.c <- Awffv.c <-NULL
  res.l <- nexus.l <- list()

  for(n in 2:MaxIts){

    per.closed=1-(j[n-1,"JawAng"]-x$ThetaClosed)/(x$Gape-x$ThetaClosed)
    setTxtProgressBar(pb,per.closed)

    if(j[n-1,"JawAng"]<=x$ThetaClosed) break


    j[n,"Time"]<- j[n,"Time"]+x$TimeStep
    j[n,"Iteration"]<- n
    #************* pressure **********
    j[n,"PressureT"]=-0.40373*(pi/2)*(x$a/x$Scale)^2*((x$JawWidth/2)/x$Scale)*((x$Gape-j[n-1,"JawAng"])/x$Gape)*x$Pressure

    if(j[n,"AngVel"]<1e-10)	j[n,"PressureT"]=0

    #************* Acceleration and velocity **********

    j[n,"AngAcc"]=(j[n-1,"AwTorque"]+j[n-1,"TendonTorque"]+j[n-1,"DragT"]+j[n-1,"PressureT"])/x$Itot

    j[n,"AngVel"]=j[n-1,"AngAcc"]*x$TimeStepSecs+j[n-1,"AngVel"];
    j[n,"dJawAng"]=j[n,"AngVel"]*x$TimeStepSecs;
    j[n,"JawAng"]=j[n-1,"JawAng"]-j[n,"dJawAng"]

    open.angle <- -1*(x$start+j[n,"JawAng"]-x$SpecAng)

    #************* Drag torques ************

    j[n,"DragT"]=(-2/15)*1000*j[n-1,"AngVel"]^2*(x$a/x$Scale)^4*(x$JawWidth/2)/x$Scale


    inlever.p <- point_ang_r(x$joint.p,x$A2Li,open.angle)
    JawTip<- point_ang_r(x$joint.p,x$JL,open.angle)

    j[n,"Time"]=j[n-1,"Time"]+x$TimeStep;

    #remove negative muscle forces
    if(j[n-1,"A2Force"]<=0) j[n-1,"A2Force"] <- 1e-15
    if(j[n-1,"A3Force"]<=0) j[n-1,"A3Force"] <- 1e-15
    if(j[n-1,"AwForce"]<=0) j[n-1,"AwForce"] <- 1e-15

    ###### cart coordinates at nth iteration (instantaneous)####

    AwOr<- point_ang_r(x$joint.p,x$AwLi,open.angle)#this is dynamic
    #which line of action inputs most torque

    ####### RESULTANTS  for Aw caclculations #############

    #****** find resultant of A2 and A3 forces
    #need ang between A2 and A3 LOA

    A2A3.ang <- cos_ang(j[n-1,"A3ML"],j[n-1,"A2ML"],dist_2d(x$A2Or,x$A3Or))

    #resultant and its angle
    res <- sqrt(j[n-1,"A2Force"]^2+j[n-1,"A3Force"]^2+2*j[n-1,"A2Force"]*j[n-1,"A3Force"]*cos(A2A3.ang))

    res.theta <- atan((j[n-1,"A2Force"]*sin(A2A3.ang))/(j[n-1,"A3Force"]+j[n-1,"A2Force"]*cos(A2A3.ang)))
    if(res==0) res.theta=A2A3.ang

    #res <- sqrt(j[n,"A2Force"]^2+j[n,"A3Force"]^2+2*j[n,"A2Force"]*j[n,"A3Force"]*cos(A2A3.ang))

    #angle between line def by A2 and A3 origs and A3ML
    A2A3.A3ang <- cos_ang(dist_2d(x$A2Or,x$A3Or),j[n-1,"A3ML"],j[n-1,"A2ML"])

    #theoretical length of resultant forces from A2 and A3
    res.ML <- sin(A2A3.A3ang)*(j[n-1,"A3ML"]/sin(pi-A2A3.A3ang-res.theta))

    #using resultant lengths, compute resultants position
    #distance btween A3 orig and resultant orig
    A3.res <- cos_side(j[n-1,"A3ML"],res.ML,res.theta)
    #cos_side <- function(l,r,ang){sqrt(l^2+ r^2-2*l*r* cos(ang))}
    resOr <- point_ang_r(x$A3Or,
                         A3.res,
                         pi-cos_ang(dist_2d(x$A3Or,c(x$A2Or[1],x$A3Or[2])),dist_2d(x$A3Or,x$A2Or),dist_2d(x$A2Or,c(x$A2Or[1],x$A3Or[2]))) #horiz ang of A3or-A2or
    )

    res.l[[n]] <- resOr

    resOr <- do.call(rbind,res.l) %>% tail(tail.n)%>% colMeans()

    res.l[[n]] <- resOr

    tendon.theta <- cos_ang(x$A2Li,dist_2d(inlever.p,resOr),dist_2d(resOr,x$joint.p)) #tendon ang based on resultants

    tendon.thetaNoAw <- tendon.theta

    nexus <- point_ang_r(inlever.p,j[1,"TendonLength"],pi+open.angle-tendon.theta)


    ######muscle lengths
    m.max <- ifelse(n==2,1,which.max(abs(c(j[n-1,"A2TendonTorque"],j[n-1,"A3TendonTorque"],j[n-1,"AwTendonTorque"]))))


    ### reposition according to sag and equilibrium imposed by Aw on resultant line.
    nexus <- point_ang_r(inlever.p,j[1,"TendonLength"],pi+open.angle-tendon.theta)

    j[n,"AwML"] <- dist_2d(nexus,x$AwOr)
    j[n,"A2ML"] <-dist_2d(nexus,x$A2Or)
    j[n,"A3ML"] <-dist_2d(nexus,x$A3Or)

    ### Add effect of Aw on tendon position, that is, sag imposed by its force

    AwTendonTheta <- cos_ang(j[n,"AwML"],j[1,"TendonLength"],dist_2d(x$AwOr,inlever.p))

    AwTendonF <- sin(AwTendonTheta)*j[n-1,"AwForce"]

    resTendonTheta <- pi-asin(AwTendonF/(AwTendonF+res))


    res.theta.sag <- asin(((resTendonTheta)*j[1,"TendonLength"])/dist_2d(inlever.p,resOr))

    TendonThetaAdj <- pi-res.theta.sag+resTendonTheta

    #### recalculate tendon theta and nexus according to sag (but only for config==3, sling in place)
    if(config==3) tendon.theta <- tendon.theta-TendonThetaAdj #tendon ang based on resultants


    #if(n>2 & tendon.theta>(pi-rad(10))) tendon.theta <- pi-rad(10) # keeps tendon from being flat, must accomodate some thickness off the fossa

    #damp large tendon position changes
    tendon.c <- c(tendon.c,tendon.theta)
    tendon.theta <- mean(tail(tendon.c,tail.n))

    nexus <- point_ang_r(inlever.p,j[1,"TendonLength"],pi+(open.angle-tendon.theta))

    #Meckelian tendon theta ignoring Aw
    noAw.TendonTheta <- ifelse(m.max==1,cos_ang(x$A2Li,dist_2d(inlever.p,x$A2Or),dist_2d(x$A2Or,x$joint.p)),cos_ang(x$A2Li,dist_2d(inlever.p,x$A3Or),dist_2d(x$A3Or,x$joint.p)))
    noAw.nexus <- point_ang_r(inlever.p,j[1,"TendonLength"], pi+open.angle-noAw.TendonTheta)

    if(n==2) {tendon.theta <- noAw.TendonTheta
    nexus <- noAw.nexus
    }

    j[n,"AwML"] <- dist_2d(nexus,AwOr)
    j[n,"A2ML"] <-dist_2d(nexus,x$A2Or)
    j[n,"A3ML"] <-dist_2d(nexus,x$A3Or)

    #}

    j[n,"AwTheta"] <- cos_ang(dist_2d(inlever.p,AwOr),j[n,"AwML"],dist_2d(inlever.p,nexus))


    #for smoothing FFV over prev 2 ms
    if(n>21){range.n <- (n-20):(n-1)}else{range.n <- n-1}

    # nexus is bound by lines of action of A2 and A3
    j[n,"TendonLength"] <-j[1,"TendonLength"]

    j[n,"TendonTheta"]=cos_ang(x$A2Li,j[1,"TendonLength"],dist_2d(x$joint.p,nexus))

    #************* A3 force ********
    j[n,"A3dML"] =j[n-1,"A3ML"]-j[n,"A3ML"];x$TimeStepSecs
    j[n,"A3v"] =j[n,"A3dML"]/x$A3MLOpen/x$TimeStepSecs;
    j[n,"A3FL"] =sqrt((sin(j[n-1,"A3Penn"])*j[n-1,"A3FL"])^2+((cos(j[n-1,"A3Penn"])*j[n-1,"A3FL"])-j[n,"A3dML"])^2);
    j[n,"A3Penn"] =asin((sin(j[n-1,"A3Penn"])*j[n-1,"A3FL"])/j[n,"A3FL"]);
    j[n,"A3Fmax"] =x$A3FmaxNoPenn*cos(j[n,"A3Penn"]);
    A3ffv.c <- c(A3ffv.c,j[n,"A3v"])
    j[n,"A3Ffv"] =vl_ffv(V=mean(tail(A3ffv.c,tail.n)))
    #j[n,"A3Ffv"] =(Vmax-mean(j[range.n,"A3v"]))/(Vmax+mean(j[range.n,"A3v"])*G);
    ifelse(j[n,"A3FL"] >x$A3FLOpt,j[n,"A3Ffl"]<- -6.25*(j[n,"A3FL"]/x$A3FLOpt)^2+12.5*(j[n,"A3FL"]/x$A3FLOpt)-5.25,j[n,"A3Ffl"]  <- 1)

    #ifelse(j[n,"A3FL"]/A3FLOpt >= 0.5 & j[n,"A3FL"]/A3FLOpt <= 1.56, j[n,"A3Ffl"]<- predict.lm(fl.poly,data.frame(p.fl=j[n,"A3FL"]/A3FLOpt)),j[n,"A3Ffl"]  <- 0)# from Porro et al. (2002)

    if(j[n,"Time"]<x$ActRiseTime){j[n,"A3Fact"]=0.5-0.5*cos(pi*j[n,"Time"]/x$ActRiseTime)} else{j[n,"A3Fact"]=1};

    ifelse(j[n,"A3ML"]>x$A3MLOpt,j[n,"A3Fpar"] <- x$A3FmaxClosed*exp(2*log(1+x$MaxFpar)*(j[n,"A3ML"]/x$A3MLOpt-1))-x$A3FmaxClosed,j[n,"A3Fpar"] <- 0)
    j[n,"A3Force"] =j[n,"A3Fmax"] *j[n,"A3Ffv"]*ifelse(j[n,"A3Ffl"]<0,0,j[n,"A3Ffl"])*j[n,"A3Fact"]+j[n,"A3Fpar"];

    j[n,"A3Work"]=j[n,"A3Force"]*	j[n,"A3dML"]/x$Scale
    j[n,"A3Power"]=j[n,"A3Work"]/x$TimeStepSecs
    j[n,"A3relL"]=j[n,"A3ML"]/x$A2MLOpt



    #************* A2 force ********
    j[n,"A2dML"] =j[n-1,"A2ML"]-j[n,"A2ML"];
    j[n,"A2v"] =j[n,"A2dML"]/x$A2MLOpen/x$TimeStepSecs;
    j[n,"A2FL"] =sqrt((sin(j[n-1,"A2Penn"])*j[n-1,"A2FL"])^2+((cos(j[n-1,"A2Penn"])*j[n-1,"A2FL"])-j[n,"A2dML"])^2);
    j[n,"A2Penn"] =asin((sin(j[n-1,"A2Penn"])*j[n-1,"A2FL"])/j[n,"A2FL"]);
    j[n,"A2Fmax"] =x$A2FmaxNoPenn*cos(j[n,"A2Penn"]);
    A2ffv.c <- c(A2ffv.c,j[n,"A2v"])
    j[n,"A2Ffv"] =vl_ffv(V=mean(tail(A2ffv.c,tail.n)))
    #j[n,"A2Ffv"]=(Vmax-mean(j[range.n,"A2v"]))/(Vmax+mean(j[range.n,"A2v"])*G);
    ifelse(j[n,"A2FL"] >x$A2FLOpt,j[n,"A2Ffl"]<- -6.25*(j[n,"A2FL"]/x$A2FLOpt)^2+12.5*(j[n,"A2FL"]/x$A2FLOpt)-5.25,j[n,"A2Ffl"]  <- 1)

    if(j[n,"Time"]<x$ActRiseTime){j[n,"A2Fact"]=0.5-0.5*cos(pi*j[n,"Time"]/x$ActRiseTime)} else{j[n,"A2Fact"]=1};
    ifelse(j[n,"A2ML"]>x$A2MLOpt,j[n,"A2Fpar"] <- x$A2FmaxClosed*exp(2*log(1+x$MaxFpar)*(j[n,"A2ML"]/x$A2MLOpt-1))-x$A2FmaxClosed,j[n,"A2Fpar"] <- 0)

    j[n,"A2Force"] =j[n,"A2Fmax"] *j[n,"A2Ffv"]*ifelse(j[n,"A2Ffl"]<0,0,j[n,"A2Ffl"])*j[n,"A2Fact"]+j[n,"A2Fpar"];

    j[n,"A2Work"]=j[n,"A2Force"]*	j[n,"A2dML"]/x$Scale
    j[n,"A2Power"]=j[n,"A2Work"]/x$TimeStepSecs
    j[n,"A2relL"]=j[n,"A2ML"]/x$A2MLOpt


    #************* Aw force ********
    j[n,"AwdML"] =j[n-1,"AwML"]-j[n,"AwML"];
    j[n,"Awv"] =j[n,"AwdML"]/x$AwMLOpen/x$TimeStepSecs;
    j[n,"AwFL"] =sqrt((sin(j[n-1,"AwPenn"])*j[n-1,"AwFL"])^2+((cos(j[n-1,"AwPenn"])*j[n-1,"AwFL"])-j[n,"AwdML"])^2);
    j[n,"AwPenn"] =asin((sin(j[n-1,"AwPenn"])*j[n-1,"AwFL"])/j[n,"AwFL"]);
    j[n,"AwFmax"] =x$AwFmaxNoPenn*cos(j[n,"AwPenn"]);
    Awffv.c <- c(Awffv.c,j[n,"Awv"])
    j[n,"AwFfv"] =vl_ffv(V=mean(tail(Awffv.c,tail.n)))
    #j[n,"AwFfv"] =(Vmax-mean(j[range.n,"Awv"]))/(Vmax+mean(j[range.n,"Awv"])*G);
    ifelse(j[n,"AwFL"] >x$AwFLOpt,j[n,"AwFfl"]<- -6.25*(j[n,"AwFL"]/x$AwFLOpt)^2+12.5*(j[n,"AwFL"]/x$AwFLOpt)-5.25,j[n,"AwFfl"]  <- 1) # length-tension relationship from Van Wassenbergh


    if(j[n,"Time"]<x$ActRiseTime){j[n,"AwFact"]=0.5-0.5*cos(pi*j[n,"Time"]/x$ActRiseTime)} else{j[n,"AwFact"]=1};
    ifelse(j[n,"AwML"]>x$AwMLOpt,j[n,"AwFpar"] <- x$AwFmaxClosed*exp(2*log(1+x$MaxFpar)*(j[n,"AwML"]/x$AwMLOpt-1))-x$AwFmaxClosed,j[n,"AwFpar"] <- 0)

    j[n,"AwForce"] =j[n,"AwFmax"] *j[n,"AwFfv"]*ifelse(j[n,"AwFfl"]<0,0,j[n,"AwFfl"])*j[n,"AwFact"]+j[n,"AwFpar"];

    j[n,"AwWork"]=j[n,"AwForce"]*	j[n,"AwdML"]/x$Scale
    j[n,"AwPower"]=j[n,"AwWork"]/x$TimeStepSecs
    j[n,"AwrelL"]=j[n,"AwML"]/x$A2MLOpt


    #************** Torques 	#Scaled to N*m  *****************
    #Angles of A2 and A3 with Aw, both constrained to conform to A2-AW line by F condition, see above for d

    j[n,"A3AwTheta"] <-  cos_ang(j[n,"AwML"],j[n,"A3ML"],dist_2d(x$A3Or,AwOr))

    j[n,"A2AwTheta"] <-cos_ang(j[n,"A2ML"],j[n,"AwML"],dist_2d(x$A2Or,AwOr))


    #Angles of A2 and A3 with tendon
    j[n,"A3TendonTheta"] <-cos_ang(j[n,"TendonLength"],j[n,"A3ML"],dist_2d(inlever.p,x$A3Or))

    j[n,"A2TendonTheta"] <-cos_ang(j[n,"TendonLength"],j[n,"A2ML"],dist_2d(inlever.p,x$A2Or))

    j[n,"AwTendonTheta"] <- cos_ang(j[n,"TendonLength"],j[n,"AwML"],dist_2d(AwOr,inlever.p))


    #************** Aw
    #remove negative muscle forces
    if(j[n,"A2Force"]<0) j[n,"A2Force"] <- 1e-15
    if(j[n,"A3Force"]<0) j[n,"A3Force"] <- 1e-15
    if(j[n,"AwForce"]<0) j[n,"AwForce"] <- 1e-15

    j[n,"AwTorque"]=min(c(
      j[n,"AwForce"]*abs(cos(j[n,"AwTheta"])), #Aw component
      j[n,"A2Force"]*abs(cos(j[n,"A2AwTheta"]))+j[n,"A3Force"]*abs(cos(j[n,"A3AwTheta"])) #facialis components
    ))*x$AwLi/x$Scale



    #************** Tendon input

    j[n,"TendonTorque"]=((j[n,"A2Force"]*abs(cos(j[n,"A2TendonTheta"])+j[n,"A3Force"]*abs(cos(j[n,"A3TendonTheta"])))))*sin(j[n,"TendonTheta"])*x$A2Li/x$Scale

    #************** Torques about tendon*****************
    #************** A2
    #Scaled to N*m
    j[n,"A3TendonTorque"]=j[n,"A3Force"]*j[1,"TendonLength"]/x$Scale*sin(j[n,"A3TendonTheta"]);#A3TendonTheta is outside?
    #************** A2
    j[n,"A2TendonTorque"]=j[n,"A2Force"]*j[1,"TendonLength"]/x$Scale*sin(j[n,"A2TendonTheta"]);#A2TendonTheta is outside?
    #************** Aw
    j[n,"AwTendonTorque"]=j[n,"AwForce"]*j[1,"TendonLength"]/x$Scale*sin(j[n,"AwTendonTheta"]);#AwTendonTheta is outside?
    #****** Tendon: Torque applied by A2, A3, Aw at tendon insertion

    #************* performance
    j[n,"AwFout"]= j[n,"AwTorque"]/(x$AwLo/x$Scale) #Aw bite force at tip of the jaw
    j[n,"TendFout"]=j[n,"TendonTorque"]/(x$A2Lo/x$Scale) #A2 bite force at tip of the jaw

    j[n,"JawVelTip"]=j[n,"AngVel"]*(x$A2Lo) #in cm/s
    j[n,"StaticBite"]=j[n,"AwFout"]+j[n,"TendFout"] #Total bite force at tip of jaw
    j[n,"EMA"]=j[n,"StaticBite"]/(j[n,"AwForce"]+j[n,"A2Force"]+j[n,"A3Force"]) #Total bite force at tip of jaw


    plot.ls <- list(joint.p=x$joint.p,A2Or=x$A2Or,A3Or=x$A3Or,JawTip=JawTip,inlever=inlever.p,AwIns=AwOr,nexus=nexus,resOr=resOr)

    plot.df<- as.data.frame(do.call(rbind,plot.ls)) %>% mutate(Time=j[n,"Time"],n=n,JawAng=j[n,"JawAng"])
    colnames(plot.df) <- c("x","y","Time","n","JawAng")
    geom.l[[n]] <- plot.df %>% mutate(name=rownames(plot.df),Time=round(Time,5))

    if(n>3 & restart==TRUE){j[2,] <- j[4,]
    n <- 2
    restart <- FALSE
    }
  }

  j <- j %>% filter(JawAng>=x$ThetaClosed) %>% mutate(Time=round(Time,5))

  #if(min(j %>% filter(JawAng>0) %>% pull(JawAng))>x$ThetaClosed) warning("Max number of interations has been reached and jaw is not closed.")

  geom.df <- as.data.frame(do.call(rbind,geom.l))

  return(list(j=j,geom=geom.df))
}


#' @export

var_df_old<- function(x){
  data.frame(
    Iteration=0,
    Time=x$Start,
    AngAcc=0,
    AngVel=0,
    JawAng=x$Gape,
    dJawAng=0,
    dTendonTheta=0,
    TendonTheta=x$TendonThetaOpen,
    TendonLength=x$TendonLengthClosed,
    TendonAngAcc=0,
    TendonAngVel=0,
    NormalDrag=0,
    LooseDrag=0,
    DragT=0,
    PreyDrag=0,
    PressureT=0,

    ###########  A3
    A3ML=x$A3MLOpen,
    A3dML=0,
    A3v=0,
    A3Penn=0,
    A3FL=x$A3FLOpen, #Instantaneous FL
    A3Fmax=0,
    A3Ffv=1, #Force-velocity factor
    A3Ffl=1, #Force-length factor
    A3Fact=0, #Activation rise-time factor
    A3Fpar=ifelse(x$A3MLOpen>x$A3MLOpt,x$A3Fmax*exp(2*log(1+x$MaxFpar)*(x$A3MLOpen/x$A3MLOpt-1))-x$A3Fmax,0), #Parallel elastic force
    A3Force=0, #Force of the A3 division
    A3Theta=x$A3MuscleThetaOpen, #A3 angle of instertion on lower jaw;
    A3TendonTorque=0, #Torque imparted by A3 on tendon;
    A3Torque=0, #Torque imparted by A3 on lower jaw;
    A3Fout=0,
    A3TendonTheta=x$A3TendonThetaOpen, #angle of A3 relative to tendon
    A3AwTheta=x$A3AwThetaOpen, #angle of A3 relative to Aw

    ################ A2 vectors
    A2ML=x$A2MLOpen,
    A2dML=0,
    A2v=0, #Muscle velocity in ml/s
    A2Penn=x$A2PennOpen, #Instantaneous pennation angle
    A2FL=x$A2FLOpen, #Instantaneous FL
    A2Fmax=0, #Max force produced by A2
    A2Ffv=1, #Force-velocity factor
    A2Ffl=1, #Force-length factor
    A2Fact=0, #Activation rise-time factor
    A2Fpar=ifelse(x$A2MLOpen>x$A2MLOpt,x$A2Fmax*exp(2*log(1+x$MaxFpar)*(x$A2MLOpen/x$A2MLOpt-1))-x$A2Fmax,0), #Parallel elastic force
    A2Force=0, #Force of the A2 division
    A2Theta=x$A2MuscleThetaOpen, #A2 angle of instertion on lower jaw;
    A2TendonTorque=0, #Torque imparted by A2 on tendon;
    A2Torque=0, #Torque imparted by A2 on tendon;
    A2Fout=0,
    A2TendonTheta=x$A2TendonThetaOpen, #angle of A2 relative to tendon
    A2AwTheta=x$A2AwThetaOpen, #angle of A3 relative to Aw

    #Aw vectorsAw
    AwML=x$AwMLOpen,
    AwdML=0,
    Awv=0, #Muscle velocity in ml/s
    AwPenn=x$AwPennOpen, #Instantaneous pennation angle
    AwFL=x$AwFLOpen, #Instantaneous FL
    AwFmax=0, #Max force produced by Aw
    AwFfv=1, #Force-velocity factor
    AwFfl=1, #Force-length factor
    AwFact=0,#Activation rise-time factor
    AwFpar=ifelse(x$AwMLOpen>x$AwMLOpt,x$AwFmax*exp(2*log(1+x$MaxFpar)*(x$AwMLOpen/x$AwMLOpt-1))-x$AwFmax,0), #Parallel elastic force
    AwForce=0,#Force of the Aw division
    AwTheta=x$AwMuscleThetaOpen, #Aw angle of instertion on lower jaw;
    AwTendonTorque=0, #Torque imparted by Aw on tendon;
    AwTorque=0, #Torque imparted by Aw on lower jaw;
    AwFout=0,
    TendFout=0,
    AwTendonTheta=x$AwTendonThetaOpen, #angle of Aw relative to tendon

    #A2/3, Aw vectors
    A2A3Fsum=0,


    ###Tendon
    TendonTorque=0,
    #vectors for static bite force and velocity
    StaticBite=0,
    JawVelTip=0
  )
}



#'    start_old(spec.n=1,config=2,pars=list(Pressure=100,CdPlate=2)) %>%
#'     close_old()
#'     cl.l <- list()
#'
#'

if(F){
  dat<-system.file("extdata", "csloani.csv", package = "JawCloseR") %>% read.csv() %>% rownames_to_column(var="spec.n")

  cl.l <- list()
  for(s in dat$spec.n %>% unique){
    for(c in 1:3 ){
      cl.i <- start_old(spec.n=s,config=c) %>%
        close_old(MaxIts = 1500,tail.n = 50)
      cl.l[[paste0(s,"-",c)]] <- cl.i$j %>% mutate(config=c,spec=s)
    }
  }
  cl <- do.call(rbind,cl.l)
  dat2 <- dat %>% select(spec.n,mass) %>%
    mutate(spec=spec.n)
  cl <- cl %>% left_join(dat2)
  cl %>%
    group_by(spec,config) %>%
    dplyr::summarise(close=max(Time),mass=max(mass)) %>%
    ggplot(aes(log(mass),close,col=as.factor(config)))+geom_point()
  se <- function(x,na.rm=T){if(na.rm==T)
  {x <- x[is.na(x)==F]}else{x <- x}
    sd(x)/sqrt(length(x))}
  cl %>%
    group_by(spec,config) %>%
    mutate(per_closed=1-round_any(JawAng/rad(110),0.02)) %>%
    mutate(per_closed=round(per_closed,2)) %>%
    filter(per_closed %in% seq(0,1,0.02)) %>%
    group_by(config,per_closed) %>%
    dplyr::summarise(A2v.m=mean(A2Power),A2v.se=se(A2Power)) %>%
    ggplot(aes(per_closed,A2v.m,col=as.factor(config)))+geom_point()+geom_errorbar(aes(x=per_closed,ymin=A2v.m-A2v.se,ymax=A2v.m+A2v.se))
}

}
