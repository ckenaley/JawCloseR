
#' @title  Initiate jaw-closing simulations

#' @description  Establishes starting values for important morphological and physiological parameters to be used in adduction simulations with \code{close_jaw()}.
#'
#' @param x a \code{data.frame}, see \code{details}.
#' @param spec.n numeric, which of the specimens in the data file should be analyzed.
#' @param config numeric, the muscle configuration. See \code{details}.
#' @param pars, a list of named numerical values. If \code{NULL}, the default values will be used. See \code{details}.
#'
#' @import tidyverse
#' @export
#'
#'
#' @details
#'This function takes morphometric data of the feeding system and prepares it for analysis with \code{close_jaw}. It essentially uses the input data to open the jaw to a specified or default gape angle and outputs important muscle and model parameters to be used by \code{close_jaw()}.
#'
#'\code{x} should be a \code{data.frame} containing only numeric values with variables/column names specified by \code{jaw_vars()}. The function \code{review_data()} will help new users test data before running \code{close_jaw()}.
#'
#'\code{config} should be 1, 2, or 3, each corresponding to a specific muscle configuration for simulation according to \insertCite{kenaley2019;textual}{JawCloseR}:
#'
#' \itemize{
#' \item \code{config=1}: The \eqn{A_\omega} mass is added to the \eqn{A_2} mass for physiological cross-sectional area calculations. This simplifies simulations, ignoring the geometry of \eqn{A_\omega}, but still accounting for it's force contribution.
#' \item \code{config=2}: The \eqn{A_\omega} mass and geometry are completely ignored.
#' \item \code{config=3}: The \eqn{A_\omega} mass and geometry are accounted for as separate muscle units. This is most biologically relevant, of course, but complicates the model considerably.
#' }
#'
#'
#' @seealso \code{\link{close_jaw()}} \code{\link{jaw_vars()}} \code{\link{review_data()}}
#' @references
#' \insertAllCited{}
#' @examples
#' #Using JawCloseR data from Kenaley et al. (2019) for Chauliodus sloani
#' cs <- start_sim()
#' print(cs)
#'
#' #Using the same, but through entry of a data.frame and the spec.n option
#' cs_dat<-read.csv(system.file("extdata", "csloani.csv", package = "JawCloseR"))
#' print(cs_dat)
#' cs_start <- start_sim(cs_dat,spec.n=2)
#' print(cs_start)
#'
#' #Works with pipe (%>%), too
#' library(tidyverse)
#'
#' cs_cl <- cs_dat %>%
#'    start_sim(spec.n=1,config=3) %>%
#'    close_jaw()
#'    cs_cl$j$Time %>% max()
#'    cs_cl$j %>% ggplot(aes(Time,AngVel %>% deg)) +geom_point()
#'
#' lm_dat<-system.file("extdata", "SMBass_1_muscle_data.csv", package = "JawCloseR") %>% read.csv() %>% rownames_to_column(var="spec.n")
#' lm_cl <- lm_dat %>%
#'    mutate(GapeWidth=5.0) %>%
#'    start_sim(spec.n=1,config=2,pars=list(Pressure=1000,CdPlate=2.4)) %>%
#'    close_jaw()
#'    lm_cl$j$Time %>% max()

start_sim <- function(x=NULL,spec.n=1,config=1,pars=NULL)	{

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

  A23wcols <- c("OrIns","OrJoint","Tendon","Li","Lo","AvgFL","AvgPenn","Mass")

  A1cols <- c("OrIns","OrJoint","Li","Lo","AvgFL","AvgPenn","TendonAntSeg","TendonPostSeg","TendonOrIn","Mass","InsTheta","TendonTheta")


  SpecimenNumber=spec.n;
  #Set morphological variables
  SL=MuscleData[SpecimenNumber,"SL"];
  JL=MuscleData[SpecimenNumber,"JL"];

  #A1 geometry and mass

  if(!any(grepl("A1", colnames(MuscleData)))){

    MuscleData[,paste0("A1",A1cols)] <- 1e-10
  }

  A1OrIns=MuscleData[SpecimenNumber,"A1OrIns"];
  A1OrJoint=MuscleData[SpecimenNumber,"A1OrJoint"]; #change w/ SpecAnge
  A1Li=MuscleData[SpecimenNumber,"A1Li"];
  A1Lo=MuscleData[SpecimenNumber,"A1Lo"];
  A1AvgFL=MuscleData[SpecimenNumber,"A1AvgFL"];#change w/ SpecAnge
  A1AvgPenn=(pi/180)*MuscleData[SpecimenNumber,"A1AvgPenn"]; #change w/ SpecAng
  A1TendonAntSeg <- MuscleData[SpecimenNumber,"A1TendonAntSeg"]  #ant segnment length
  A1TendonPostSeg <- MuscleData[SpecimenNumber,"A1TendonPostSeg"] #post segnment length
  A1TendonOrIn <- MuscleData[SpecimenNumber,"A1TendonOrIn"] #distance betwween two insertions of A1 Tendon
  A1Mass <- MuscleData[SpecimenNumber,"A1Mass"]
  A1InsTheta <- MuscleData[SpecimenNumber,"A1InsTheta"] #insertion angle betwee post seg of A1 Tendon
  A1TendonTheta <- MuscleData[SpecimenNumber,"A1TendonTheta"] #Ins angle of A1 tendon on lower jaw

  if(!any(grepl("A2", colnames(MuscleData)))){

    MuscleData[,paste0("A2",A23wcols)] <- 1e-10
  }
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
  if(!any(grepl("A3", colnames(MuscleData)))){

    MuscleData[,paste0("A3",A23wcols)] <- 1e-10
  }
  A3OrIns =MuscleData[SpecimenNumber,"A3OrIns"]
  A3OrJoint = MuscleData[SpecimenNumber,"A3OrJoint"]
  A3Tendon= MuscleData[SpecimenNumber,"A3Tendon"]
  A3Li =MuscleData[SpecimenNumber,"A3Li"]
  A3Lo <- MuscleData[SpecimenNumber,"A3Lo"]
  A3AvgFL <- MuscleData[SpecimenNumber,"A3AvgFL"]
  A3AvgPenn <- (pi/180)*MuscleData[SpecimenNumber,"A3AvgPenn"]
  A3Mass <- MuscleData[SpecimenNumber,"A3Mass"]

  if(!any(grepl("Aw", colnames(MuscleData)))){

    MuscleData[,paste0("Aw",A23wcols)] <- 1e-10
  }
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


  #fiber length-tension relationship from Porro et al. (2011)
  p.fl <- c(0.508, 0.68, 0.92, 1.0, 1.56) #fl/flo
  p.f <- c(0, 0.84, 1, 1, 0) #%force
  fl.poly <- lm(p.f ~ poly(p.fl, 4, raw=TRUE))


  ThetaClosed <-ParameterValues %>% get_var(v="GapeClosed") %>% rad
  Scale <-  ParameterValues %>% get_var(v="Scale")
  CdPlate <-  ParameterValues %>% get_var(v="CdPlate")
  Start <- TimeStep * -1


  # #prey parameters
  # PreyProp <- prey.per;
  # PreyStrikeDist <- 0.85*A2Lo#MuscleData[SpecimenNumber,"strike"];
  # PreyStrikeAng <-prey.strike.ang*(pi/180);
  #
  # #Set type of Simulation
  # Loosejaw <- Loose;
  # ThinFlat <- prey.pos;

  #Calculate Fmax of each muscle w/ out cos(pennationangle)
  A3FmaxNoPenn <-(A3Mass)/(MuscleDen*A3AvgFL)*2*MaxIso/10;
  A2FmaxNoPenn <-(A2Mass)/(MuscleDen*A2AvgFL)*2*MaxIso/10;
  AwFmaxNoPenn <-(AwMass)/(MuscleDen*AwAvgFL)*2*MaxIso/10;
  A1FmaxNoPenn <-(A1Mass)/(MuscleDen*A1AvgFL)*2*MaxIso/10;

  #Calculate PCSA closed
  A3PCSAclosed <-(A3Mass*cos(A3AvgPenn))/(MuscleDen*A3AvgFL)*2
  A2PCSAclosed <-(A2Mass*cos(A2AvgPenn))/(MuscleDen*A2AvgFL)*2
  AwPCSAclosed <-(AwMass*cos(AwAvgPenn))/(MuscleDen*AwAvgFL)*2
  A1PCSAclosed <-(A1Mass*cos(A1AvgPenn))/(MuscleDen*A1AvgFL)*2

  #Calculate MA
  A1MA <- A1Li/A1Lo
  A3MA <- A3Li/A3Lo
  A2MA <- A2Li/A2Lo
  AwMA <- AwLi/AwLo
  #
  # #Calculate Fmax (N) in closed position
  A1FmaxClosed <-MaxIso/10*A1PCSAclosed
  A3FmaxClosed <-MaxIso/10*A3PCSAclosed
  A2FmaxClosed <-MaxIso/10*A2PCSAclosed
  AwFmaxClosed <-MaxIso/10*AwPCSAclosed


  #*************** Muscle geometry *******************

  #Angle of line between origin-joint and lower jaw when closed
  A3ThetaJointClosed <- acos((A3Li^2+A3OrJoint^2-A3OrIns^2)/(2*A3Li*A3OrJoint))
  A2ThetaJointClosed <- acos((A2Li^2+A2OrJoint^2-A2OrIns^2)/(2*A2Li*A2OrJoint))
  AwThetaJointClosed <- acos((AwLi^2+AwOrJoint^2-AwOrIns^2)/(2*AwLi*AwOrJoint))
  A1ThetaJointClosed <- acos((A1Li^2+A1OrJoint^2-A1OrIns^2)/(2*A1Li*A1OrJoint))

  #Angle of line between origin-joint and lower jaw when open
  A3ThetaJointOpen <- A3ThetaJointClosed+Gape-SpecAng
  A2ThetaJointOpen <- A2ThetaJointClosed+Gape-SpecAng
  AwThetaJointOpen <- AwThetaJointClosed+Gape-SpecAng
  A1ThetaJointOpen <- A1ThetaJointClosed+Gape-SpecAng

  #Angle of insertion on line of actions when jaw is closed
  A3MuscleThetaClosed <-acos((A3Li^2+A3OrIns^2-A3OrJoint^2)/(2*A3Li*A3OrIns))
  A2MuscleThetaClosed <-acos((A2Li^2+A2OrIns^2-A2OrJoint^2)/(2*A2Li*A2OrIns))
  AwMuscleThetaClosed <-acos((AwLi^2+AwOrIns^2-AwOrJoint^2)/(2*AwLi*AwOrIns))
  A1MuscleThetaClosed <-acos((A1Li^2+A1OrIns^2-A1OrJoint^2)/(2*A1Li*A1OrIns))


  ########## cartesian coordinates closed #######
  joint.p <- c(0,0)
  start <- rad(0)
  open.ang <- -1*(start+Gape-SpecAng)
  A2inlever <- point_ang_r(joint.p,A2Li,SpecAng-start) #meckelian tendo ins point shared by A2/A3
  A1inlever <- point_ang_r(joint.p,A1Li,SpecAng-start)

  JawTip <- point_ang_r(joint.p,JL,-start)
  A2Or <- point_ang_r(joint.p,A2OrJoint,A2ThetaJointClosed-start)
  A3Or <- point_ang_r(joint.p,A3OrJoint,A3ThetaJointClosed-start)
  #figure our Aw angles
  AwOr <- point_ang_r(joint.p,AwLi,-start)
  A1Or <- point_ang_r(joint.p,A1OrJoint,A1ThetaJointClosed-start)
  #insertion point on tendon based on muscle theta and orig-ins length (so distal tendon)
  TendonLengthClosed <- sin(AwMuscleThetaClosed)*AwOrIns

  #force Aw to insert on middle of tendon according to muscle theta
  AwIns<- point_ang_r(AwOr,AwLi-A2Li,-rad(180)-AwMuscleThetaClosed-start)


  #####  A1 position closed
  ###  A1 insert on A1 Tendon position closed
  A1Ins <- point_ang_r(A1inlever,A1TendonPostSeg,(rad(180))-A1InsTheta-start) #insertion on A1 tendon (will be dynamics)
  ###  A1 tendon anterior segment position closed
  A1AntSeg <- point_ang_r(A1Ins,A1TendonAntSeg,A1TendonTheta-A1InsTheta)

  ################################################################

  #Length of Muscle (origin to insertion on tendon + proximal tendon) when closed, based on cartesian coords
  A3MLClosed <- dist_2d(AwIns,A3Or)
  A2MLClosed<- dist_2d(AwIns,A2Or)
  AwMLClosed<- dist_2d(AwIns,AwOr)
  A1MLClosed<- dist_2d(A1Ins,A1Or)


  #Length along line of action when open

  A3OrInsOpen <-(A3Li^2+A3OrJoint^2-cos(A3ThetaJointClosed+Gape-SpecAng)*A3OrJoint*A3Li*2)^(0.5)
  A2OrInsOpen <-(A2Li^2+A2OrJoint^2-cos(A2ThetaJointClosed+Gape-SpecAng)*A2OrJoint*A2Li*2)^(0.5)
  AwOrInsOpen <-(AwLi^2+AwOrJoint^2-cos(AwThetaJointClosed+Gape-SpecAng)*AwOrJoint*AwLi*2)^(0.5)
  A1OrInsOpen <-(A1Li^2+A1OrJoint^2-cos(A1ThetaJointClosed+Gape-SpecAng)*A1OrJoint*A1Li*2)^(0.5)



  # geom.df.open <- t(data.frame(joint.p,A3Or,JawTip.open, A2inlever.open,nexus));colnames(geom.df.open) <- c("x","y")
  # geom.df.open.melt <- as.data.frame(geom.df.open) %>% mutate(name=rownames( geom.df.open))
  # p.open <- ggplot(data=geom.df.open.melt,aes(x=x,y=y))+geom_point(aes(colour=name),size=6)+coord_fixed()
  # p.open


  #Angle of along line of action relative to lowerjaw  when jaw is open
  A3MuscleThetaOpen <-acos((A3Li^2+A3OrInsOpen^2-A3OrJoint^2)/(2*A3Li*A3OrInsOpen))
  A2MuscleThetaOpen <-acos((A2Li^2+A2OrInsOpen^2-A2OrJoint^2)/(2*A2Li*A2OrInsOpen))
  A1MuscleThetaOpen <-acos((A1Li^2+A1OrInsOpen^2-A1OrJoint^2)/(2*A1Li*A1OrInsOpen))

  ########## cartesian coordinates open #######

  A2inlever.open <- point_ang_r(joint.p,A2Li,open.ang)
  A1inlever.open <- point_ang_r(joint.p,A1Li,open.ang)
  JawTip.open <- point_ang_r(joint.p,JL,open.ang)


  #### A1 tendon positions, assume straighline between two insertions.

  A1TendonOrIn.open <- dist_2d(A1inlever.open,A1AntSeg)

  A1Ins.open <- point_ang_r(
    A1inlever.open,
    A1TendonPostSeg,
    asin((A1AntSeg[2]-A1inlever.open[2])/ A1TendonOrIn.open)
  )

  #recalc A1 Ant segment
  A1TendonAntSeg <- A1TendonOrIn.open-A1TendonPostSeg

  ### ins angle of post seg of A1 tendon on jaw
  A1InsThetaOpen <- cos_ang(A1Li,A1TendonPostSeg,dist_2d(A1Ins.open,A1inlever.open))

  # geom.df.open <- t(data.frame(joint.p,A2Or,A3Or,JawTip.open,A1inlever.open,A1Or,A1AntSeg,A1Ins.open));colnames(geom.df.open) <- c("x","y")
  #
  #
  # geom.df.open.melt <- as.data.frame(geom.df.open) %>% mutate(name=rownames( geom.df.open))
  # p.open <- ggplot(data=geom.df.open.melt,aes(x=x,y=y))+geom_point(aes(colour=name),size=6)+coord_fixed()
  # p.open

  #figure out Aw angles
  AwOr.open <- point_ang_r(joint.p,AwLi,open.ang)

  #force Aw to insert on middle of tendon according to muscle theta
  #mean of two muscle insertions scaled around 180
  TendonThetaOpen <-mean(c(A2MuscleThetaOpen,A3MuscleThetaOpen))
  AwIns.open<- point_ang_r(A2inlever.open,TendonLengthClosed,(rad(180)-(Gape-SpecAng))-TendonThetaOpen-start) #really the origin and the nexus

  nexus <- AwIns.open

  #where A3,A2, Aw meet at aponeurosis

  A3MLOpen<- dist_2d(AwIns.open,A3Or)
  A2MLOpen<- dist_2d(AwIns.open,A2Or)
  AwMLOpen<- dist_2d(AwIns.open,AwOr)

  #A1 muscle length open
  A1MLOpen<- dist_2d(A1Ins.open,A1Or)

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

  # geom.df.open <- t(data.frame(joint.p,A2Or,A3Or,JawTip.open,A1inlever.open,A1Or,A1AntSeg,A1Ins.open,nexus));colnames(geom.df.open) <- c("x","y")
  # geom.df.open.melt <- as.data.frame(geom.df.open) %>% mutate(name=rownames( geom.df.open))
  # p.open <- ggplot(data=geom.df.open.melt,aes(x=x,y=y))+geom_point(aes(colour=name),size=6)+coord_fixed()
  # p.open

  #Aw insertion angle calculated from open geometry
  AwMuscleThetaOpen <- asin(sin(pi-TendonThetaOpen)/AwMLOpen*TendonLengthClosed)
  #Aw Origin Angle calculated from ^ and subtraction
  AwTendonThetaOpen <-pi-(pi-TendonThetaOpen+AwMuscleThetaOpen)

  #Total muscle contraction
  A3Con <-A3MLOpen-A3MLClosed
  A2Con <-A2MLOpen-A2MLClosed
  AwCon <-AwMLOpen-AwMLClosed
  A1Con <-A1MLOpen-A1MLClosed

  #Fiber length when jaw is open
  A1FLOpen <-sqrt((sin(A1AvgPenn)*A1AvgFL)^2+((cos(A1AvgPenn)*A1AvgFL)+A1Con)^2)
  A3FLOpen <-sqrt((sin(A3AvgPenn)*A3AvgFL)^2+((cos(A3AvgPenn)*A3AvgFL)+A3Con)^2)
  A2FLOpen <-sqrt((sin(A2AvgPenn)*A2AvgFL)^2+((cos(A2AvgPenn)*A2AvgFL)+A2Con)^2)
  AwFLOpen <-sqrt((sin(AwAvgPenn)*AwAvgFL)^2+((cos(AwAvgPenn)*AwAvgFL)+AwCon)^2)

  #Pennation angle open
  A1PennOpen <-asin((sin(A1AvgPenn)*A1AvgFL)/A1FLOpen)
  A3PennOpen <-asin((sin(A3AvgPenn)*A3AvgFL)/A3FLOpen)
  A2PennOpen <-asin((sin(A2AvgPenn)*A2AvgFL)/A2FLOpen)
  AwPennOpen <-asin((sin(AwAvgPenn)*AwAvgFL)/AwFLOpen)

  #Optimal muscle length
  A1MLOpt <-A1MLOpen-(A1MLOpen-A1MLClosed)*MLOpt
  A3MLOpt <-A3MLOpen-(A3MLOpen-A3MLClosed)*MLOpt
  A2MLOpt <- A2MLOpen-(A2MLOpen-A2MLClosed)*MLOpt
  AwMLOpt <- AwMLOpen-(AwMLOpen-AwMLClosed)*MLOpt


  #Optimal fiber length
  A1FLOpt <- sqrt((sin(A1AvgPenn)*A1FLOpen)^2+(cos(A1AvgPenn)*A1FLOpen-(A1MLOpen-A1MLOpt))^2)
  A3FLOpt <- sqrt((sin(A3AvgPenn)*A3FLOpen)^2+(cos(A3AvgPenn)*A3FLOpen-(A3MLOpen-A3MLOpt))^2)
  A2FLOpt <-  sqrt((sin(A2AvgPenn)*A2FLOpen)^2+(cos(A2AvgPenn)*A2FLOpen-(A2MLOpen-A2MLOpt))^2)
  AwFLOpt <-  sqrt((sin(AwAvgPenn)*AwFLOpen)^2+(cos(AwAvgPenn)*AwFLOpen-(AwMLOpen-AwMLOpt))^2)

  #set initial Fmaxs to 0
  A3Fmax=0;A2Fmax=0;AwFmax=0;A1Fmax=0



  #**************** jaw length from the center of the of the axis ****************
  a <- sqrt(A2Lo^2-(JawWidth/2)^2)
  #
  # #**************** Lower-jaw dimensions for loosejaws****************
  # bprime <-SymphLength/a*(JawWidth/2)
  # s <-sqrt(SymphLength^2+bprime^2)
  # l <-a-SymphLength
  # MandRadius <-(MandWidth+MandDepth)/4

  #******************** Mass components for loosejaws *******************
  #
  # #Mass of rami minus symphyseal membrane
  # MandMass <-(MandRadius/Scale)^2*pi*((A3Lo/Scale)-(s/Scale))*1000
  #
  # #Cd of rami ellipsoids after Blevins (1994)
  # CdRami <-0.966*(MandDepth/MandWidth)^(-0.746)

  #Moment of inertia of half ellipse
  Inormal <-(2/15)*pi*1000*(a/Scale)^3*((JawWidth/2)/Scale)^2

  # #Moment of interia of partial ellipse (loosejaw)
  # Ipartellip <-(((JawWidth/2)/Scale)^2*sqrt((a/Scale)^2-(l/Scale)^2)*(2*(a/Scale)^4+(a/Scale)^2*(l/Scale)^2-3*(l/Scale)^4)*pi*1000)/(15*(a/Scale)^2)

  #Moment of inertial of rami (loosejaw)
  # Irami <-(2*(0.25*MandMass*(MandRadius/Scale)^2+(1/3)*MandMass*((A3Lo-s)/Scale)^2))

  #Total moment of inertia (loosejaw)
  # Iloose <-Ipartellip+Irami

  #Moment of inertia of distal tendon, distal to muscle insertions (I=(m*L^2)/3), a thin rigid wire, scale to mass of A3 which is about as thick
  Itendon <- ((TendonLengthClosed/A3MLClosed)*A3Mass/1000*(TendonLengthClosed/Scale)^2)/3

  # #Which jaw mass to use?
  # if(Loosejaw==T){Ijaw <- Iloose}else{if(Loosejaw==F){Ijaw <- Inormal}else{stop("Loose or Normal simulation input unrecognized. Only 'T' or 'F' accepted.")}}
  #
  # if(Prey==T){Itot <- Ijaw+Iprey}else{if(Prey==F){Itot <- Ijaw}else{stop("Prey input unrcognized. Only 'T' or 'F' accepted.")}}

 Ijaw <- Inormal


  geom.df.open <- t(data.frame(joint.p,A2Or,A3Or,JawTip.open,A1inlever.open,A1Ins.open,A1Or,A1AntSeg));colnames(geom.df.open) <- c("x","y")


  geom.df.open.melt <- as.data.frame(geom.df.open) %>% mutate(name=rownames( geom.df.open))
  p.open <- ggplot(data=geom.df.open.melt,aes(x=x,y=y))+geom_point(aes(colour=name),size=6)+coord_fixed()
  p.open

  start.l <-
    list(
      a=a,
      A2AvgFL = A2AvgFL,
      A2AvgPenn = A2AvgPenn,
      A2AwThetaOpen = A2AwThetaOpen,
      A2Con = A2Con,
      A2FLOpen = A2FLOpen,
      A2FLOpt = A2FLOpt,
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
      A1AvgFL = A1AvgFL,
      A1AvgPenn = A1AvgPenn,
      A1Con = A1Con,
      A1FLOpen = A1FLOpen,
      A1FLOpt = A1FLOpt,
      A1FmaxClosed = A1FmaxClosed,
      A1FmaxNoPenn = A1FmaxNoPenn,
      A1Li = A1Li,
      A1Lo = A1Lo,
      A1MA = A1MA,
      A1Mass = A1Mass,
      A1MLClosed = A1MLClosed,
      A1MLOpen = A1MLOpen,
      A1MLOpt = A1MLOpt,
      A1MuscleThetaClosed = A1MuscleThetaClosed,
      A1MuscleThetaOpen = A1MuscleThetaOpen,
      A1Or = A1Or,
      A1Ins=A1Ins.open,
      A1OrIns = A1OrIns,
      A1OrInsOpen = A1OrInsOpen,
      A1OrJoint = A1OrJoint,
      A1PCSAclosed = A1PCSAclosed,
      A1PennOpen = A1PennOpen,
      A1ThetaJointClosed = A1ThetaJointClosed,
      A1ThetaJointOpen = A1ThetaJointOpen,
      A1TendonAntSeg=A1TendonAntSeg,
      A1TendonPostSeg=A1TendonPostSeg,
      A1InsThetaOpen= A1InsThetaOpen,
      A1AntSeg=A1AntSeg,
      CdPlate = CdPlate,
      config=config,
      G = G,
      Gape = Gape,
      Ijaw=Ijaw,
      JawWidth=JawWidth,
      A2inlever.open = A2inlever.open,
      A2inlever = A2inlever,
      A1inlever=A1inlever.open,
      Inormal = Inormal,
      Itendon = Itendon,
      JawTip.open = JawTip.open,
      JawTip = JawTip,
      JawWidth = JawWidth,
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

#' @title  Jaw closing simulations

#' @description  Establishes starting values for important morphological and physiological parameters to be used in adduction simulations with \code{close_jaw()}.
#'
#' @param x a list obtained by \code{start_sim}
#' @param MaxIts numeric, how long should the simulation proceed.
#' @param progress logical, .....
#'
#' @export
#' @import tidyverse
#'
#'
#' @details
#'This function takes morphometric data of the feeding system and prepares it for analysis with \code{close_jaw}. It essentially uses the input data to open the jaw to a specified or default gape angle and outputs important muscle and model parameters to be used by \code{close_jaw}.
#'
#'\code{x} should be a \code{data.frame} or  containing the the variables specified by \code{jaw_vars}.
#'
#' @seealso \code{\link{close_jaw}}
#' @examples
#'
#'cs_dat<-read.csv(system.file("extdata", "csloani.csv", package ="JawCloseR"))
#'cs_dat %>% start_sim()
#' sm <- start_sim(config=1,pars=list(Pressure=100)) %>% close_jaw()
#' #print(sm$ani)
#' sm$j %>% select(Time,DragT,PressureT) %>% pivot_longer(-Time) %>% ggplot(aes(Time,value))+geom_point()+facet_wrap(~name)

close_jaw <- function(x,MaxIts=1000,progress=F){

  #if(any(sapply(x, function(x) length(x)==0))) stop("Start values passed to `x` contain 0-length variables.")

  open <- x
  j.var <- var_df(x=x)
  j <- data.frame(matrix(0,ncol=ncol(j.var),nrow=MaxIts))
  j[1,] <- j.var
  colnames(j) <- colnames(j.var)

  tendon.c <-NULL

  pb = txtProgressBar(min = 2, max = MaxIts, initial = 2,style=3)




  geom.l <- list()

  for(n in 2:MaxIts){
    if(n>2&j[n-1,"JawAng"]<=x$ThetaClosed) break


    setTxtProgressBar(pb,n)
    j[n,"Time"]<- j[n-1,"Time"]+x$TimeStep
    j[n,"Iteration"]<- n
    #************* pressure **********
    j[n,"PressureT"]=-0.40373*(pi/2)*(x$a/x$Scale)^2*((x$JawWidth/2)/x$Scale)*((x$Gape-j[n-1,"JawAng"])/x$Gape)*x$Pressure

    if(j[n-1,"AngVel"]<1e-10)	j[n,"PressureT"]=0

    #************* Acceleration and velocity **********

    j[n,"AngAcc"]=(j[n-1,"A1Torque"]+j[n-1,"AwTorque"]+j[n-1,"TendonTorque"]+j[n-1,"DragT"]+j[n-1,"PressureT"])/x$Ijaw

    j[n,"AngVel"]=j[n-1,"AngAcc"]*x$TimeStepSecs+j[n-1,"AngVel"];
    j[n,"dJawAng"]=j[n,"AngVel"]*x$TimeStepSecs;
    j[n,"JawAng"]=j[n-1,"JawAng"]-j[n,"dJawAng"]

    open.angle <- -1*(x$start+j[n,"JawAng"]-x$SpecAng)

    #************* Drag torques ************

    j[n,"DragT"]=(-2/15)*1000*j[n-1,"AngVel"]^2*(x$a/x$Scale)^4*(x$JawWidth/2)/x$Scale*x$CdPlate*0.5

    #https://www.symbolab.com/solver/integral-calculator/%5Cint_%7B0%7D%5E%7Ba%7D%20%5Cleft(%5Cfrac%7Bb%7D%7Ba%7D%5Cright)%5Csqrt%7B%5Cleft(a%5E%7B2%7D-x%5E%7B2%7D%5Cright)%7D%5Cleft(w%5Ccdot%20x%5Cright)%5E%7B2%7D%5Ccdot%20C%5Ccdot%20p%5Ccdot%20x%5Ccdot%20dx%5Ccdot0.5?or=input


    A2inlever <- point_ang_r(x$joint.p,x$A2Li,open.angle)
    A1inlever <- point_ang_r(x$joint.p,x$A1Li,open.angle)
    JawTip<- point_ang_r(x$joint.p,x$JL,open.angle)

    j[n,"Time"]=j[n-1,"Time"]+x$TimeStep;

    #remove negative muscle forces
    if(j[n-1,"A1Force"]<=0) j[n-1,"A1Force"] <- 1e-15
    if(j[n-1,"A2Force"]<=0) j[n-1,"A2Force"] <- 1e-15
    if(j[n-1,"A3Force"]<=0) j[n-1,"A3Force"] <- 1e-15
    if(j[n-1,"AwForce"]<=0) j[n-1,"AwForce"] <- 1e-15

    ###### cart coordinates at nth iteration (instantaneous)####

    AwOr<- point_ang_r(x$joint.p,x$AwLi,open.angle)#this is dynamic
    #which line of action inputs most torque

    ####### RESULTANTS  #############

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

    tendon.theta <- cos_ang(x$A2Li,dist_2d(A2inlever,resOr),dist_2d(resOr,x$joint.p)) #tendon ang based on resultants

    tendon.thetaNoAw <- tendon.theta

    nexus <- point_ang_r(A2inlever,j[1,"TendonLength"],pi+open.angle-tendon.theta)

    ###### A1

    A1_close <- function(ant_p,inlev_p,ant,post){
      TendonL <- dist_2d(ant_p,inlev_p)
      dTheta <- cos_ang(post,TendonL,ant)
      A1_ins_p <- point_ang_r(inlev_p,post,asin((ant_p[2]-inlev_p[2])/TendonL)+dTheta)

      return(A1_ins_p)
    }

    ### ins of A1 on tendon
    A1Ins <- A1_close(ant_p=x$A1AntSeg,inlev_p=A1inlever,ant=x$A1TendonAntSeg,post=x$A1TendonPostSeg)

    ######muscle lengths
    m.max <- ifelse(n==2,1,which.max(abs(c(j[n-1,"A2TendonTorque"],j[n-1,"A3TendonTorque"],j[n-1,"AwTendonTorque"]))))


    ### reposition according to sag and equilibrium imposed by Aw on resultant line.
    nexus <- point_ang_r(A2inlever,j[1,"TendonLength"],pi+open.angle-tendon.theta)

    j[n,"AwML"] <- dist_2d(nexus,AwOr)
    j[n,"A2ML"] <-dist_2d(nexus,x$A2Or)
    j[n,"A3ML"] <-dist_2d(nexus,x$A3Or)
    j[n,"A1ML"] <-dist_2d(A1Ins,x$A1Or)


    ### Add effect of Aw on tendon position, that is, sag imposed by its force

    AwTendonTheta <- cos_ang(j[n,"AwML"],j[1,"TendonLength"],dist_2d(AwOr,A2inlever))

    AwTendonF <- sin(AwTendonTheta)*j[n-1,"AwForce"]

    resTendonTheta <- pi-asin(AwTendonF/(AwTendonF+res))


    res.theta.sag <- asin(((resTendonTheta)*j[1,"TendonLength"])/dist_2d(A2inlever,resOr))

    TendonThetaAdj <- pi-res.theta.sag+resTendonTheta

    #### recalculate tendon theta and nexus according to sag (but only for config==3, sling in place)
    if(x$config==3) tendon.theta <- tendon.theta-TendonThetaAdj #tendon ang based on resultants


    if(n>2 & tendon.theta>(pi-rad(10))) tendon.theta <- pi-rad(10) # keeps tendon from being flat, must accomodate some thickness off the fossa

    #damp large tendon position changes
    tendon.c <- c(tendon.c,tendon.theta)
    tendon.theta <- mean(tail(tendon.c,20))

    nexus <- point_ang_r(A2inlever,j[1,"TendonLength"],pi+(open.angle-tendon.theta))

    #Meckelian tendon theta ignoring Aw
    noAw.TendonTheta <- ifelse(m.max==1,cos_ang(x$A2Li,dist_2d(A2inlever,x$A2Or),dist_2d(x$A2Or,x$joint.p)),cos_ang(x$A2Li,dist_2d(A2inlever,x$A3Or),dist_2d(x$A3Or,x$joint.p)))
    noAw.nexus <- point_ang_r(A2inlever,j[1,"TendonLength"], pi+open.angle-noAw.TendonTheta)

    if(n==2) {tendon.theta <- noAw.TendonTheta
    nexus <- noAw.nexus
    }

    j[n,"AwML"] <- dist_2d(nexus,AwOr)
    j[n,"A2ML"] <-dist_2d(nexus,x$A2Or)
    j[n,"A3ML"] <-dist_2d(nexus,x$A3Or)

    #need post of A1 ins on A1 tendon




    #}



    if(x$config!=3){j[n,"AwTheta"] <- 0}else
      {j[n,"AwTheta"] <- cos_ang(dist_2d(A2inlever,AwOr),j[n,"AwML"],j[1,"TendonLength"])}#ignore Aw more or less if not config

    #for smoothing FFV over prev 2 ms
    if(n>21){range.n <- (n-20):(n-1)}else{range.n <- n-1}

    # nexus is bound by lines of action of A2 and A3
    j[n,"TendonLength"] <-j[1,"TendonLength"]

    j[n,"TendonTheta"]=cos_ang(x$A2Li,j[1,"TendonLength"],dist_2d(x$joint.p,nexus))

    #************* A3 force ********
    j[n,"A3dML"] =j[n-1,"A3ML"]-j[n,"A3ML"]
    j[n,"A3v"] =j[n,"A3dML"]/x$A3MLOpen/x$TimeStepSecs;
    j[n,"A3FL"] =sqrt((sin(j[n-1,"A3Penn"])*j[n-1,"A3FL"])^2+((cos(j[n-1,"A3Penn"])*j[n-1,"A3FL"])-j[n,"A3dML"])^2);
    j[n,"A3Penn"] =asin((sin(j[n-1,"A3Penn"])*j[n-1,"A3FL"])/j[n,"A3FL"]);
    j[n,"A3Fmax"] =x$A3FmaxNoPenn*cos(j[n,"A3Penn"]);
    j[n,"A3Ffv"] =vl_ffv(V=mean(j[range.n,"A3v"]),Vmax = x$Vmax)
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
    j[n,"A2Ffv"] =vl_ffv(V=mean(j[range.n,"A2v"]),Vmax = x$Vmax);
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
    j[n,"AwFfv"] =vl_ffv(V=mean(j[range.n,"Awv"]),Vmax = x$Vmax)
    #j[n,"AwFfv"] =(Vmax-mean(j[range.n,"Awv"]))/(Vmax+mean(j[range.n,"Awv"])*G);
    ifelse(j[n,"AwFL"] >x$AwFLOpt,j[n,"AwFfl"]<- -6.25*(j[n,"AwFL"]/x$AwFLOpt)^2+12.5*(j[n,"AwFL"]/x$AwFLOpt)-5.25,j[n,"AwFfl"]  <- 1) # length-tension relationship from Van Wassenbergh

    if(j[n,"Time"]<x$ActRiseTime){j[n,"AwFact"]=0.5-0.5*cos(pi*j[n,"Time"]/x$ActRiseTime)} else{j[n,"AwFact"]=1};
    ifelse(j[n,"AwML"]>x$AwMLOpt,j[n,"AwFpar"] <- x$AwFmaxClosed*exp(2*log(1+x$MaxFpar)*(j[n,"AwML"]/x$AwMLOpt-1))-x$AwFmaxClosed,j[n,"AwFpar"] <- 0)

    j[n,"AwForce"] =j[n,"AwFmax"] *j[n,"AwFfv"]*ifelse(j[n,"AwFfl"]<0,0,j[n,"AwFfl"])*j[n,"AwFact"]+j[n,"AwFpar"];

    j[n,"AwWork"]=j[n,"AwForce"]*	j[n,"AwdML"]/x$Scale
    j[n,"AwPower"]=j[n,"AwWork"]/x$TimeStepSecs
    j[n,"AwrelL"]=j[n,"AwML"]/x$A2MLOpt


    #************* A1 force ********
    j[n,"A1dML"] =j[n-1,"A1ML"]-j[n,"A1ML"];
    j[n,"A1v"] =j[n,"A1dML"]/x$A1MLOpen/x$TimeStepSecs;
    j[n,"A1FL"] =sqrt((sin(j[n-1,"A1Penn"])*j[n-1,"A1FL"])^2+((cos(j[n-1,"A1Penn"])*j[n-1,"A1FL"])-j[n,"A1dML"])^2);
    j[n,"A1Penn"] =asin((sin(j[n-1,"A1Penn"])*j[n-1,"A1FL"])/j[n,"A1FL"]);
    j[n,"A1Fmax"] =x$A1FmaxNoPenn*cos(j[n,"A1Penn"]);
    j[n,"A1Ffv"] =vl_ffv(V=mean(j[range.n,"A1v"]),Vmax = x$Vmax);
    #j[n,"A1Ffv"]=(Vmax-mean(j[range.n,"A1v"]))/(Vmax+mean(j[range.n,"A1v"])*G);
    ifelse(j[n,"A1FL"] >x$A1FLOpt,j[n,"A1Ffl"]<- -6.25*(j[n,"A1FL"]/x$A1FLOpt)^2+12.5*(j[n,"A1FL"]/x$A1FLOpt)-5.25,j[n,"A1Ffl"]  <- 1)

    if(j[n,"Time"]<x$ActRiseTime){j[n,"A1Fact"]=0.5-0.5*cos(pi*j[n,"Time"]/x$ActRiseTime)} else{j[n,"A1Fact"]=1};
    ifelse(j[n,"A1ML"]>x$A1MLOpt,j[n,"A1Fpar"] <- x$A1FmaxClosed*exp(2*log(1+x$MaxFpar)*(j[n,"A1ML"]/x$A1MLOpt-1))-x$A1FmaxClosed,j[n,"A1Fpar"] <- 0)

    j[n,"A1Force"] =j[n,"A1Fmax"] *j[n,"A1Ffv"]*ifelse(j[n,"A1Ffl"]<0,0,j[n,"A1Ffl"])*j[n,"A1Fact"]+j[n,"A1Fpar"];

    j[n,"A1Work"]=j[n,"A1Force"]*	j[n,"A1dML"]/x$Scale
    j[n,"A1Power"]=j[n,"A1Work"]/x$TimeStepSecs
    j[n,"A1relL"]=j[n,"A1ML"]/x$A1MLOpt



    #************** Torques 	#Scaled to N*m  *****************
    #Angles of A2 and A3 with Aw, both constrained to conform to A2-AW line by F condition, see above for d

    j[n,"A3AwTheta"] <-  cos_ang(j[n,"AwML"],j[n,"A3ML"],dist_2d(x$A3Or,AwOr))

    j[n,"A2AwTheta"] <-cos_ang(j[n,"A2ML"],j[n,"AwML"],dist_2d(x$A2Or,AwOr))


    #Angles of A2 and A3 with tendon
    j[n,"A3TendonTheta"] <-cos_ang(j[n,"TendonLength"],j[n,"A3ML"],dist_2d(A2inlever,x$A3Or))

    j[n,"A2TendonTheta"] <-cos_ang(j[n,"TendonLength"],j[n,"A2ML"],dist_2d(A2inlever,x$A2Or))

    j[n,"AwTendonTheta"] <- cos_ang(j[n,"TendonLength"],j[n,"AwML"],dist_2d(AwOr,A2inlever))

    #here

    #************** Aw
    #remove negative muscle forces
    if(j[n,"A2Force"]<0) j[n,"A2Force"] <- 1e-15
    if(j[n,"A3Force"]<0) j[n,"A3Force"] <- 1e-15
    if(j[n,"AwForce"]<0) j[n,"AwForce"] <- 1e-15
    if(j[n,"A1Force"]<0) j[n,"A1Force"] <- 1e-15


    j[n,"AwTorque"]=min(c(
      j[n,"AwForce"]*abs(cos(j[n,"AwTheta"])), #Aw component
      j[n,"A2Force"]*abs(cos(j[n,"A2AwTheta"]))+j[n,"A3Force"]*abs(cos(j[n,"A3AwTheta"])) #facialis components
    ))*x$AwLi/x$Scale



    #************** Tendon inputs
    #*A2/3 and w input

    j[n,"TendonTorque"]=((j[n,"A2Force"]*abs(cos(j[n,"A2TendonTheta"])+j[n,"A3Force"]*abs(cos(j[n,"A3TendonTheta"])))))*sin(j[n,"TendonTheta"])*x$A2Li/x$Scale



    #************** Torques about tendon*****************
    #************** A2
    #Scaled to N*m
    j[n,"A3TendonTorque"]=j[n,"A3Force"]*j[1,"TendonLength"]/x$Scale*sin(j[n,"A3TendonTheta"]);
    #************** A2
    j[n,"A2TendonTorque"]=j[n,"A2Force"]*j[1,"TendonLength"]/x$Scale*sin(j[n,"A2TendonTheta"]);
    #************** Aw
    j[n,"AwTendonTorque"]=j[n,"AwForce"]*j[1,"TendonLength"]/x$Scale*sin(j[n,"AwTendonTheta"]);
    #****** Tendon: Torque applied by A2, A3, Aw at tendon insertion
    #*
    #************** A1 Tendon input *****************
    #* insertion angle of post seg of a
    j[n,"A1TendonTheta"]=cos_ang(x$A1Li,x$A1TendonPostSeg,dist_2d(x$joint.p,A1Ins))

    A1_tension <- function(orig_p,ins_p,inlev_p,ant_p,f){
      alpha <- abs(cos_ang(dist_2d(inlev_p,ins_p),dist_2d(orig_p,ins_p),dist_2d(inlev_p,orig_p)))
      if(alpha>pi/2) alpha <- alpha-pi/2

      beta <- cos_ang(dist_2d(orig_p,ins_p),dist_2d(ant_p,ins_p),dist_2d(orig_p,ant_p)) %>% abs
      if(beta>pi/2) beta <- beta-pi/2

      tension <- f/(cos(alpha)*sin(beta)/cos(beta)+sin(alpha))

      return(tension)
    }


    j[n,"A1TendonForce"] <- A1_tension(x$A1Or,A1Ins,A1inlever, x$A1AntSeg,f=j[n,"A1Force"])
    j[n,"A1Torque"]=j[n,"A1TendonForce"]*sin(j[n,"A1TendonTheta"]);

    #************* performance
    j[n,"AwFout"]= j[n,"AwTorque"]/(x$AwLo/x$Scale) #Aw bite force at tip of the jaw
    j[n,"TendFout"]=j[n,"TendonTorque"]/(x$A2Lo/x$Scale) #A2 bite force at tip of the jaw

    j[n,"JawVelTip"]=j[n,"AngVel"]*(x$A2Lo) #in cm/s
    j[n,"StaticBite"]=j[n,"AwFout"]+j[n,"TendFout"] #Total bite force at tip of jaw
    j[n,"EMA"]=j[n,"StaticBite"]/(j[n,"AwForce"]+j[n,"A2Force"]+j[n,"A3Force"]) #Total bite force at tip of jaw

    #for plotting
    plot.ls <- list(joint.p=x$joint.p,A2Or=x$A2Or,A3Or=x$A3Or,JawTip=JawTip,inlever=A2inlever,AwIns=AwOr,nexus=nexus,resOr=resOr,A1Ins=A1Ins,AntSeg=x$A1AntSeg)

   plot.df<- as.data.frame(do.call(rbind,plot.ls)) %>% mutate(Time=j[n,"Time"],n=n,JawAng=j[n,"JawAng"])
   colnames(plot.df) <- c("x","y","Time","n","JawAng")
    geom.l[[n]] <- plot.df %>% mutate(name=rownames(plot.df),Time=round(Time,5))

    j[n,"A1Force"] <- 0

  }

  j <- j %>% filter(JawAng>=x$ThetaClosed) %>% mutate(Time=round(Time,5))



  geom.df <- as.data.frame(do.call(rbind,geom.l))

  #geom.out[[paste(n)]] <- data.frame(ms=j[n,"Time"],ang=j[n,"JawAng"],var=unlist(plot.ls))



  p <- geom.df %>%
    ggplot(aes(x,y,col=name))+geom_point()+coord_fixed()+transition_time(as.integer(n))+labs(title = "frame: {m} ")+xlim(range(geom.df$x))
p


  return(list(open=open,j=j,pos=geom.df,ani=p))
}


#' @title  View variables needed for simulation
#' @description  View variables needed for simulation and their description
#'L
#' @details
#'This function produces a \code{data.frame} of the names and descriptions of the variables needed for \code{start_sim}. Users of \code{JawCloseR} should reference this before inputing data to \code{start_sim}.
#'
#' @export
#'
#' @seealso \code{\link{start_sime}}
#' @examples
#'
#'
#'
jaw_vars <- function(){
  data.frame(name=c(
  "JL",
  "A2OrIns",
  "A2OrJoint",
  "A2Tendon",
  "A2Li",
  "A2Lo",
  "A2AvgFL",
  "A2AvgPenn",
  "A2Mass",
  "A3OrIns",
  "A3OrJoint",
  "A3Tendon",
  "A3Li",
  "A3Lo",
  "A3AvgFL",
  "A3AvgPenn",
  "A3Mass",
  "AwOrIns",
  "AwOrJoint",
  "AwTendon",
  "AwLi",
  "AwLo",
  "AwAvgFL",
  "AwAvgPenn",
  "AwMass",
  "SpecAng",
  "Gape",
  "GapeWidth"))


}


# `_old` functions . . .
#code from github repo
# produces results similar to MS Feb 23.

#' @export

start_old <- function(x=NULL,spec.n=3,config=1,pars=NULL)	{
  #*********** File Control ***************
  if(!is.null(x)){ w <- tryCatch(review_data(x),warning = function(w) w)

  if("warning" %in% class(w)) warning("Input data passed to `x` are broken. They may have missing variables, non-numeric data, etc. Please use `reveiew_data()`.")
  }
  if(!is.null(x) & !is.data.frame(x)) stop("'x' and data passed to `start_sim()` must be a data.frame")

  #Load specimen data and parameter values
  if(!is.null(x)) MuscleData<-x


  if(is.null(x)) MuscleData<-system.file("extdata", "csloani.csv", package = "JawCloseR") %>% read.csv() %>% rownames_to_column(var="spec.n")

  ParameterValues <-system.file("extdata", "par_val.csv", package = "JawCloseR") %>% read.csv()



  if(!is.null(pars)){  if (any(!names(pars) %in% ParameterValues$Parameter))
    stop(
      '`pars` must be either NULL or a named list with names in c("TimeStepMS","Density","MaxIso","MaxFpar","Vmax","Gcon","k",ActivTime","MLOptimum","GapeClosed","Scale","CdPlate","MLmin","Pressure")'
    )
  }

  if (!is.null(pars))
    for (p in names(pars))
      ParameterValues[ParameterValues$Parameter == p, ]$Value <- pars[[which(names(pars) == p)]]

  if(nrow(MuscleData)==1) spec.n <- 1
  SpecimenNumber=spec.n;
  #Set morphological variables
  SL=MuscleData[SpecimenNumber,"SL"];
  JL=MuscleData[SpecimenNumber,"JL"];

  A23wcols <- c("OrIns","OrJoint","Tendon","Li","LiTheta","Lo","AvgFL","AvgPenn","Mass")

  A1cols <- c("OrIns","OrJoint","Li","LiTheta","Lo","AvgFL","AvgPenn","TendonAntSeg","TendonPostSeg","TendonOrIn","Mass","InsTheta","TendonTheta")


  #A1 geometry and mass

  if(!any(grepl("A1", colnames(MuscleData)))){

    MuscleData[,paste0("A1",A1cols)] <- 1e-10
  }

  A1OrIns=MuscleData[SpecimenNumber,"A1OrIns"];
  A1OrJoint=MuscleData[SpecimenNumber,"A1OrJoint"]; #change w/ SpecAnge
  A1Li=MuscleData[SpecimenNumber,"A1Li"];
  A1LiTheta=MuscleData[SpecimenNumber,"A1LiTheta"]; #ang between orgin, li, and horiz
  A1Lo=MuscleData[SpecimenNumber,"A1Lo"];
  A1AvgFL=MuscleData[SpecimenNumber,"A1AvgFL"];#change w/ SpecAnge
  A1AvgPenn=(pi/180)*MuscleData[SpecimenNumber,"A1AvgPenn"]; #change w/ SpecAng
  A1TendonAntSeg <- MuscleData[SpecimenNumber,"A1TendonAntSeg"]  #ant segnment length
  A1TendonPostSeg <- MuscleData[SpecimenNumber,"A1TendonPostSeg"] #post segnment length
  A1TendonOrIn <- MuscleData[SpecimenNumber,"A1TendonOrIn"] #distance betwween two insertions of A1 Tendon
  A1Mass <- MuscleData[SpecimenNumber,"A1Mass"]
  A1InsTheta <- MuscleData[SpecimenNumber,"A1InsTheta"] #insertion angle between post seg of A1 Tendon
  A1TendonTheta <- MuscleData[SpecimenNumber,"A1TendonTheta"] #Ins angle of A1 tendon on lower jaw

  #A2 geometry and mass

  if(!any(grepl("A2", colnames(MuscleData)))){

    MuscleData[,paste0("A2",A23wcols)] <- 1e-10
  }


  A2OrIns=MuscleData[SpecimenNumber,"A2OrIns"];
  A2OrJoint=MuscleData[SpecimenNumber,"A2OrJoint"]; #change w/ SpecAnge
  A2Tendon=MuscleData[SpecimenNumber,"A2Tendon"]
  A2Li=MuscleData[SpecimenNumber,"A2Li"];
  A2LiTheta=MuscleData[SpecimenNumber,"A2LiTheta"];
  A2Lo=MuscleData[SpecimenNumber,"A2Lo"];
  A2AvgFL=MuscleData[SpecimenNumber,"A2AvgFL"];#change w/ SpecAnge
  A2AvgPenn=(pi/180)*MuscleData[SpecimenNumber,"A2AvgPenn"]; #change w/ SpecAng
  A2Mass=ifelse(config==1,MuscleData[SpecimenNumber,"A2Mass"]+MuscleData[SpecimenNumber,"AwMass"],MuscleData[SpecimenNumber,"A2Mass"]);


  #A3 geometry and mass

  if(!any(grepl("A3", colnames(MuscleData)))){

    MuscleData[,paste0("A3",A23wcols)] <- 1e-10
  }

  A3OrIns =MuscleData[SpecimenNumber,"A3OrIns"]
  A3OrJoint = MuscleData[SpecimenNumber,"A3OrJoint"]
  A3Tendon= MuscleData[SpecimenNumber,"A3Tendon"]
  A3Li =MuscleData[SpecimenNumber,"A3Li"]
  A3LiTheta=MuscleData[SpecimenNumber,"A3LiTheta"];
  A3Lo <- MuscleData[SpecimenNumber,"A3Lo"]
  A3AvgFL <- MuscleData[SpecimenNumber,"A3AvgFL"]
  A3AvgPenn <- (pi/180)*MuscleData[SpecimenNumber,"A3AvgPenn"]
  A3Mass <- MuscleData[SpecimenNumber,"A3Mass"]

  #Aw geometry and mass
  if(!any(grepl("Aw", colnames(MuscleData)))){

    MuscleData[,paste0("Aw",A23wcols)] <- 1e-10


   if(config==1) warning("Aw components are incomplete and ignored")

  }

  AwOrIns <- MuscleData[SpecimenNumber,"AwOrIns"]
  AwOrJoint <- MuscleData[SpecimenNumber,"AwOrJoint"]
  AwTendon <- MuscleData[SpecimenNumber,"AwTendon"]
  AwLi <- MuscleData[SpecimenNumber,"AwLi"]
  AwLiTheta=MuscleData[SpecimenNumber,"AwLiTheta"];
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
  k <- ParameterValues %>% get_var(v="k")
  #Activation tise time
  ActRiseTime <- ParameterValues %>% get_var(v="ActivTime")
  MLOpt <- ParameterValues %>% get_var(v="MLOptimum")/100

  Pressure <-ParameterValues %>% get_var(v="Pressure")


  ThetaClosed <-ParameterValues %>% get_var(v="GapeClosed") %>% rad
  Scale <-  ParameterValues %>% get_var(v="Scale")
  CdPlate <-  ParameterValues %>% get_var(v="CdPlate")
  Start <- TimeStep * -1


  #Calculate Fmax of each muscle w/ out cos(pennationangle)
  A1FmaxNoPenn <-(A1Mass)/(MuscleDen*A1AvgFL)*2*MaxIso/10;
  A3FmaxNoPenn <-(A3Mass)/(MuscleDen*A3AvgFL)*2*MaxIso/10;
  A2FmaxNoPenn <-(A2Mass)/(MuscleDen*A2AvgFL)*2*MaxIso/10;
  AwFmaxNoPenn <-(AwMass)/(MuscleDen*AwAvgFL)*2*MaxIso/10;

  #Calculate PCSA closed
  A1PCSAclosed <-(A1Mass*cos(A1AvgPenn))/(MuscleDen*A1AvgFL)*2
  A3PCSAclosed <-(A3Mass*cos(A3AvgPenn))/(MuscleDen*A3AvgFL)*2
  A2PCSAclosed <-(A2Mass*cos(A2AvgPenn))/(MuscleDen*A2AvgFL)*2
  AwPCSAclosed <-(AwMass*cos(AwAvgPenn))/(MuscleDen*AwAvgFL)*2

  #Calculate MA
  A1MA <- A1Li/A3Lo
  A3MA <- A3Li/A3Lo
  A2MA <- A2Li/A2Lo
  AwMA <- AwLi/AwLo
  #
  # #Calculate Fmax (N) in closed position
  A1FmaxClosed <-MaxIso/10*A1PCSAclosed
  A3FmaxClosed <-MaxIso/10*A3PCSAclosed
  A2FmaxClosed <-MaxIso/10*A2PCSAclosed
  AwFmaxClosed <-MaxIso/10*AwPCSAclosed


  #*************** Muscle geometry *******************

  #Angle of line between origin-joint and inlever when closed
  A1ThetaJointClosed <- acos((A1Li^2+A1OrJoint^2-A1OrIns^2)/(2*A1Li*A1OrJoint))
  A3ThetaJointClosed <- acos((A3Li^2+A3OrJoint^2-A3OrIns^2)/(2*A3Li*A3OrJoint))
  A2ThetaJointClosed <- acos((A2Li^2+A2OrJoint^2-A2OrIns^2)/(2*A2Li*A2OrJoint))
  AwThetaJointClosed <- acos((AwLi^2+AwOrJoint^2-AwOrIns^2)/(2*AwLi*AwOrJoint))

  #Angle of line between origin-joint and inlever w when open
  A1ThetaJointOpen <- A1ThetaJointClosed+Gape-SpecAng
  A3ThetaJointOpen <- A3ThetaJointClosed+Gape-SpecAng
  A2ThetaJointOpen <- A2ThetaJointClosed+Gape-SpecAng
  AwThetaJointOpen <- AwThetaJointClosed+Gape-SpecAng

  #Angle of insertion on line of actions when jaw is closed
  A1MuscleThetaClosed <-acos((A1Li^2+A1OrIns^2-A1OrJoint^2)/(2*A1Li*A1OrIns))
  A3MuscleThetaClosed <-acos((A3Li^2+A3OrIns^2-A3OrJoint^2)/(2*A3Li*A3OrIns))
  A2MuscleThetaClosed <-acos((A2Li^2+A2OrIns^2-A2OrJoint^2)/(2*A2Li*A2OrIns))
  AwMuscleThetaClosed <-acos((AwLi^2+AwOrIns^2-AwOrJoint^2)/(2*AwLi*AwOrIns))



  ########## cartesian coordinates closed #######
  joint.p <- c(0,0)
  start <- rad(-0)


  JawTip <- point_ang_r(joint.p,JL,-start)

  A1Or <- point_ang_r(joint.p,A1OrJoint,A1ThetaJointClosed-start+A1LiTheta)
  A2Or <- point_ang_r(joint.p,A2OrJoint,A2ThetaJointClosed-start+A2LiTheta)
  A3Or <- point_ang_r(joint.p,A3OrJoint,A3ThetaJointClosed-start+A3LiTheta)


  #force Aw to insert on middle of tendon according to muscle theta
  if(config %in% c(1:3)){AwIns<- point_ang_r(AwOr,AwLi-A2Li,-rad(180)-AwMuscleThetaClosed-start)}else{
    AwIns <- point_ang_r(joint.p,AwOrJoint,AwThetaJointClosed-start+AwLiTheta)
  }


  ## Inlever positions openw
  A3inlever <- point_ang_r(joint.p,A3Li,A3LiTheta)
  A2inlever <- point_ang_r(joint.p,A2Li,A2LiTheta) #meckelian tendo ins point shared by A2/A3
  A1inlever <- point_ang_r(joint.p,A1Li,A1LiTheta)
  Awinlever <- point_ang_r(joint.p,AwLi,AwLiTheta)

  #figure our Aw angles
  AwOr <- point_ang_r(joint.p,AwLi,AwLiTheta)
  #insertion point on tendon based on muscle theta and orig-ins length (so distal tendon)
  TendonLengthClosed <- sin(AwMuscleThetaClosed)*AwOrIns


  #####  A1 position closed
  ###  A1 insert on A1 Tendon position closed
  A1Ins <- point_ang_r(A1inlever,A1TendonPostSeg,(rad(180))-A1InsTheta-start+A1LiTheta) #insertion on A1 tendon (will be dynamics)
  ###  A1 tendon anterior segment position closed
  A1AntSeg <- point_ang_r(A1Ins,A1TendonAntSeg,A1TendonTheta-A1InsTheta+A1LiTheta)

#   geom.df <- t(data.frame(joint.p=joint.p,A1AntSeg=A1AntSeg,A1Ins=A1Ins,A1Or=A1Or,A3Or=A3Or,JawTip=JawTip,A1inlever=A1inlever,A3inlever=A3inlever,AwIns=AwOr,AwOr=AwIns));colnames(geom.df) <- c("x","y")
#   geom.df.melt <- as.data.frame(geom.df) %>% mutate(name=rownames(geom.df))
#   p <- ggplot(data=geom.df.melt,aes(x=x,y=y))+geom_point(aes(colour=name),size=6)+coord_fixed()
# print(p)




  ################################################################

  #Length of Muscle (origin to insertion on tendon + proximal tendon) when closed, based on cartesian coords
  A1MLClosed<- dist_2d(A1Ins,A1Or)


  if(config %in% 1:3){
    A2MLClosed<- dist_2d(AwIns,A2Or)
  A3MLClosed <- dist_2d(AwIns,A3Or)
  AwMLClosed<- dist_2d(AwIns,AwOr)
  }else{
    A2MLClosed<- dist_2d(A2inlever,A2Or)
    A3MLClosed <- dist_2d(A3inlever,A3Or)
    AwMLClosed<- dist_2d(Awinlever,AwIns)
  }


  #Length along line of action when open
  A1OrInsOpen <-(A1Li^2+A1OrJoint^2-cos(A1ThetaJointClosed+Gape-SpecAng)*A1OrJoint*A1Li*2)^(0.5)
  A3OrInsOpen <-(A3Li^2+A3OrJoint^2-cos(A3ThetaJointClosed+Gape-SpecAng)*A3OrJoint*A3Li*2)^(0.5)
  A2OrInsOpen <-(A2Li^2+A2OrJoint^2-cos(A2ThetaJointClosed+Gape-SpecAng)*A2OrJoint*A2Li*2)^(0.5)
  AwOrInsOpen <-(AwLi^2+AwOrJoint^2-cos(AwThetaJointClosed+Gape-SpecAng)*AwOrJoint*AwLi*2)^(0.5)

  #Angle of along line of action relative to lowerjaw  when jaw is open
  A1MuscleThetaOpen <-acos((A1Li^2+A1OrInsOpen^2-A1OrJoint^2)/(2*A1Li*A1OrInsOpen))
  A3MuscleThetaOpen <-acos((A3Li^2+A3OrInsOpen^2-A3OrJoint^2)/(2*A3Li*A3OrInsOpen))
  A2MuscleThetaOpen <-acos((A2Li^2+A2OrInsOpen^2-A2OrJoint^2)/(2*A2Li*A2OrInsOpen))

  ########## cartesian coordinates open #######
  open.ang <- -1*(start+Gape-SpecAng)
  A2inlever.open <- rotate_coords(A2inlever,open.ang,center = joint.p) #point_ang_r(joint.p,A2Li,open.ang)
  A1inlever.open <- rotate_coords(A1inlever,open.ang,center = joint.p)
  A3inlever.open <- rotate_coords(A3inlever,open.ang,center = joint.p)
  Awinlever.open <- rotate_coords(Awinlever,open.ang,center = joint.p)
  JawTip.open <- rotate_coords(JawTip,open.ang,center = joint.p)

  #### A1 tendon positions, assume straighline between two insertions.

  A1TendonOrIn.open <- dist_2d(A1inlever.open,A1AntSeg)

  A1Ins.open <- point_ang_r(
    A1inlever.open,
    A1TendonPostSeg,
    asin((A1AntSeg[2]-A1inlever.open[2])/ A1TendonOrIn.open)
)




  #recalc A1 Ant segment
  A1TendonAntSeg <- A1TendonOrIn.open-A1TendonPostSeg

  ### ins angle of post seg of A1 tendon on jaw
  A1InsThetaOpen <- cos_ang(A1Li,A1TendonPostSeg,dist_2d(A1Ins.open,A1inlever.open))



  #figure out Aw angles
  AwOr.open <- Awinlever.open


  # geom.df.open <- t(data.frame(joint.p,A2Or,A3Or,AwOr.open,AwIns,A1Or,JawTip.open,A1inlever.open,A1AntSeg,A1Ins.open));colnames(geom.df.open) <- c("x","y")
  #
  #
  # geom.df.open.melt <- as.data.frame(geom.df.open) %>% mutate(name=rownames( geom.df.open))
  # p.open <- ggplot(data=geom.df.open.melt,aes(x=x,y=y))+geom_point(aes(colour=name),size=6)+coord_fixed()
  # p.open

  #force Aw to insert on middle of tendon according to muscle theta
  #mean of two muscle insertions scaled around 180
  TendonThetaOpen <-mean(c(A2MuscleThetaOpen,A3MuscleThetaOpen))
  AwIns.open<- point_ang_r(A2inlever.open,TendonLengthClosed,(rad(180)-(Gape-SpecAng))-TendonThetaOpen-start) #really the origin and the nexus

  nexus <- AwIns.open


  A1MLOpen<- dist_2d(A1Ins.open,A1Or)


  if(config %in% 1:3){#where A3,A2, Aw meet at aponeurosis, A1 on A1 tendon
  A3MLOpen<- dist_2d(AwIns.open,A3Or)
  A2MLOpen<- dist_2d(AwIns.open,A2Or)
  AwMLOpen<- dist_2d(AwIns.open,AwOr)
  }else{
    A3MLOpen<- dist_2d(A3inlever.open,A3Or)
    A2MLOpen<- dist_2d(A3inlever.open,A2Or)
    AwMLOpen<- dist_2d(Awinlever.open,AwIns.open)
  }


  A3TendonThetaOpen <- TendonThetaOpen-A3MuscleThetaOpen
  A2TendonThetaOpen <- A2MuscleThetaOpen-TendonThetaOpen

  #angle between MLs and Aw MLs
  if(config %in% c(1:3)){
    A3AwThetaOpen <-cos_ang(AwMLOpen,A3MLOpen,dist_2d(AwOr.open,A3Or))
    A2AwThetaOpen <- cos_ang(A2MLOpen,AwMLOpen,dist_2d(AwOr.open,A2Or))
  }else{
      A3AwThetaOpen <- NA
      A2AwThetaOpen <- NA}




  #get signs right
  geom.df.open <- t(data.frame(joint.p,A3Or,A1Or,JawTip.open,A1inlever.open,A3inlever.open,AwOr.open,
                               AwIns.open,A1Ins.open,A1AntSeg));colnames(geom.df.open) <- c("x","y")
  geom.df.open.melt <- as.data.frame(geom.df.open) %>% mutate(name=rownames(geom.df.open))
  p.open <- ggplot(data=geom.df.open.melt,aes(x=x,y=y))+geom_point(aes(colour=name),size=6)+coord_fixed()
#print(p.open)

  #png(filename=paste(graph.dir,"/geom_open.j.png",sep=""),width=700,height=700)
  #print(p)
  #dev.off()
  ##############


  #Aw insertion angle calculated from open geometry

  AwMuscleThetaOpen <- asin(sin(pi-TendonThetaOpen)/AwMLOpen*TendonLengthClosed)
  #Aw Origin Angle calculated from ^ and subtraction
  AwTendonThetaOpen <-pi-(pi-TendonThetaOpen+AwMuscleThetaOpen)

  #Total muscle contraction
  A1Con <-A1MLOpen-A1MLClosed
  A3Con <-A3MLOpen-A3MLClosed
  A2Con <-A2MLOpen-A2MLClosed
  AwCon <-AwMLOpen-AwMLClosed

  #Fiber length when jaw is open
  A1FLOpen <-sqrt((sin(A1AvgPenn)*A1AvgFL)^2+((cos(A1AvgPenn)*A1AvgFL)+A1Con)^2)
  A3FLOpen <-sqrt((sin(A3AvgPenn)*A3AvgFL)^2+((cos(A3AvgPenn)*A3AvgFL)+A3Con)^2)
  A2FLOpen <-sqrt((sin(A2AvgPenn)*A2AvgFL)^2+((cos(A2AvgPenn)*A2AvgFL)+A2Con)^2)
  AwFLOpen <-sqrt((sin(AwAvgPenn)*AwAvgFL)^2+((cos(AwAvgPenn)*AwAvgFL)+AwCon)^2)

  #Pennation angle open
  A1PennOpen <-asin((sin(A1AvgPenn)*A1AvgFL)/A1FLOpen)
  A3PennOpen <-asin((sin(A3AvgPenn)*A3AvgFL)/A3FLOpen)
  A2PennOpen <-asin((sin(A2AvgPenn)*A2AvgFL)/A2FLOpen)
  AwPennOpen <-asin((sin(AwAvgPenn)*AwAvgFL)/AwFLOpen)

  #Optimal muscle length
  A1MLOpt <-A1MLOpen-(A1MLOpen-A1MLClosed)*MLOpt
  A3MLOpt <-A3MLOpen-(A3MLOpen-A3MLClosed)*MLOpt
  A2MLOpt <- A2MLOpen-(A2MLOpen-A2MLClosed)*MLOpt
  AwMLOpt <- AwMLOpen-(AwMLOpen-AwMLClosed)*MLOpt


  #Optimal fascicle length
  A1FLOpt <- sqrt((sin(A1AvgPenn)*A1FLOpen)^2+(cos(A1AvgPenn)*A1FLOpen-(A1MLOpen-A1MLOpt))^2)
  A3FLOpt <- sqrt((sin(A3AvgPenn)*A3FLOpen)^2+(cos(A3AvgPenn)*A3FLOpen-(A3MLOpen-A3MLOpt))^2)
  A2FLOpt <-  sqrt((sin(A2AvgPenn)*A2FLOpen)^2+(cos(A2AvgPenn)*A2FLOpen-(A2MLOpen-A2MLOpt))^2)
  AwFLOpt <-  sqrt((sin(AwAvgPenn)*AwFLOpen)^2+(cos(AwAvgPenn)*AwFLOpen-(AwMLOpen-AwMLOpt))^2)

  #set initial Fmaxs to 0
  A3Fmax=0;A2Fmax=0;AwFmax=0;A1Fmax=0
  #**************** Lower-jaw dimensions for loosejaws****************

  a <- sqrt(A3Lo^2-(JawWidth/2)^2)
  bprime <-SymphLength/a*(JawWidth/2)
  s <-sqrt(SymphLength^2+bprime^2)
  l <-a-SymphLength
  MandRadius <-(MandWidth+MandDepth)/4

  #******************** Mass components *******************


  #Moment of inertia of half ellipse
  Inormal <-(2/15)*pi*1000*(a/Scale)^3*((JawWidth/2)/Scale)^2

  Itot <- Inormal
  #plot start

  start.l <-
    list(
      a=a,
      A1AvgFL = A1AvgFL,
      A1AvgPenn = A1AvgPenn,
      A1Con = A1Con,
      A1FLOpen = A1FLOpen,
      A1FLOpt = A1FLOpt,
      A1FmaxClosed = A1FmaxClosed,
      A1FmaxNoPenn = A1FmaxNoPenn,
      A1Li = A1Li,
      A1Lo = A1Lo,
      A1MA = A1MA,
      A1Mass = A1Mass,
      A1MLClosed = A1MLClosed,
      A1MLOpen = A1MLOpen,
      A1MLOpt = A1MLOpt,
      A1MuscleThetaClosed = A1MuscleThetaClosed,
      A1MuscleThetaOpen = A1MuscleThetaOpen,
      A1Or = A1Or,
      A1Ins=A1Ins.open,
      A1OrIns = A1OrIns,
      A1OrInsOpen = A1OrInsOpen,
      A1OrJoint = A1OrJoint,
      A1PCSAclosed = A1PCSAclosed,
      A1PennOpen = A1PennOpen,
      A1ThetaJointClosed = A1ThetaJointClosed,
      A1ThetaJointOpen = A1ThetaJointOpen,
      A1TendonAntSeg=A1TendonAntSeg,
      A1TendonPostSeg=A1TendonPostSeg,
      A1InsThetaOpen= A1InsThetaOpen,
      A1AntSeg=A1AntSeg,
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
      k=k,
      A2inlever.open = A2inlever.open,
      A2inlever = A2inlever,
      A1inlever.open=A1inlever.open,
      A1inlever=A1inlever,
      A3inlever.open=A3inlever.open,
      A3inlever=A3inlever,
      Awinlever.open=Awinlever.open,
      Awinlever=Awinlever,
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
  A2ffv.c <- A3ffv.c <- Awffv.c <- A1ffv.c<-0
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

    j[n,"DragT"]=(-2/15)*1000*j[n-1,"AngVel"]^2*(x$a/x$Scale)^4*(x$JawWidth/2)/x$Scale*x$CdPlate*0.5

    A3inlever <- point_ang_r(x$joint.p,x$A3Li,open.angle)
    A2inlever <- point_ang_r(x$joint.p,x$A2Li,open.angle)
    A1inlever <- point_ang_r(x$joint.p,x$A1Li,open.angle)
    JawTip<- point_ang_r(x$joint.p,x$JL,open.angle)

    j[n,"Time"]=j[n-1,"Time"]+x$TimeStep;

    #remove negative muscle forces
    if(j[n-1,"A1Force"]<=0) j[n-1,"A1Force"] <- 1e-15
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

    #Angle of Meckelian tendon
    Mtendon.theta <- cos_ang(x$A2Li,dist_2d(A2inlever,resOr),dist_2d(resOr,x$joint.p)) #tendon ang based on resultants

    tendon.thetaNoAw <- tendon.theta

    nexus <- point_ang_r(A2inlever,j[1,"TendonLength"],pi+open.angle-tendon.theta)


    ######muscle lengths
    m.max <- ifelse(n==2,1,which.max(abs(c(j[n-1,"A2TendonTorque"],j[n-1,"A3TendonTorque"],j[n-1,"AwTendonTorque"]))))


    ### reposition according to sag and equilibrium imposed by Aw on resultant line.
    nexus <- point_ang_r(A2inlever,j[1,"TendonLength"],pi+open.angle-tendon.theta)


    ###### A1 insertion on A1 tendon

    A1_close <- function(ant_p,inlev_p,ant,post){
      TendonL <- dist_2d(ant_p,inlev_p)
      dTheta <- cos_ang(post,TendonL,ant)
      A1_ins_p <- point_ang_r(inlev_p,post,asin((ant_p[2]-inlev_p[2])/TendonL)+dTheta)

      return(A1_ins_p)
    }

    ### ins of A1 on tendon
    if(x$config==4) A1Ins <- A1_close(ant_p=x$A1AntSeg,inlev_p=A1inlever,ant=x$A1TendonAntSeg,post=x$A1TendonPostSeg) else A1Ins <- x$joint.p

    ######muscle lengths
    j[n,"A1ML"] <-dist_2d(A1Ins,x$A1Or)
    j[n,"AwML"] <- dist_2d(nexus,x$AwOr)
    j[n,"A2ML"] <-dist_2d(nexus,x$A2Or)
    j[n,"A3ML"] <-dist_2d(nexus,x$A3Or)

    ### Add effect of Aw on tendon position, that is, sag imposed by its force

    AwTendonTheta <- cos_ang(j[n,"AwML"],j[1,"TendonLength"],dist_2d(x$AwOr,A2inlever))

    AwTendonF <- sin(AwTendonTheta)*j[n-1,"AwForce"]

    resMtendonTheta <- pi-asin(AwTendonF/(AwTendonF+res))


    res.theta.sag <- asin(((resTendonTheta)*j[1,"TendonLength"])/dist_2d(A2inlever,resOr))

    MtendonThetaAdj <- pi-res.theta.sag+resMtendonTheta

    #### recalculate tendon theta and nexus according to sag (but only for config==3, sling in place)
    if(config==3) Mtendon.theta <- Mtendon.theta-MtendonThetaAdj #tendon ang based on resultants


    #if(n>2 & tendon.theta>(pi-rad(10))) Mtendon.theta <- pi-rad(10) # keeps tendon from being flat, must accomodate some thickness off the fossa

    #damp large tendon position changes
    tendon.c <- c(tendon.c,Mtendon.theta)
    Mtendon.theta <- mean(tail(tendon.c,tail.n))
    nexus <- point_ang_r(A2inlever,j[1,"TendonLength"],pi+(open.angle-Mtendon.theta))

    #Meckelian tendon theta ignoring Aw
    noAw.MtendonTheta <- ifelse(m.max==1,cos_ang(x$A2Li,dist_2d(A2inlever,x$A2Or),dist_2d(x$A2Or,x$joint.p)),cos_ang(x$A2Li,dist_2d(A2inlever,x$A3Or),dist_2d(x$A3Or,x$joint.p)))

    noAw.nexus <- point_ang_r(A2inlever,j[1,"TendonLength"], pi+open.angle-noAw.TendonTheta)

    if(n==2) {Mtendon.theta <- noAw.MtendonTheta
    nexus <- noAw.nexus
    }


    j[n,"AwML"] <- dist_2d(nexus,AwOr)
    j[n,"A2ML"] <-dist_2d(nexus,x$A2Or)
    j[n,"A3ML"] <-dist_2d(nexus,x$A3Or)

    #}

    j[n,"AwTheta"] <- cos_ang(dist_2d(A2inlever,AwOr),j[n,"AwML"],dist_2d(A2inlever,nexus))

    # nexus is bound by lines of action of A2 and A3

    j[n,"TendonLength"] <-j[1,"TendonLength"]
    j[n,"TendonTheta"]=cos_ang(x$A2Li,j[1,"TendonLength"],dist_2d(x$joint.p,nexus))

    #************* A1 force ********
    j[n,"A1dML"] =j[n-1,"A1ML"]-j[n,"A1ML"];
    j[n,"A1v"] =j[n,"A1dML"]/x$A1MLOpen/x$TimeStepSecs;
    j[n,"A1FL"] =sqrt((sin(j[n-1,"A1Penn"])*j[n-1,"A1FL"])^2+((cos(j[n-1,"A1Penn"])*j[n-1,"A1FL"])-j[n,"A1dML"])^2);
    j[n,"A1Penn"] =asin((sin(j[n-1,"A1Penn"])*j[n-1,"A1FL"])/j[n,"A1FL"]);
    j[n,"A1Fmax"] =x$A1FmaxNoPenn*cos(j[n,"A1Penn"]);
    A1ffv.c <- c(A1ffv.c,j[n,"A1v"])
    j[n,"A1Ffv"] =vl_ffv(V=mean(tail(A1ffv.c,tail.n)),Vmax=x$Vmax,k=x$k)
    #j[n,"A1Ffv"]=(Vmax-mean(j[range.n,"A1v"]))/(Vmax+mean(j[range.n,"A1v"])*G);
    ifelse(j[n,"A1FL"] >x$A1FLOpt,j[n,"A1Ffl"]<- -6.25*(j[n,"A1FL"]/x$A1FLOpt)^2+12.5*(j[n,"A1FL"]/x$A1FLOpt)-5.25,j[n,"A1Ffl"]  <- 1)

    if(j[n,"Time"]<x$ActRiseTime){j[n,"A1Fact"]=0.5-0.5*cos(pi*j[n,"Time"]/x$ActRiseTime)} else{j[n,"A1Fact"]=1};
    ifelse(j[n,"A1ML"]>x$A1MLOpt,j[n,"A1Fpar"] <- x$A1FmaxClosed*exp(2*log(1+x$MaxFpar)*(j[n,"A1ML"]/x$A1MLOpt-1))-x$A1FmaxClosed,j[n,"A1Fpar"] <- 0)

    j[n,"A1Force"] =j[n,"A1Fmax"] *j[n,"A1Ffv"]*ifelse(j[n,"A1Ffl"]<0,0,j[n,"A1Ffl"])*j[n,"A1Fact"]+j[n,"A1Fpar"];

    j[n,"A1Work"]=j[n,"A1Force"]*	j[n,"A1dML"]/x$Scale
    j[n,"A1Power"]=j[n,"A1Work"]/x$TimeStepSecs
    j[n,"A1relL"]=j[n,"A1ML"]/x$A1MLOpt


    #************* A3 force ********
    j[n,"A3dML"] =j[n-1,"A3ML"]-j[n,"A3ML"];x$TimeStepSecs
    j[n,"A3v"] =j[n,"A3dML"]/x$A3MLOpen/x$TimeStepSecs;
    j[n,"A3FL"] =sqrt((sin(j[n-1,"A3Penn"])*j[n-1,"A3FL"])^2+((cos(j[n-1,"A3Penn"])*j[n-1,"A3FL"])-j[n,"A3dML"])^2);
    j[n,"A3Penn"] =asin((sin(j[n-1,"A3Penn"])*j[n-1,"A3FL"])/j[n,"A3FL"]);
    j[n,"A3Fmax"] =x$A3FmaxNoPenn*cos(j[n,"A3Penn"]);
    A3ffv.c <- c(A3ffv.c,j[n,"A3v"])
    j[n,"A3Ffv"] =vl_ffv(V=mean(tail(A3ffv.c,tail.n)),Vmax=x$Vmax,k=x$k)
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
    j[n,"A2Ffv"] =vl_ffv(V=mean(tail(A2ffv.c,tail.n)),Vmax=x$Vmax,k=x$k)
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
    j[n,"AwFfv"] =vl_ffv(V=mean(tail(Awffv.c,tail.n)),Vmax=x$Vmax,k=x$k)
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
    j[n,"A3TendonTheta"] <-cos_ang(j[n,"TendonLength"],j[n,"A3ML"],dist_2d(A2inlever,x$A3Or))

    j[n,"A2TendonTheta"] <-cos_ang(j[n,"TendonLength"],j[n,"A2ML"],dist_2d(A2inlever,x$A2Or))

    j[n,"AwTendonTheta"] <- cos_ang(j[n,"TendonLength"],j[n,"AwML"],dist_2d(AwOr,A2inlever))


    #************** Aw
    #remove negative muscle forces
    if(j[n,"A1Force"]<0) j[n,"A1Force"] <- 1e-15
    if(j[n,"A2Force"]<0) j[n,"A2Force"] <- 1e-15
    if(j[n,"A3Force"]<0) j[n,"A3Force"] <- 1e-15
    if(j[n,"AwForce"]<0) j[n,"AwForce"] <- 1e-15

    j[n,"AwTorque"]=min(c(
      j[n,"AwForce"]*abs(cos(j[n,"AwTheta"])), #Aw component
      j[n,"A2Force"]*abs(cos(j[n,"A2AwTheta"]))+j[n,"A3Force"]*abs(cos(j[n,"A3AwTheta"])) #facialis components
    ))*x$AwLi/x$Scale



    #************** A2/3/2 Tendon input

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


    #************** A1 Tendon input *****************
    #* insertion angle of post seg of a
    #*
    #*
    j[n,"A1TendonTheta"]=cos_ang(x$A1Li,x$A1TendonPostSeg,dist_2d(x$joint.p,A1Ins))

    # calculate tension in post seg of A1 tendon
    A1_tension <- function(orig_p,ins_p,inlev_p,ant_p,f){
      alpha <- abs(cos_ang(dist_2d(inlev_p,ins_p),dist_2d(orig_p,ins_p),dist_2d(inlev_p,orig_p)))
      if(alpha>pi/2) alpha <- alpha-pi/2
      beta <- cos_ang(dist_2d(orig_p,ins_p),dist_2d(ant_p,ins_p),dist_2d(orig_p,ant_p)) %>% abs
      if(beta>pi/2) beta <- beta-pi/2
      tension <- f/(cos(alpha)*sin(beta)/cos(beta)+sin(alpha))
      return(tension)
    }


    j[n,"A1TendonForce"] <- A1_tension(x$A1Or,A1Ins,A1inlever, x$A1AntSeg,f=j[n,"A1Force"])
    j[n,"A1Torque"]=j[n,"A1TendonForce"]*sin(j[n,"A1TendonTheta"]);
    j[n,"A1Fout"]= j[n,"A1Torque"]/(x$A1Lo/x$Scale)

    #************* performance
    j[n,"AwFout"]= j[n,"AwTorque"]/(x$AwLo/x$Scale) #Aw bite force at tip of the jaw
    j[n,"TendFout"]=j[n,"TendonTorque"]/(x$A2Lo/x$Scale) #A2 bite force at tip of the jaw

    j[n,"JawVelTip"]=j[n,"AngVel"]*(x$A2Lo) #in cm/s
    j[n,"StaticBite"]=j[n,"AwFout"]+j[n,"TendFout"]+j[n,"A1Fout"] #Total bite force at tip of jaw
    j[n,"EMA"]=j[n,"StaticBite"]/(j[n,"AwForce"]+j[n,"A2Force"]+j[n,"A3Force"]+j[n,"A1Force"]) #Total bite force at tip of jaw


    plot.ls <- list(joint.p=x$joint.p,A2Or=x$A2Or,A3Or=x$A3Or,JawTip=JawTip,inlever=A2inlever,AwIns=AwOr,nexus=nexus,resOr=resOr,A1Ins=A1Ins,AntSeg=x$A1AntSeg)

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

close_old_config4 <- function(x,MaxIts=1000,tail.n=20){

  open <- x
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
  A2ffv.c <- A3ffv.c <- Awffv.c <- A1ffv.c<-0
  res.l <- nexus.l <- list()

  for(n in 2:MaxIts){

    per.closed=1-(j[n-1,"JawAng"]-x$ThetaClosed)/(x$Gape-x$ThetaClosed)
    setTxtProgressBar(pb,per.closed)

    if(j[n-1,"JawAng"]<=x$ThetaClosed) break


    j[n,"Time"]<- j[n,"Time"]+x$TimeStep
    j[n,"Iteration"]<- n
    #************* pressure **********
    j[n,"PressureT"]=-0.40373*(pi/2)*(x$a/x$Scale)^2*((x$JawWidth/2)/x$Scale)*((x$Gape-j[n-1,"JawAng"])/x$Gape)*x$Pressure

    if(j[n-1,"AngVel"]<1e-10)	j[n,"PressureT"]=0

    #************* Acceleration and velocity **********

    j[n,"AngAcc"]=(j[n-1,"AwTorque"]+j[n-1,"A3Torque"]+j[n-1,"A1Torque"]+j[n-1,"DragT"]+j[n-1,"PressureT"])/x$Itot

    j[n,"AngVel"]=j[n-1,"AngAcc"]*x$TimeStepSecs+j[n-1,"AngVel"];
    j[n,"dJawAng"]=j[n,"AngVel"]*x$TimeStepSecs;
    j[n,"JawAng"]=j[n-1,"JawAng"]-j[n,"dJawAng"]

    open.angle <- -1*(x$start+j[n,"JawAng"]-x$SpecAng)

    #************* Drag torques ************

    j[n,"DragT"]=(-2/15)*1000*j[n-1,"AngVel"]^2*(x$a/x$Scale)^4*(x$JawWidth/2)/x$Scale*x$CdPlate*0.5

    #Inlever positions

    A3inlever <- rotate_coords(x$A3inlever,open.angle,x$joint.p)
    A1inlever <- rotate_coords(x$A1inlever,open.angle,x$joint.p)
    Awinlever <- rotate_coords(x$Awinlever,open.angle,x$joint.p)

    JawTip<- rotate_coords(x$JawTip,open.angle,x$joint.p)

    j[n,"Time"]=j[n-1,"Time"]+x$TimeStep;

    #remove negative muscle forces
    if(j[n-1,"A1Force"]<=0) j[n-1,"A1Force"] <- 1e-15
    if(j[n-1,"A3Force"]<=0) j[n-1,"A3Force"] <- 1e-15
    if(j[n-1,"AwForce"]<=0) j[n-1,"AwForce"] <- 1e-15

    ###### cart coordinates at nth iteration (instantaneous)####

    ###### A1 insertion on A1 tendon

    A1_close <- function(ant_p,inlev_p,ant,post){
      TendonL <- dist_2d(ant_p,inlev_p)
      dTheta <- cos_ang(post,TendonL,ant)
      A1_ins_p <- point_ang_r(inlev_p,post,asin((ant_p[2]-inlev_p[2])/TendonL)+dTheta)

      return(A1_ins_p)
    }


    ### ins of A1 on tendon
   A1Ins <- A1_close(ant_p=x$A1AntSeg,inlev_p=A1inlever,ant=x$A1TendonAntSeg,post=x$A1TendonPostSeg)

    ######muscle lengths
    j[n,"A1ML"] <-dist_2d(A1inlever,x$A1Or)
    j[n,"AwML"] <- dist_2d(Awinlever,x$AwIns)
    j[n,"A3ML"] <-dist_2d(A3inlever,x$A3Or)

    ### insertion angles A3 and Aw

    A3Theta <- cos_ang(x$A3Li,j[n,"A3ML"],x$A3OrJoint)
    AwTheta <- cos_ang(x$AwLi,j[n,"AwML"],x$AwOrJoint)
    j[n,"AwTheta"] <- AwTheta
    j[n,"A3Theta"] <- A3Theta
    #for plotting
    plot.ls <- list(joint.p=x$joint.p,A1Or=x$A1Or,A3Or=x$A3Or,JawTip=JawTip,A3inlever=A3inlever,AwIns=Awinlever,AwOr=x$AwIns,A1Ins=A1Ins,AntSeg=x$A1AntSeg,A1inlever)

    # plot.df<- as.data.frame(do.call(rbind,plot.ls)) %>% mutate(Time=j[n,"Time"],n=n,JawAng=j[n,"JawAng"])
    # colnames(plot.df) <- c("x","y","Time","n","JawAng")


    # plot.df %>%
    #   mutate(name=rownames(plot.df),Time=round(Time,5)) %>%
    # ggplot(aes(x,y,col=name))+geom_point()+coord_fixed()



    #************* A1 force ********
    j[n,"A1dML"] =j[n-1,"A1ML"]-j[n,"A1ML"];
    j[n,"A1v"] =j[n,"A1dML"]/x$A1MLOpen/x$TimeStepSecs;
    j[n,"A1FL"] =sqrt((sin(j[n-1,"A1Penn"])*j[n-1,"A1FL"])^2+((cos(j[n-1,"A1Penn"])*j[n-1,"A1FL"])-j[n,"A1dML"])^2);
    j[n,"A1Penn"] =asin((sin(j[n-1,"A1Penn"])*j[n-1,"A1FL"])/j[n,"A1FL"]);
    j[n,"A1Fmax"] =x$A1FmaxNoPenn*cos(j[n,"A1Penn"]);
    A1ffv.c <- c(A1ffv.c,j[n,"A1v"])
    j[n,"A1Ffv"] =vl_ffv(V=mean(tail(A1ffv.c,tail.n)),Vmax=x$Vmax,k=x$k)
    #j[n,"A1Ffv"]=(Vmax-mean(j[range.n,"A1v"]))/(Vmax+mean(j[range.n,"A1v"])*G);
    ifelse(j[n,"A1FL"] >x$A1FLOpt,j[n,"A1Ffl"]<- -6.25*(j[n,"A1FL"]/x$A1FLOpt)^2+12.5*(j[n,"A1FL"]/x$A1FLOpt)-5.25,j[n,"A1Ffl"]  <- 1)

    if(j[n,"Time"]<x$ActRiseTime){j[n,"A1Fact"]=0.5-0.5*cos(pi*j[n,"Time"]/x$ActRiseTime)} else{j[n,"A1Fact"]=1};
    ifelse(j[n,"A1ML"]>x$A1MLOpt,j[n,"A1Fpar"] <- x$A1FmaxClosed*exp(2*log(1+x$MaxFpar)*(j[n,"A1ML"]/x$A1MLOpt-1))-x$A1FmaxClosed,j[n,"A1Fpar"] <- 0)

    j[n,"A1Force"] =j[n,"A1Fmax"] *j[n,"A1Ffv"]*ifelse(j[n,"A1Ffl"]<0,0,j[n,"A1Ffl"])*j[n,"A1Fact"]+j[n,"A1Fpar"];

    j[n,"A1Work"]=j[n,"A1Force"]*	j[n,"A1dML"]/x$Scale
    j[n,"A1Power"]=j[n,"A1Work"]/x$TimeStepSecs
    j[n,"A1relL"]=j[n,"A1ML"]/x$A1MLOpt


    #************* A3 force ********
    j[n,"A3dML"] =j[n-1,"A3ML"]-j[n,"A3ML"];x$TimeStepSecs
    j[n,"A3v"] =j[n,"A3dML"]/x$A3MLOpen/x$TimeStepSecs;
    j[n,"A3FL"] =sqrt((sin(j[n-1,"A3Penn"])*j[n-1,"A3FL"])^2+((cos(j[n-1,"A3Penn"])*j[n-1,"A3FL"])-j[n,"A3dML"])^2);
    j[n,"A3Penn"] =asin((sin(j[n-1,"A3Penn"])*j[n-1,"A3FL"])/j[n,"A3FL"]);
    j[n,"A3Fmax"] =x$A3FmaxNoPenn*cos(j[n,"A3Penn"]);
    A3ffv.c <- c(A3ffv.c,j[n,"A3v"])
    j[n,"A3Ffv"] =vl_ffv(V=mean(tail(A3ffv.c,tail.n)),Vmax=x$Vmax,k=x$k)
    #j[n,"A3Ffv"] =(Vmax-mean(j[range.n,"A3v"]))/(Vmax+mean(j[range.n,"A3v"])*G);
    ifelse(j[n,"A3FL"] >x$A3FLOpt,j[n,"A3Ffl"]<- -6.25*(j[n,"A3FL"]/x$A3FLOpt)^2+12.5*(j[n,"A3FL"]/x$A3FLOpt)-5.25,j[n,"A3Ffl"]  <- 1)

    #ifelse(j[n,"A3FL"]/A3FLOpt >= 0.5 & j[n,"A3FL"]/A3FLOpt <= 1.56, j[n,"A3Ffl"]<- predict.lm(fl.poly,data.frame(p.fl=j[n,"A3FL"]/A3FLOpt)),j[n,"A3Ffl"]  <- 0)# from Porro et al. (2002)

    if(j[n,"Time"]<x$ActRiseTime){j[n,"A3Fact"]=0.5-0.5*cos(pi*j[n,"Time"]/x$ActRiseTime)} else{j[n,"A3Fact"]=1};

    ifelse(j[n,"A3ML"]>x$A3MLOpt,j[n,"A3Fpar"] <- x$A3FmaxClosed*exp(2*log(1+x$MaxFpar)*(j[n,"A3ML"]/x$A3MLOpt-1))-x$A3FmaxClosed,j[n,"A3Fpar"] <- 0)
    j[n,"A3Force"] =j[n,"A3Fmax"] *j[n,"A3Ffv"]*ifelse(j[n,"A3Ffl"]<0,0,j[n,"A3Ffl"])*j[n,"A3Fact"]+j[n,"A3Fpar"];

    j[n,"A3Work"]=j[n,"A3Force"]*	j[n,"A3dML"]/x$Scale
    j[n,"A3Power"]=j[n,"A3Work"]/x$TimeStepSecs
    j[n,"A3relL"]=j[n,"A3ML"]/x$A2MLOpt


    #************* Aw force ********
    j[n,"AwdML"] =j[n-1,"AwML"]-j[n,"AwML"];
    j[n,"Awv"] =j[n,"AwdML"]/x$AwMLOpen/x$TimeStepSecs;
    j[n,"AwFL"] =sqrt((sin(j[n-1,"AwPenn"])*j[n-1,"AwFL"])^2+((cos(j[n-1,"AwPenn"])*j[n-1,"AwFL"])-j[n,"AwdML"])^2);
    j[n,"AwPenn"] =asin((sin(j[n-1,"AwPenn"])*j[n-1,"AwFL"])/j[n,"AwFL"]);
    j[n,"AwFmax"] =x$AwFmaxNoPenn*cos(j[n,"AwPenn"]);
    Awffv.c <- c(Awffv.c,j[n,"Awv"])
    j[n,"AwFfv"] =vl_ffv(V=mean(tail(Awffv.c,tail.n)),Vmax=x$Vmax,k=x$k)
    #j[n,"AwFfv"] =(Vmax-mean(j[range.n,"Awv"]))/(Vmax+mean(j[range.n,"Awv"])*G);
    ifelse(j[n,"AwFL"] >x$AwFLOpt,j[n,"AwFfl"]<- -6.25*(j[n,"AwFL"]/x$AwFLOpt)^2+12.5*(j[n,"AwFL"]/x$AwFLOpt)-5.25,j[n,"AwFfl"]  <- 1) # length-tension relationship from Van Wassenbergh


    if(j[n,"Time"]<x$ActRiseTime){j[n,"AwFact"]=0.5-0.5*cos(pi*j[n,"Time"]/x$ActRiseTime)} else{j[n,"AwFact"]=1};
    ifelse(j[n,"AwML"]>x$AwMLOpt,j[n,"AwFpar"] <- x$AwFmaxClosed*exp(2*log(1+x$MaxFpar)*(j[n,"AwML"]/x$AwMLOpt-1))-x$AwFmaxClosed,j[n,"AwFpar"] <- 0)

    j[n,"AwForce"] =j[n,"AwFmax"] *j[n,"AwFfv"]*ifelse(j[n,"AwFfl"]<0,0,j[n,"AwFfl"])*j[n,"AwFact"]+j[n,"AwFpar"];

    j[n,"AwWork"]=j[n,"AwForce"]*	j[n,"AwdML"]/x$Scale
    j[n,"AwPower"]=j[n,"AwWork"]/x$TimeStepSecs
    j[n,"AwrelL"]=j[n,"AwML"]/x$A2MLOpt


    #************** Torques 	#Scaled to N*m  *****************



    #remove negative muscle forces
    if(j[n,"A1Force"]<0) j[n,"A1Force"] <- 1e-15
    if(j[n,"A3Force"]<0) j[n,"A3Force"] <- 1e-15
    if(j[n,"AwForce"]<0) j[n,"AwForce"] <- 1e-15

    #************** Aw
    j[n,"AwTorque"]=j[n,"AwForce"]*abs(sin(j[n,"AwTheta"]))*x$AwLi/x$Scale
    #************** A3
    j[n,"A3Torque"]=j[n,"A3Force"]*abs(sin(j[n,"A3Theta"]))*x$A3Li/x$Scale

    #************** A1 Tendon input *****************
    #* insertion angle of post seg of a
    #*
    #*
    j[n,"A1TendonTheta"]=cos_ang(x$A1Li,x$A1TendonPostSeg,dist_2d(x$joint.p,A1Ins))

    # calculate tension in post seg of A1 tendon
    A1_tension <- function(orig_p,ins_p,inlev_p,ant_p,f){
      alpha <- abs(cos_ang(dist_2d(inlev_p,ins_p),dist_2d(orig_p,ins_p),dist_2d(inlev_p,orig_p)))
      if(alpha>pi/2) alpha <- alpha-pi/2
      beta <- cos_ang(dist_2d(orig_p,ins_p),dist_2d(ant_p,ins_p),dist_2d(orig_p,ant_p)) %>% abs
      if(beta>pi/2) beta <- beta-pi/2
      tension <- f/(cos(alpha)*sin(beta)/cos(beta)+sin(alpha))
      return(tension)
    }


    j[n,"A1TendonForce"] <- A1_tension(x$A1Or,A1Ins,A1inlever, x$A1AntSeg,f=j[n,"A1Force"])
    j[n,"A1Torque"]=j[n,"A1TendonForce"]*sin(j[n,"A1TendonTheta"]);
    j[n,"A1Fout"]= j[n,"A1Torque"]/(x$A1Lo/x$Scale)
    j[n,"A2Fout"]= j[n,"A2Torque"]/(x$A2Lo/x$Scale)
    #************* performance
    j[n,"AwFout"]= j[n,"AwTorque"]/(x$AwLo/x$Scale) #Aw bite force at tip of the jaw
    j[n,"JawVelTip"]=j[n,"AngVel"]*(x$A2Lo) #in cm/s
    j[n,"StaticBite"]=j[n,"A1Fout"]+j[n,"A3Fout"]+j[n,"AwFout"] #Total bite force at tip of jaw
    j[n,"EMA"]=j[n,"StaticBite"]/(j[n,"AwForce"]+j[n,"A1Force"]+j[n,"A3Force"]) #Total bite force at tip of jaw


    plot.ls <- list(joint=x$joint.p,A1Or=x$A1Or,A3Or=x$A3Or,JawTip=JawTip,A3inlever=A3inlever,AwOr=x$AwIns,Awinlever=Awinlever,A1Ins=A1Ins,AntSeg=x$A1AntSeg,A1inlever=A1inlever)

    plot.df<- as.data.frame(do.call(rbind,plot.ls)) %>% mutate(Time=j[n,"Time"],n=n,JawAng=j[n,"JawAng"])
    colnames(plot.df) <- c("x","y","Time","n","JawAng")
    geom.l[[n]] <- plot.df %>% mutate(name=rownames(plot.df),Time=round(Time,5))

    if(n>2 & restart==TRUE){j[2,] <- j[3,]
    n <- 2
    j[n,"Time"] <- 0
    j[n,"Iteration"] <- 1

    restart <- FALSE
    }
  }

  j <- j %>% filter(JawAng>=x$ThetaClosed) %>% mutate(Time=round(Time,5))


  #if(min(j %>% filter(JawAng>0) %>% pull(JawAng))>x$ThetaClosed) warning("Max number of interations has been reached and jaw is not closed.")
closed <- x$ThetaClosed
  geom.df <- as.data.frame(do.call(rbind,geom.l)) %>% filter(JawAng>=closed) %>% mutate(Time=round(Time,5))


  return(list(j=j,pos=geom.df,open=open))
}


# p.anim <- geom.df %>% ggplot()+geom_point(aes(x,y,col=name))+theme_void()+transition_time(as.integer(Time))+labs(title = "Time: {frame_time} ms")+xlim(range(geom.df$x))+ylim(range(geom.df$y))+ coord_fixed()
#
#
#
#   p <- gganimate::animate(p.anim, fps=10)
#   print(p)




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

    ###########  A1
    A1ML=x$A1MLOpen,
    A1dML=0,
    A1v=0,
    A1Penn=0,
    A1FL=x$A1FLOpen, #Instantaneous FL
    A1Fmax=0,
    A1Ffv=0, #Force-velocity factor
    A1Ffl=1, #Force-length factor
    A1Fact=0, #Activation rise-time factor
    A1Fpar=ifelse(x$A1MLOpen>x$A1MLOpt,x$A1Fmax*exp(2*log(1+x$MaxFpar)*(x$A1MLOpen/x$A1MLOpt-1))-x$A1Fmax,0), #Parallel elastic force
    A1Force=0, #Force of the A1 division
    A1Theta=x$A1MuscleThetaOpen, #angle of A1 relative to tendon
    A1TendonTorque=0, #Torque imparted by A1 on tendon;
    A1Torque=0, #Torque imparted by A1 on lower jaw;
    A1Fout=0,
    A1InsTheta=x$A1InsThetaOpen, #angle post seg of A1 ins on lower jaw
    A1TendonAntSeg=x$A1TendonAntSeg,
    A1TendonAntSeg=x$A1TendonPostSeg,

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

    #Aw vectors
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

lm_dat<-system.file("extdata", "LMBass_1_muscle_data.csv", package = "JawCloseR") %>% read.csv() %>% rownames_to_column(var="spec.n")

sm_dat<-system.file("extdata", "SMBass_1_muscle_data.csv", package = "JawCloseR") %>% read.csv() %>% rownames_to_column(var="spec.n")

ps_dat<-system.file("extdata", "Pseed_1_muscle_data.csv", package = "JawCloseR") %>% read.csv() %>% rownames_to_column(var="spec.n")


ps <- ps_dat%>% mutate_at(.vars = vars(ends_with("Mass")),function(x)x/1) %>% mutate(Gape=55)%>% start_old(pars = list(Pressure=4000,Vmax=15,CdPlate=2.5),config=4) %>% close_old_config4()
ps$j$Time %>% max()

x=lm_dat
review_data(lm_dat)

lm <- lm_dat%>% mutate_at(.vars = vars(ends_with("Mass")),function(x)x/1) %>% mutate(Gape=80)%>% start_old(pars = list(Pressure=4000,Vmax=15,CdPlate=2.5),config=4) %>% close_old_config4()
lm$j$Time %>% max()

lm$j %>%
  ggplot(aes(JawAng %>% deg,A1v))+geom_point()+scale_x_reverse()


sm <- sm_dat%>% mutate_at(.vars = vars(ends_with("Mass")),function(x)x/1)%>% mutate(Gape=60)%>%start_old(pars = list(Pressure=4000,Vmax=15,CdPlate=2.5),config=4) %>% close_old_config4()

lm$j %>%
  ggplot(aes(JawAng %>% deg,DragT))+geom_point()+scale_x_reverse()


sm$j$Time %>% max()
cl.l <- list()

for(s in dat$spec.n %>% unique){
  for(c in 1:3 ){
    cl.i <- start_old(spec.n=s,config=c,pars=list(Vmax=30)) %>%
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
  mutate(per_closed=1-plyr::round_any(JawAng/rad(110),0.02)) %>%
  mutate(per_closed=round(per_closed,2)) %>%
  filter(per_closed %in% seq(0,1,0.02)) %>%
  group_by(config,per_closed) %>%
  summarise(A2v.m=mean(A2v),A2v.se=se(A2v)) %>%
  ggplot(aes(per_closed,A2v.m,col=as.factor(config)))+geom_point()+geom_errorbar(aes(x=per_closed,ymin=A2v.m-A2v.se,ymax=A2v.m+A2v.se))
}




