
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
#' cs_dat %>%
#'    start_sim(spec.n=2) %>%
#'    print()


start_simOLD <- function(x=NULL,spec.n=1,config=1,pars=NULL)	{

 if(!is.null(x)){ w <- tryCatch(review_data(x),warning = function(w) w)

  if("warning" %in% class(w)) stop("Input data passed to `x` are broken. Thry may have missing variables, non-numeric data, etc. Please use `reveiew_data()`.")
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


 # MuscleData<-system.file("extdata", "SMBass_1_muscle_data.csv", package = "JawCloseR") %>% read.csv() %>% rownames_to_column(var="spec.n")
  SpecimenNumber=spec.n;
  #Set morphological variables
  SL=MuscleData[SpecimenNumber,"SL"];
  JL=MuscleData[SpecimenNumber,"JL"];

  #A1 geometry and mass
  A1OrIns=MuscleData[SpecimenNumber,"A1OrIns"];
  A1OrJoint=MuscleData[SpecimenNumber,"A1OrJoint"]; #change w/ SpecAnge
  A1Tendon=MuscleData[SpecimenNumber,"A1Tendon"]
  A1Li=MuscleData[SpecimenNumber,"A1Li"];
  A1Lo=MuscleData[SpecimenNumber,"A1Lo"];
  A1AvgFL=MuscleData[SpecimenNumber,"A1AvgFL"];#change w/ SpecAnge
  A1AvgPenn=(pi/180)*MuscleData[SpecimenNumber,"A1AvgPenn"]; #change w/ SpecAng
  A1TendonAntSeg <- MuscleData[SpecimenNumber,"A1TendonAntSeg"]
  A1TendonPostSeg <- MuscleData[SpecimenNumber,"A1TendonPostSeg"]
  A1TendonOrIn <- MuscleData[SpecimenNumber,"A1TendonOrIn"]
  A1Mass <- MuscleData[SpecimenNumber,"A1Mass"]
  A1InsTheta <- MuscleData[SpecimenNumber,"A1InsTheta"] #insertion angle of post seg of A1 Tendon
  A1TendonTheta <- MuscleData[SpecimenNumber,"A1TendonTheta"]

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
  #MandWidth <- MuscleData[SpecimenNumber,"MandWidth"]
  # MandDepth <- MuscleData[SpecimenNumber,"MandDepth"]
  # SymphLength <- MuscleData[SpecimenNumber,"SymphLength"]
  Gape <- (pi/180)*MuscleData[SpecimenNumber,"Gape"]
  SpecAng<- (pi/180)*MuscleData[SpecimenNumber,"SpecAng"] #Jaw open angle of specimen in


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
  A1MA <- A1Li/A1Lo
  A3MA <- A3Li/A3Lo
  A2MA <- A2Li/A2Lo
  AwMA <- AwLi/AwLo

  # #Calculate Fmax (N) in closed position
  A1FmaxClosed <-MaxIso/10*A1PCSAclosed
  A3FmaxClosed <-MaxIso/10*A3PCSAclosed
  A2FmaxClosed <-MaxIso/10*A2PCSAclosed
  AwFmaxClosed <-MaxIso/10*AwPCSAclosed


  #*************** Muscle geometry *******************

  #Angle of line between origin-joint and lower jaw when closed
  A1ThetaJointClosed <- acos((A1Li^2+A1OrJoint^2-A1OrIns^2)/(2*A1Li*A1OrJoint))
  A3ThetaJointClosed <- acos((A3Li^2+A3OrJoint^2-A3OrIns^2)/(2*A3Li*A3OrJoint))
  A2ThetaJointClosed <- acos((A2Li^2+A2OrJoint^2-A2OrIns^2)/(2*A2Li*A2OrJoint))
  AwThetaJointClosed <- acos((AwLi^2+AwOrJoint^2-AwOrIns^2)/(2*AwLi*AwOrJoint))

  #Angle of line between origin-joint and lower jaw when open
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
  start <- rad(0)
  open.ang <- -1*(start+Gape-SpecAng)
  A2inlever <- point_ang_r(joint.p,A2Li,SpecAng-start) #meckelian tendo ins point shared by A2/A3
  A1inlever <- point_ang_r(joint.p,A1Li,SpecAng-start)
  JawTip <- point_ang_r(joint.p,JL,SpecAng-start)
  A1Or <- point_ang_r(joint.p,A1OrJoint,A1ThetaJointClosed-start)
  A2Or <- point_ang_r(joint.p,A2OrJoint,A2ThetaJointClosed-start)
  A3Or <- point_ang_r(joint.p,A3OrJoint,A3ThetaJointClosed-start)
  #figure our Aw angles
  AwOr <- point_ang_r(joint.p,AwLi,SpecAng-start)
  #insertion point on tendon based on muscle theta and orig-ins length (so distal tendon)
  TendonLengthClosed <- sin(AwMuscleThetaClosed)*AwOrIns

  #force Aw to insert on middle of tendon according to muscle theta
  AwIns<- point_ang_r(AwOr,AwLi-A2Li,-rad(180)-AwMuscleThetaClosed-start)

  dist
  #A1 post segment insertion angle


  A1Ins <- point_ang_r(A1inlever,A1TendonPostSeg,(rad(180))-A1InsTheta-start) #insertion on A1 tendon (will be dynamics)



  #A1 ant segment insertion on maxillary

  A1AntSeg <- point_ang_r(A1Ins,MuscleData$A1TendonAntSeg,A1TendonTheta-A1InsTheta)


  geom.df.open <- t(data.frame(joint.p,A2Or,A3Or,JawTip,A1inlever,A1Ins,A1Or,A1AntSeg));colnames(geom.df.open) <- c("x","y")
  geom.df.open.melt <- as.data.frame(geom.df.open) %>% mutate(name=rownames( geom.df.open))
  p.open <- ggplot(data=geom.df.open.melt,aes(x=x,y=y))+geom_point(aes(colour=name),size=6)+coord_fixed()
p.open

  ################################################################

  #Length of Muscle (origin to insertion on tendon + proximal tendon) when closed, based on cartesian coords
  A1MLClosed <- dist_2d(A1Ins,A1Or)

  A3MLClosed <- dist_2d(AwIns,A3Or)
  A2MLClosed<- dist_2d(AwIns,A2Or)
  AwMLClosed<- dist_2d(AwIns,AwOr)



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

  A2inlever.open <- point_ang_r(joint.p,A2Li,open.ang)
  A1inlever.open <- point_ang_r(joint.p,A1Li,open.ang)

  JawTip.open <- point_ang_r(joint.p,JL,open.ang)
  #figure out Aw angles
  AwOr.open <- point_ang_r(joint.p,AwLi,open.ang)

  #force Aw to insert on middle of tendon according to muscle theta
  #mean of two muscle insertions scaled around 180
  TendonThetaOpen <-mean(c(A2MuscleThetaOpen,A3MuscleThetaOpen))
  AwIns.open<- point_ang_r(A2inlever.open,TendonLengthClosed,(rad(180)-(Gape-SpecAng))-TendonThetaOpen-start) #really the origin and the nexus



  #Dist between A1 tendon ant and post ins/or when open
  A1TendonOrIn.open <- dist_2d(A1inlever.open,A1AntSeg)

  #Angle between ant and post seg of A1 tendon
  A1TendonTheta.open<- cos_ang(A1TendonPostSeg,A1TendonAntSeg,A1TendonOrIn.open)

  #Adjust A1 Ant segment to that it's just long enough at open angle to make it straight, Start sims from there

  ## Recalc angle
  A1TendonTheta.open<- cos_ang(A1TendonPostSeg,A1TendonAntSeg,A1TendonOrIn.open)

  A1AntSeg.open <- point_ang_r(A1Ins.open,A1TendonAntSeg,A1TendonTheta.open-A1InsTheta-Gape-SpecAng)

  A1TendonOrIn.open <- dist_2d(A1inlever.open,A1AntSeg.open)

  ### ins angle of post seg of A1 tendon on jaw
  A1InsThetaOpen <- cos_ang(A1Li,A1TendonPostSeg,dist_2d(A1Ins.open,A1inlever.open))

  #open position of A1 ins on jaw
  A1Ins.open<- point_ang_r(A1inlever.open,A1TendonPostSeg,(rad(180)-(Gape-SpecAng))- A1InsThetaOpen -start)

  A1TendonAntSeg <- A1TendonOrIn.open-A1TendonPostSeg


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
      A1OrIns = A1OrIns,
      A1OrInsOpen = A1OrInsOpen,
      A1OrJoint = A1OrJoint,
      A1PCSAclosed = A1PCSAclosed,
      A1PennOpen = A1PennOpen,
      A1TendonTheta = A1Ins.open,
      A1ThetaJointClosed = A1ThetaJointClosed,
      A1ThetaJointOpen = A1ThetaJointOpen,
      A1TendonAntSeg=A1TendonAntSeg,
      A1TendonPostSeg=A1TendonPostSeg,
      A1InsThetaOpen= A1InsThetaOpen,
      A1AntSeg=A1AntSeg.open,
      CdPlate = CdPlate,
      config=config,
      G = G,
      Gape = Gape,
      Ijaw = Ijaw,
      Itot=Itot,
      A2inlever.open = A2inlever.open,
      A2inlever = A2inlever,
      A1inlever=A1inlever.open,
      A1Ins=A1Ins.open,
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

#' x=start.l

close_jawOLD <- function(x,MaxIts=1000,progress=F){

  if(any(sapply(x, function(x) length(x)==0))) stop("Start values passed to `x` contain 0-length variables.")

  open <- x
  j.var <- var_df(x=x)
  j <- data.frame(matrix(0,ncol=ncol(j.var),nrow=MaxIts))
  j[1,] <- j.var
  colnames(j) <- colnames(j.var)

  tendon.c <-NULL

  pb = txtProgressBar(min = 2, max = MaxIts, initial = 2,style=3)

  if(progress) pdf(file=paste(graph.dir,"/geom_closed.pdf",sep=""),width=7,height=7)


  geom.l <- list()

  for(n in 2:MaxIts){
    if(j[n-1,"JawAng"]<=x$ThetaClosed) break


    setTxtProgressBar(pb,n)
    j[n,"Time"]<- j[n-1,"Time"]+x$TimeStep
    j[n,"Iteration"]<- n-1
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


    ######muscle lengths
    m.max <- ifelse(n==2,1,which.max(abs(c(j[n-1,"A2TendonTorque"],j[n-1,"A3TendonTorque"],j[n-1,"AwTendonTorque"]))))


    ### reposition according to sag and equilibrium imposed by Aw on resultant line.
    nexus <- point_ang_r(A2inlever,j[1,"TendonLength"],pi+open.angle-tendon.theta)

    j[n,"AwML"] <- dist_2d(nexus,AwOr)
    j[n,"A2ML"] <-dist_2d(nexus,x$A2Or)
    j[n,"A3ML"] <-dist_2d(nexus,x$A3Or)


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
    tendon.theta <- mean(tail(tendon.c,10))

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

    if(n>=2){ A1Ins <-x$A1Ins

    A1TendonOrIn <- dist_2d(A1inlever,x$A1AntSeg)
    }




   A1_ins<- function(pt1,pt2,l,r){
      il=A1inlever
      ant_ins=x$A1AntSeg
      pt=A1Ins

      theta1=atan2(diff(c(il[2],ant_ins[2])),diff(c(il[1],ant_ins[1])))
      theta2=cos_ang(x$A1TendonPostSeg,A1TendonOrIn,x$A1TendonAntSeg)
      point_ang_r(p = il,theta = theta1,r = x$A1TendonPostSeg)

      dist_2d(il,A1Ins)

        x$A1TendonAntSeg
      r=x$A1TendonPostSeg
      pt1=x$A1AntSeg
      pt2=A1inlever
      o=dist_2d(pt1,pt2)
      y.pt3=(o^2+l^2-r^2)/(2*o)
      x.pt3=sqrt(l^2-y.pt3^2)

      x.pt3 <- pt2[1]+x.pt3
      y.pt3 <- y.pt3-pt1[2]
      return(c(x.pt3,y.pt3))

    }

    point_sides(x$A1TendonAntSeg,x$A1TendonPostSeg,A1TendonOrIns)
    A1TendonOrIns <- dist_2d(A1inlever,x$A1AntSeg)
    A1TendonTheta <- cos_ang(x$A1TendonAntSeg,x$A1TendonPostSeg)


    j[n,"A1ML"] <- dist_2d(,x$A1Or)
    #}

    j[n,"AwTheta"] <- cos_ang(dist_2d(A2inlever,AwOr),j[n,"AwML"],j[1,"TendonLength"])


    #for smoothing FFV over prev 2 ms
    if(n>21){range.n <- (n-20):(n-1)}else{range.n <- n-1}

    # nexus is bound by lines of action of A2 and A3
    j[n,"TendonLength"] <-j[1,"TendonLength"]

    j[n,"TendonTheta"]=cos_ang(x$A2Li,j[1,"TendonLength"],dist_2d(x$joint.p,nexus))

    #************* A3 force ********
    j[n,"A3dML"] =j[n-1,"A3ML"]-j[n,"A3ML"];x$imeStepSecs
    j[n,"A3v"] =j[n,"A3dML"]/x$A3MLOpen/x$TimeStepSecs;
    j[n,"A3FL"] =sqrt((sin(j[n-1,"A3Penn"])*j[n-1,"A3FL"])^2+((cos(j[n-1,"A3Penn"])*j[n-1,"A3FL"])-j[n,"A3dML"])^2);
    j[n,"A3Penn"] =asin((sin(j[n-1,"A3Penn"])*j[n-1,"A3FL"])/j[n,"A3FL"]);
    j[n,"A3Fmax"] =x$A3FmaxNoPenn*cos(j[n,"A3Penn"]);
    j[n,"A3Ffv"] =vl_ffv(V=mean(j[range.n,"A3v"]))
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
    j[n,"A2Ffv"] =vl_ffv(V=mean(j[range.n,"A2v"]));
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
    j[n,"AwFfv"] =vl_ffv(V=mean(j[range.n,"Awv"]))
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

    #for plotting
    plot.ls <- list(joint.p=x$joint.p,A2Or=x$A2Or,A3Or=x$A3Or,JawTip=JawTip,inlever=A2inlever,AwIns=AwOr,nexus=nexus,resOr=resOr)

   plot.df<- as.data.frame(do.call(rbind,plot.ls)) %>% mutate(Time=j[n,"Time"],JawAng=j[n,"JawAng"])
   colnames(plot.df) <- c("x","y","Time","JawAng")
    geom.l[[n]] <- plot.df %>% mutate(name=rownames(plot.df),Time=round(Time,5))
  }

  j <- j %>% filter(JawAng>=x$ThetaClosed) %>% mutate(Time=round(Time,5))




  geom.df <- as.data.frame(do.call(rbind,geom.l))

  #geom.out[[paste(n)]] <- data.frame(ms=j[n,"Time"],ang=j[n,"JawAng"],var=unlist(plot.ls))


  return(list(open=open,j=j,pos=geom.df))
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
jaw_varsOLD <- function(){
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




