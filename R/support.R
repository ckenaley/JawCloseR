library(cowplot)
library(tidyverse)

#' @title  Radians to degrees
#' @param x numeric, degree value.
#' @return A single value of the angle in degrees
#' @export
deg <- function(x){180/pi*x}

#' @title  Degrees to radians
#' @param x numeric, radian value.
#' @return A single value of the angle in radians.
#' @export
rad <- function(x){pi/180*x}


#' @title  Find coordinates of point some path away from origin

#' @description  Internal function used to find x,y positions of paths given and origin, radius, and angle.
#'
#' @param p numeric, the origin. A vector or two values, x and y of a point in cartesian space.
#' @param r numeric, the path length or radius.
#' @param theta numeric, the angle of the path.
#' @return A vector of length 2 represent x and y position
#' @export
#'
point_ang_r <- function(p,r,theta){x.2 <-p[1]+r*cos(theta); y.2 <- p[2]+r*sin(theta);return(c(x.2,y.2))}


# #add line between two points labelled with first col
# line_2pt <- function(a,b,dat,col="black"){geom_line(data=subset(dat,dat[,3]==a|dat[,3]==b),aes(x,y,group = 1),col=col,size=2,alpha=0.7)}


#' @title  Compute angle between two segments of a triangle using law of cosines.

#' @description  Computes angle between two segments of a triangle using law of cosines.
#'
#' @param l numeric,  length of segment to the left.
#' @param r length of segment to the right.
#' @param o length of opposite segment.
#' @return A single value of the angle in radians.
#' @export
#'
cos_ang <- function(l,r,o) {
  if(abs((l+r)-o<1e-10)){return(pi)}else{
    a <- acos((o^2-l^2-r^2)/(-2*l*r))
    return(a)}
}

# #law of cosines to find length of opposite given angle, and lengths to left and right
#' @title  Compute length of as segment using law of
#' #' @description  Computes length of a segment given the length of two adjacent segments and their angle using law of cosines.
#' @param l numeric,  length of segment to the left.
#' @param r length of segment to the right.
#' @param and numeric, the angle shared by the segments.
#' @return A single value the length.
#' @export

cos_side <- function(l,r,ang){sqrt(l^2+ r^2-2*l*r* cos(ang))}

#' @title  Compute distance between two points in Cartesian space

#' @description  Computes distance between two points in Cartesian space using simple trigonometry functions
#'
#' @param pt1 numeric, vector of length 2 representing x,y position of the first point.
#' @param pt2 numeric, vector of length 2 representing x,y position of the second point.
#' @return A single value of the length.
#' @export
#'
dist_2d <- function(pt1,pt2){sqrt((pt2[1]-pt1[1])^2+(pt1[2]-pt2[2])^2)}


# #2d distance between point (pt0) and line defined by two points, neg value reports to positive values for x and y relative to line
#
# line.point.2d <- function(pt1,pt2,pt0){
#   ((pt2[2]-pt1[2])*pt0[1]-(pt2[1]-pt1[1])*pt0[2]+pt2[1]*pt1[2]-pt2[2]*pt1[1])/sqrt((pt2[2]-pt1[2])^2+(pt2[1]-pt1[1])^2)}

#st error
std <- function(x) sd(x)/sqrt(length(x))

#' @title  Compute the force-velocity factor from Van Leuwen et al. (1990) using a modified Hill model

#' @description  Compute the force-velocity factor from \insertCite{van1990function;textual}{JawCloseR}. Used internally by \code{close_jaw} to scale muscle-force output according to shortening velocity.
#'
#' @param V numeric, he instantaneous contraction velocity (in muscle lengths, ML, per second).
#' @param k numeric, the constant determining the shape of the Hill curve.
#' @param Vmax numeric, the maximum contraction velocity (in ML per second).
#' @return A single numeric value of the factor.
#'
#' @details  \loadmathjax This function returns the force-velocity factor given by a modified Hill  model \insertCite{hill1938heat}{JawCloseR} under both concentric \eqn{V>0} and eccentric \eqn{V<0} conditions published by Van Leuwen et al. 1990:
#'
#' \mjseqn{ F_{FV}=\cfrac{V_{max}-V}{V_{max}+V/k}, \text{when } V>0}
#'
#' \mjseqn{ F_{FV}=1.8-0.8 \times \cfrac{1+V/V_{max}}{1-7.56 \times V/(k \times V_{max})}, \text{when }  V<0}
#'
#' The defualt value for \eqn{k} or 0.25 corresponds to vertebrate white muscle according to \insertCite{alexander2003principles;textual}{JawCloseR}.
#'
#' @import mathjaxr
#' @importFrom Rdpack reprompt
#' @export
#'
#' @references
#' \insertAllCited{}
#'
vl_ffv <- function(V,k=0.25,Vmax=15){
  sapply(V,function(x){
    Fv.con <-  (Vmax-x)/(Vmax+x/k)
    #Fv.con <- (1-x/Vmax)/(1+x/(k*Vmax))
    Fv.ecc <- 1.6-0.6*(1+x/Vmax)/(1-7.56*x/(k*Vmax))
    if(x<0) Fv <- Fv.ecc
    if(x>0) Fv <- Fv.con
    if(x==0) Fv <- 1
    #if(Fv<0) Fv <-  1e-5
    return(Fv)
  })
}

# #plot t see v-tension relationship
# v=seq(-20,25,0.01)
# ffv <- vl_ffv(V=v)
# min(ffv)
# #plot(v,ffv)

# #parallel passive force-not used, but could be
# Fpar <- function(ml,ml.opt=1.0,fmax.par=0.2){
#   exp(2*log(1+fmax.par)*(ml/ml.opt-1))-1
# }

# #plot t see length-tension relationshi
# l <- seq(1.2,0.75,-0.01)
# f <- unlist(sapply(l, function(x) Fpar(ml=x)))
# #plot(l,f)

#' @title  Establish a simulation data.frame

#' @description  Internal function used to establish  a \code{data.frame} where results are stored by \code{close_jaw}.
#'
#' @param x a list obtained by \code{start_sim}
#' @export
#'
#'
var_df<- function(x){
d <- data.frame(
  Start=x$Start,
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
  A2Ffl=0, #Force-length factor
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
  return(d)
}


#' @title  Review input data

#' @description  Review input data, making sure it is appropriate for \code{start_sim}.
#' @details This functions tests a \code{data.frame} a user intends to pass to \code{start_sim()}, evaluating if the required variables are present and indicating if the data pass or what variables are missing.
#'
#' @param x a \code{data.frame} a user intends to pass to \code{start_sim()}
#' @param missing logical, should missing variables be returned
#' @return If \code{missing=TRUE}, missing variables are returned as vector.
#' @export
#'
#' @examples
#' library(tidyverse)
#' dat <- system.file("extdata", "csloani.csv", package = "JawCloseR") %>% read.csv()
#' dat <- dat[-4]
#' miss_vars <- dat %>% review_data(missing=TRUE)
#' print(miss_vars)
#'
review_data <- function(x,missing=FALSE){
  df.n <- jaw_vars() %>% pull(name)
  n <- colnames(x)
  extra <- n[!n%in%df.n] %>% unlist
  needed <- df.n[!df.n%in%n]

  no.num <- names(x)[!sapply(x, is.numeric)]

  no.num.needed <- no.num[no.num%in%df.n]

  if(!is.data.frame(x)) stop("'x' and data passed to `start_sim()` must be a data.frame")

  if(any(no.num%in%needed)) warning(paste0("The following required variables/columns are in the input data, but not numeric: ",paste(needed,collapse=", ")))

  if(length(needed)>0) warning(paste0("The following required variables/columns are not in the input data: ",paste(needed,collapse=", ")))

  if(length(extra)>0) message(paste0("The following variables/columns are in the input data, but are unneeded and will be ignored: ",paste(extra,collapse=", ")))

  if(missing) return(needed)
}




#' @title  Plot a start of jaw closing simulation

#' @description  Plots the output of important positions established by \code{start_sim}.
#' @details This function returns a plot of the origin and insertions positions as well as the jaw based on the output of data passed through \code{start_sim}.
#'
#' @param x a \code{data.frame} a list returned by \code{start_sim()}
#' @return a ggplot object
#' @export
#'
#' @examples
#' start_sim() %>%
#'    plot_open()
#'

plot_open <- function(x){

  geom.df.open <- t(data.frame(joint=x$joint.p,A2Or=x$A2Or,A3Or=x$A3Or,JawTip=x$JawTip.open,Inlever=x$inlever.p.open,AwIns=x$AwOr.open,nexus=x$AwIns.open)) %>% data.frame
  colnames(geom.df.open) <- c("x","y")
  geom.df.open$name <- rownames(geom.df.open)

  p <- geom.df.open %>%
    ggplot(aes(x=x,y=y))+line_2pt("joint","JawTip",dat = geom.df.open)+
    line_2pt("A2Or","nexus",dat = geom.df.open,col="gray")+
    line_2pt("A3Or","nexus",dat = geom.df.open,col="gray")+
    line_2pt("AwIns","nexus",dat = geom.df.open,col="gray")+
    line_2pt("AwIns","nexus",dat = geom.df.open,col="gray")+
    geom_point(aes(colour=name),size=3,alpha=0.7)+
    theme_minimal()
  return(p)
}





#' @title  Animates the simulated closing of a jaw

#' @description  Plots the output of important positions established simulated by \code{close_jaw()} throughout adduction.
#' @details This function returns a plot of the origin and insertions positions as well as the jaw based on the output of data passed through \code{start_sim}.
#'
#' @param x a \code{list} returned by \code{close_jaw()}.
#' @param addjaw logical should a cartoon of a jaw and head be plotted
#' @param jawsp character, head and jaw species. See \code{details}.
#' @param render should plot be rendered with \code{gganimate::animate()}. See \code{details}.
#' @param fps numeric, frames per second controlling playback speed of animation.
#' @param size numeric, width and height in pixels of animation.
#' @param shadow logical, should shadow of previous positions be plotted.
#' @param inset logical, should muscle data be printed to animation images.
#' @return a gif \code{magick-image} object if \code{render=TRUE}, a \code{gganim} object if \code{render=FALSE}. See \code{details}.
#'
#' @details Plotted data include positions origins, insertions, and lines of action/muscle lengths for the adductor muscle divisions, the Meckelian tendon, and jaw. If \code{addjaw=TRUE} a cartoon of  the head and jaw of a limited number of species indicated by \code{jawsp} is also plotted beneath these data. Current species option include a viperfish (\code{jawsp="Viperfish"}), a largemouth bass (\code{jawp="LMBass"}), wrasse (\code{jawsp="Wrasse"}), or gar (\code{jawsp="Gar"}).
#'
#'
#' If \code{render=TRUE}, a \code{magick-image} with playback speed and a shadow determined by the arguments \code{fps} and \code{shadow}. Users may want more flexibility in controlling \code{gganimate} options and therefore may use \code{render=FALSE} to return a \code{gganim} object, add options and pass it to \code{gganimate:animate()}. The data passed from \code{x} are identical in both cases except that muscle data are not inset.
#'
#' If \code{inset=TRUE}, instantaneous values for muscle velocity and power are included in the animation.
#'
#'
#'
#' @export
#' @import gganimate viridis cowplot gifski
#'
#' @examples
#'
#' #without a jaw
#' start_sim(config=3) %>%
#'    close_jaw() %>%
#'     {. ->> cs } %>%
#'    animate_close(addjaw=TRUE,jawsp="LMBass")
#'
#'#with a jaw (oddly, a wrasse jaw)
#'
#' cs %>%
#'    animate_close(addjaw=TRUE,jawsp="LMBass")
#'
#'cs_an <- cs %>%
#'    animate_close(render=FALSE)
#'
#'#looks terrible
#'cs_an2 <- cs_an+shadow_wake(0.1, size = 2, alpha = FALSE, colour = 'grey92')
#'animate(cs_an2)
#'
#'x=sm
#'lm %>% animate_close(inset=T)
animate_close <- function(x,addjaw=FALSE,jawsp=NULL,render=TRUE,fps=10,shadow=FALSE,inset=FALSE,size=c(800,800)){

  # if(!is.list(x)|!all(names(x) %in% c("open","j","pos"))) stop("looks like `x` is not a list produced by `close_jaw()`")

  if(!is.numeric(fps)) stop("`fps` must be a numeric value")

  if(!is.logical(shadow)) stop("`shadow must be logical")
  gape <- x$open$Gape
  config <- x$open$config


  if(config>3)  perf <- x$j %>% select(Time,JawAng,StaticBite,A1v,A3v,Awv,ends_with("Power"))

  if(config<=3)  perf <- x$j %>% select(Time,JawAng,StaticBite,A1v,A2v,A3v,Awv,ends_with("Power"))



  ta <- x$j %>% select(Time,dJawAng,JawAng)
  x <- x$pos %>% left_join(ta)


  spec.axis <- x %>% filter(name %in% c("JawTip","joint") & Time==first(Time))


  if(!is.null(jawsp)&addjaw)  {
    jaw.f <- system.file("extdata/jaws", package = "JawCloseR") %>% list.files(full.names = TRUE)

    jaw.n <- gsub(".txt","",basename(jaw.f))
    jaw.n <- jaw.n[!grepl("head",jaw.n)]
    if(!jawsp%in%jaw.n) stop(paste0("`addjaw` name is not available. Please choose from: `",paste0(jaw.n,collapse = "','"),"'"))

    head.n <- basename(jaw.f)[grepl(paste0(jawsp,"_head"),jaw.f)]


    jaw <- system.file("extdata/jaws", paste0(jawsp,".txt"), package = "JawCloseR") %>% read_tsv(col_types = cols()) %>% rownames_to_column(var="n") %>% mutate(n=as.numeric(n)) %>% mutate(part=ifelse(n==last(n),"JawTip","Jaw")) %>% mutate(part=ifelse(n==(last(n)-1),"joint",part))

    head <- system.file("extdata/jaws", head.n, package = "JawCloseR") %>% read_tsv(col_types = cols()) %>% rownames_to_column(var="n") %>% mutate(n=as.numeric(n)) %>% mutate(part=ifelse(n==last(n),"HeadTip","Head")) %>% mutate(part=ifelse(n==(last(n)-1),"joint",part))

    jaw.axis <- jaw %>% filter(part!="Jaw")
    jaw <- jaw %>% filter(part=="Jaw")
    head.axis <- head %>% filter(part!="Head")
    head <- head %>% filter(part=="Head")

    if(nrow(jaw)>200) jaw <- jaw[round(seq(1,nrow(jaw),length.out=200)),]
    if(nrow(head)>200) head <- head[round(seq(1,nrow(head),length.out=200)),]
    jaw.l <- dist_2d(c(jaw.axis[1,]$x,jaw.axis[1,]$y),c(jaw.axis[2,]$x,jaw.axis[2,]$y))

    spec.jaw.l <- dist_2d(c(spec.axis[1,]$x,spec.axis[1,]$y),c(spec.axis[2,]$x,spec.axis[2,]$y))

    jaw.scale <- spec.jaw.l/jaw.l

    jaw <- jaw %>%
      mutate(x=x*jaw.scale,y=y*jaw.scale)
    jaw.axis <- jaw.axis %>%
      mutate(x=x*jaw.scale,y=y*jaw.scale)

    head <- head %>%
      mutate(x=x*jaw.scale,y=y*jaw.scale)
    head.axis <- head.axis %>%
      mutate(x=x*jaw.scale,y=y*jaw.scale)

    jaw <- jaw %>%
      mutate(x=x-jaw.axis[1,]$x,y=y-jaw.axis[1,]$y)

    head <- head %>%
      mutate(x=x-head.axis[1,]$x,y=y-head.axis[1,]$y)


    jaw.axis <- jaw.axis %>%
      mutate(x=x-jaw.axis[1,]$x,y=y-jaw.axis[1,]$y)

    head.axis <- head.axis %>%
      mutate(x=x-head.axis[1,]$x,y=y-head.axis[1,]$y)

    spec.ang <- atan2(spec.axis[2,]$y - spec.axis[1,]$y, spec.axis[2,]$x - spec.axis[1,]$x)
    jaw.ang <- atan2(jaw.axis[2,]$y - jaw.axis[1,]$y, jaw.axis[2,]$x - jaw.axis[1,]$x)
    head.ang <- atan2(head.axis[2,]$y - head.axis[1,]$y, head.axis[2,]$x - head.axis[1,]$x)

    rot.ang <- spec.ang-jaw.ang
    head.rot.ang <- spec.ang-head.ang



    new.jaw <- rotate_coords(jaw[,c("x","y")],rot.ang,c(0,0))
    new.head <- rotate_coords(head[,c("x","y")],head.rot.ang+gape,c(0,0)) %>% data.frame


    jaw.df <-ta %>%
      mutate(dJawAng=abs(JawAng-lag(JawAng))) %>%
      na.omit %>%
      mutate(newAng=cumsum(dJawAng)) %>%
      group_by(Time,JawAng) %>%
      summarize(rotate_coords(new.jaw[,c("x","y")],newAng,c(0,0)))

    head.df <-ta %>%
      mutate(dJawAng=abs(JawAng-lag(JawAng))) %>%
      na.omit %>%
      mutate(newAng=cumsum(dJawAng)) %>%
      group_by(Time,JawAng) %>%
      summarize(rotate_coords(new.head[,c("x","y")],0,c(0,0)))

  }



  if(!addjaw) {jaw.df <- data.frame(x=0,y=0,Time=x$Time)
  head.df <- data.frame(x=0,y=0,Time=x$Time)
  }

  x.lim <- range(c(range(x$x),range(jaw.df$x),range(head.df$x)))
  y.lim <- range(c(range(x$y),range(jaw.df$y),range(head.df$y)))

  grps13 <- data.frame(
    a=c("A2Or","A3Or","AwIns","joint","nexus"),
    b=c("nexus","nexus","nexus","JawTip","inlever"),
    path=c("A2","A3","Aw","jaw","tendon"),
    type=c("muscle","muscle",'muscle',"jaw","Meck. tend"),
    col=c("red","red","red","black","gray")
  ) %>% pivot_longer(cols = a:b,values_to="name",names_to = "position")

  grps4 <- data.frame(
    a=c("A1Or","A3Or","Awinlever","joint","AntSeg","A1Ins"),
    b=c("A1Ins","A3inlever","AwOr","JawTip","A1Ins","A1inlever"),
    path=c("A1","A3","Aw","jaw","AntA1Tendon","PostA1Tendon"),
    type=c("muscle","muscle",'muscle',"jaw","A1 Tendon","A1 Tendon"),
    col=c("red","red","red","black","gray","gray")
  ) %>% pivot_longer(cols = a:b,values_to="name",names_to = "position")

  if(config<4) {grps <- grps13} else{grps <- grps4}

  x <- x %>%
    left_join(grps) %>%
    filter(name!="resOr")



  sh <- shadow_null()

  if(shadow) sh<- gganimate::shadow_mark(alpha=0.05)

  p.anim <- ggplot()+geom_path(data=jaw.df,aes(x,y))+geom_path(data=head.df,aes(x,y))+geom_point(data=x,aes(x,y))+geom_line(data=x,aes(x,y,group=path,color=type))+theme_void(20)+transition_time(Time)+labs(title = "Time: {frame_time %>% round(.,1)} ms")+xlim(x.lim)+ylim(y.lim)+ coord_fixed()

  if(inset){
  vals <- perf %>%
    pivot_longer(c(starts_with("A"),starts_with("Static")))

  max.val <- perf %>%
    na.omit %>%
    summarise_all(range)

  max.val[1,] <- 0

  max <- max.val %>%
    mutate(m=c("min","max")) %>%
    select(-Time,-JawAng) %>%
    pivot_longer(c(starts_with("A"),starts_with("Static"))) %>%
    pivot_wider(names_from=m)

  vals <- vals %>%
    left_join(max) %>%
    mutate(per=value/max)

  p.perf <- vals %>%  ggplot()+
    geom_bar(aes(name,per,fill=per), stat="identity",width=1,position=position_dodge(width = 0.1))+
    scale_fill_viridis(option = "H") +
    theme_minimal(30)+ylim(c(0,1))+theme(legend.position = "none",axis.title = element_blank(),axis.text.x = element_blank())+coord_flip()+transition_time(Time)

  perf.gif <- animate(p.perf,renderer = magick_renderer(), height = size[1])
  close.gif <- animate(p.anim, renderer = magick_renderer(), height = size[2], width =size[1])

  dims <- magick::image_info(close.gif[1]) %>% select(width,height) %>% unlist

  tdir <- tempdir()

  ggsave(
    filename = file.path(tdir, "out_001.png"),
    plot = ggdraw(xlim = c(0,1000),ylim=c(0,1000)),  device = "png",width = dims[1],height = dims[2],units = "px",dpi = "retina",bg = "white")


  for(i in 1:length( perf.gif)){
    new_gif <- plot_grid(ggdraw() + draw_image(close.gif[i], scale = 1.1,halign=.3),
                         ggdraw() +draw_image(perf.gif[i], scale = 0.3,valign = 0.7),
                         ncol=2)
    ggsave(
      filename = file.path(tdir, paste0("out_", sprintf("%03d", i+1), ".png")),
      plot = new_gif,  device = "png",width = dims[1],height = dims[2],units = "px",dpi = "retina")

  }

  gif_file <- file.path(tempdir(), 'final.gif')

  png_files <- sort(list.files(path = tdir, pattern = "out_", full.names = TRUE))
  gifski(png_files, gif_file = gif_file, width = dims[1], height = dims[2], delay = .1,
         progress = FALSE)

  viewer <- getOption("viewer")
  if (!is.null(viewer))
    viewer(gif_file)
  else
    utils::browseURL(gif_file)
  }

  if(render&!inset){
    FPS <- fps
   gganimate::animate(p.anim+sh, fps=FPS,start_pause=2,end_pause=2, height = size[2], width =size[1])
  }

  if(!render&!inset) return(p.anim)
}


if(F){
ggplot()+geom_path(data=jaw.df %>% filter(Time==0.1),aes(x,y))+geom_path(data=head.df %>% filter(Time==0.1),aes(x,y))+geom_point(data=x,aes(x,y))+geom_line(data=x,aes(x,y,group=path,color=type))+theme_void()+transition_time(as.integer(Time))+labs(title = "Time: {frame_time} ms")+xlim(x.lim)+ylim(y.lim)+ coord_fixed()
}

#' @title  Plot and print descriptions of configurations supported by \code{JawCloseR}

#' @description  Plots approximate positions of muscle divisions in the different configurations supported by JawCloseR.  Prints a description of each muscle division to console.
#'
#' @param config numeric, the configuration.

#' @details This function is intended to help the user choose which configuration to implement in \code{start_sim()}. The choice will, of course, depend on the species being modeled. For the most part, the names of muscle divisions are arbitrary with regards to the model. Much ink has been spilled over the nomenclature and homology of fish jaw adductors. The \eqn{A_2} and \eqn{A_3} divisions can represent separate adductors applying torque to the lower jaw, however, at the same inlever position in configurations 1 and 5. The implementation of the \eqn{A_1} and \eqn{A_\omega} divisions (configurations 4 and 5), however, represent special cases. \eqn{A_1} is implemented as a division applying torque to the lower jaw via a tendon with inserting on a stationary dorsal point (i.e., the origin of the \eqn{A_1} tendon on the maxilla). \eqn{A_\omega} represents an intramandibular division, originating along the lower jaw and either inserting on the Meckelian tendon (along with \eqn{A_2} and \eqn{A_3} in configuration 3) or on suspensorium (configuration 5).
#'
#'@return A text description of the configurations is printed to the console along with a \code{ggplot} object.
#'
#'
#'
#'
#' @export
#' @import ggplot2
#'
#' @examples
#'
#'print_configs(4)
#'
print_configs <- function(config=NULL){


  jaw <- system.file("extdata/jaws","LMBass.txt", package = "JawCloseR") %>% read_tsv(col_types = cols()) %>% rownames_to_column(var="n") %>% mutate(n=as.numeric(n)) %>% mutate(part=ifelse(n==last(n),"JawTip","Jaw")) %>% mutate(part=ifelse(n==(last(n)-1),"joint.p",part)) %>% filter(part=="Jaw")

  upperjaw <- system.file("extdata/jaws","LMBass_max.txt", package = "JawCloseR") %>% read_tsv(col_types = cols())

  head <- system.file("extdata/jaws", "LMBass_head.txt", package = "JawCloseR") %>% read_tsv(col_types = cols()) %>% rownames_to_column(var="n") %>% mutate(n=as.numeric(n)) %>% mutate(part=ifelse(n==last(n),"HeadTip","Head")) %>% mutate(part=ifelse(n==(last(n)-1),"joint.p",part))%>% filter(part=="Head")

  cons <- system.file("extdata", "configurations.csv", package = "JawCloseR") %>% read_csv(col_types = cols())


  config1 <- data.frame(
    a=c("A1Or","AntSeg","PostSeg","A2Or","A3Or","AwIns","joint.p"),
    b=c("A1Ins","A1Ins","A1Ins", "A2Ins","A3Ins","AwOr","JawTip"),
    path=c("A1","A1AntSeg","A1PostSeg","A2","A3","Aw","jaw"),
    type=c("muscle","A1 tend","A1 tend","muscle","muscle",'muscle',"jaw"),
    col=c("red","gray","gray","red","red","red","black")
  ) %>% pivot_longer(cols = a:b,values_to="name",names_to = "position") %>%
    mutate(x=c(450,550,605,550,550,550,455,560,460,560,560,605,530,655),
           y=c(150,120,175,120, 60,120,100, 85, 105, 85, 85, 65, 68, 68),
           config=1) %>%
    filter(!grepl("A1",path) & !grepl("Aw",path)) %>%
    mutate(path=ifelse(path=="A2","A2 \n (A2 mass +Aw mass)",path))



  config2 <- data.frame(
    a=c("A1Or","AntSeg","PostSeg","A2Or","A3Or","AwIns","joint.p"),
    b=c("A1Ins","A1Ins","A1Ins", "A2Ins","A3Ins","AwOr","JawTip"),
    path=c("A1","A1AntSeg","A1PostSeg","A2","A3","Aw","jaw"),
    type=c("muscle","A1 tend","A1 tend","muscle","muscle",'muscle',"jaw"),
    col=c("red","gray","gray","red","red","red","black")
  ) %>% pivot_longer(cols = a:b,values_to="name",names_to = "position") %>%
    mutate(x=c(450,550,605,550,550,550,455,560,460,560,560,605,530,655),
           y=c(150,120,175,120, 60,120,100, 85, 105, 85, 85, 65, 68, 68),
           config=2) %>%
    filter(!grepl("A1",path) & !grepl("Aw",path))



  config3 <- data.frame(
    a=c("A1Or","AntSeg","PostSeg","A2Or","A3Or","AwIns","joint.p"),
    b=c("A1Ins","A1Ins","A1Ins", "A2Ins","A3Ins","AwOr","JawTip"),
    path=c("A1","A1AntSeg","A1PostSeg","A2","A3","Aw","jaw"),
    type=c("muscle","A1 tend","A1 tend","muscle","muscle",'muscle',"jaw"),
    col=c("red","gray","gray","red","red","red","black")
  ) %>% pivot_longer(cols = a:b,values_to="name",names_to = "position") %>%
    mutate(x=c(450,550,605,550,550,550,455,560,460,560,560,605,530,655),
           y=c(150,120,175,120, 60,120,100, 85, 105, 85, 85, 65, 68, 68),
           config=3) %>%
    filter(!grepl("A1",path))


  #largemouth
  config4 <- data.frame(
    a=c("A1Or","AntSeg","PostSeg","A2Or","A3Or","AwIns","joint.p"),
    b=c("A1Ins","A1Ins","A1Ins", "A2Ins","A3Ins","AwOr","JawTip"),
    path=c("A1","A1AntSeg","A1PostSeg","A2","A3","Aw","jaw"),
    type=c("muscle","A1 tend","A1 tend","muscle","muscle",'muscle',"jaw"),
    col=c("red","gray","gray","red","red","red","black")
  ) %>% pivot_longer(cols = a:b,values_to="name",names_to = "position") %>%
    mutate(x=c(450,550,585,550,550,550,455,550,460,560,510,580,530,655),
           y=c(150,105,140,105, 60,105,100,105,105, 85, 80, 65, 68, 68),
           config=4) %>%
    filter(path!="A2") %>%
    mutate(path=ifelse(path=="A1","A1/A2",path))



config5 <- data.frame(
    a=c("A1Or","AntSeg","PostSeg","A2Or","A3Or","AwIns","joint.p"),
    b=c("A1Ins","A1Ins","A1Ins", "A2Ins","A3Ins","AwOr","JawTip"),
    path=c("A1","A1AntSeg","A1PostSeg","A2","A3","Aw","jaw"),
    type=c("muscle","A1 tend","A1 tend","muscle","muscle",'muscle',"jaw"),
    col=c("red","gray","gray","red","red","red","black")
  ) %>% pivot_longer(cols = a:b,values_to="name",names_to = "position") %>%
    mutate(x=c(450,550,585,550,550,550,455,560,460,560,510,605,530,655),
           y=c(150,105,140,105, 60,105,100, 85, 105, 85, 80, 65, 68, 68),
           config=5)


configs <- list(config1,config2,config3,config4,config5)
config.df <- configs[[config]]

muscles <- config.df %>%
  filter(type=="muscle") %>%
  group_by(path) %>%
  mutate(OrIn=ifelse(grepl("Or",name),"or","ins")) %>%
  select(-position,-name) %>%
  pivot_wider(values_from=c("x","y"),names_from=OrIn)


muscle_tri<- function(x_or,y_or,x_ins,y_ins,config=NULL,musc=NULL){
  fact <- 2
  theta=atan((y_or-y_ins)/(x_or-x_ins))

  ad <- pi
  if(!is.null(config) & config==3 & musc=="Aw"){ ad=-pi/10;fact=1.25}
  if(!is.null(config) & config %in% c(4,5) & musc=="Aw"){ ad=-pi/30;fact=1.25}
  pt1=point_ang_r(c(x_ins,y_ins),dist_2d(c(x_or,y_or),c(x_ins,y_ins)),ad+theta*fact)#upper
  pt2=point_ang_r(c(x_ins,y_ins),dist_2d(c(x_or,y_or),c(x_ins,y_ins)), ad-theta/fact)#lower
  xy <- data.frame(x=c(pt1[1],pt2[1],x_or,x_ins),y=c(pt1[2],pt2[2],y_or,y_ins))
  #xy <- data.frame(x=c(pt1[1],pt2[1]),y=c(pt1[2],pt2[2]))
  return(xy[chull(xy),])
}

tris <- muscles %>%
  group_by(path) %>%
  summarize(muscle_tri(x_or,y_or,x_ins,y_ins,config = config,musc=path))

labs <- tris %>%
  group_by(path) %>%
  summarize_at(c("x","y"),mean)

labs2 <- config.df %>%
  filter(type!="muscle") %>%
  group_by(type) %>%
  summarize_at(c("x","y"),mean)

p <- ggplot()+geom_path(data=jaw,aes(x,y),col="gray60")+geom_path(data=head,aes(x,y),col="gray60")+geom_path(data=upperjaw,aes(x,y),col="gray60",alpha=0.6)+geom_line(data=config.df %>% filter(type!="muscle"),aes(x,y,group=path),col="gray10")+coord_fixed()+geom_point(data=config.df,aes(x,y),alpha=0.5)+geom_polygon(data=tris,aes(x,y,col=path),alpha=0.1)+geom_text(data=labs,aes(x,y,label=path,col=path))+geom_text(data=labs2,aes(x,y,label=type),nudge_x = 25,nudge_y = 10)+ggtitle(label=paste0("configuration ",config.df$config[1]))+theme_void()+theme(legend.position = "none")

print(p)
con2 <- cons %>%
  filter(configuration==config)

for(i in con2$division) {
  or.i <-con2 %>% filter(division==i) %>% select(-configuration,-division)
  row.i <- paste0(names(or.i),":", c("\t\t","\t","\t","\t\t","\t","\t"))
  cat(paste("________________",i,"_______________"))
  cat("\n")
  cat(paste(row.i,or.i) %>% unlist, sep = "\n")
  cat("\n")
  cat("\n")
}

}


#' @title  Plot and print descriptions of configurations supported by \code{JawCloseR}

#' @description  Plots approximate positions of muscle divisions in the different configurations supported by JawCloseR.  Prints a description of each muscle division to console.
#'
#' @param config \code{data.frame} of x and y coordinates in the first and second columns, respectively, or a 1D vector of x and y in that order.
#' @param rot numeric, the rotation angle in radians
#' @param center numeric, a vector of x and y, the orig about which to rotate the points
#'
#'@return a \code{data.frame} or vector, depending on what is passed to \code{xy}
#'
#'
#'
#'
#' @export



rotate_coords <- function(xy, rot, center) {

  if(!is.data.frame(xy) & length(xy)==2) {xy2 <- data.frame(x=xy[1],y=xy[2])}else{xy2 <- xy}
  co <- cos(rot)
  si <- sin(rot)
  adj <- matrix(rep(center, nrow(xy2)), ncol = 2, byrow = TRUE)
  xy2 <- xy2 - adj
  res <- cbind(co * xy2[, 1] - si * xy2[, 2], si * xy2[, 1] + co *
               xy2[, 2]) + adj
  res <- res %>% data.frame()
  colnames(res) <- c("x","y")
  if(!is.data.frame(xy) & length(xy)==2) res <- c(res[,1],res[,2])
  return(res)
}

if(F){

perf_inset <- function(x){

  x=perf
  vals <- x %>%
  pivot_longer(c(starts_with("A"),starts_with("Static")))


max.val <- x %>%
  na.omit %>%
  summarise_all(range)

max.val[1,] <- 0

max <- max.val %>%
  mutate(m=c("min","max")) %>%
  select(-Time,-JawAng) %>%
  pivot_longer(c(starts_with("A"),starts_with("Static"))) %>%
  pivot_wider(names_from=m)

vals <- vals %>%
  filter(Time>0) %>%
  left_join(max) %>%
  mutate(per=value/max)


vals %>%  ggplot()+
  geom_bar(aes(name,per,fill=per), stat="identity")+
  scale_fill_viridis(option = "H") +
  facet_wrap(~name,scales = "free_x",nrow=2)+
  theme_void()+transition_time(Time)+labs(title = "Time: {frame_time %>% round(.,1)} ms")+ylim(c(0,1))

tdir <- tempdir()
for(i in 1:100){
  new_gif <- plot_grid(ggdraw() + draw_image(a_gif[i], scale = 0.9),
                       ggdraw() + draw_image(d_gif[i], scale = 0.9),
                       ggdraw() + draw_image(b_gif[i], scale = 0.9),
                       ggdraw(),
                       ggdraw() + draw_image(c_gif[i], scale = 0.9),
                       ncol=2)
  ggsave(
    filename = file.path(tdir, paste0("out_", sprintf("%03d", i), ".png")),
    plot = new_gif, width = 2.4, height = 3.6, device = "png")
}

png_files <- sort(list.files(path = tdir, pattern = "out_", full.names = TRUE))
gifski(png_files, gif_file = "out.gif", width = 480, height = 720, delay = .1,
       progress = TRUE)
angvel.exp <- expression(radians ~s^-1)
vel.exp <- expression(ML ~s^-1)


stack <- expand_grid(name=vals$name %>% unique,per=seq(0,1,0.1))

vals$name <- factor(vals$name,levels=c("A2v","A3v","Awv","AngVel","A2Power","A3Power","AwPower","StaticBite"))

stack$name<- factor(stack$name,levels=c("A2v","A3v","Awv","AngVel","A2Power","A3Power","AwPower","StaticBite"))

exps <- tibble(
  exp=c(rep(as.character(vel.exp),3),as.expression(angvel.exp),rep("W",3),"N"),
  name=levels(vals$name)
)


vals2 <- vals %>%
  group_by(Time,name) %>%
  na.omit %>%
  filter(Time>0) %>%
  mutate(per=ifelse(per<0,0,per)) %>%
  summarize(per=seq(0,round(per,1),0.1))


p <- ggplot()+
  geom_bar(data=stack,aes(name,per),position="stack", stat="identity",col="gray90",fill="white")+
  geom_bar(data=vals2,aes(name,per,fill=per),position="stack", stat="identity")+
  scale_fill_viridis(option = "H") +
  facet_wrap(~name,scales = "free",nrow=2)+
  theme_void()+
  theme(strip.background = element_rect(
   fill="gray"),strip.text = element_text(margin = margin(b = 5, t = 5)),legend.position = "none")+
  geom_text(data=vals,aes(x=name,y=6,label=signif(value,2)))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+ylab(angvel.exp)+transition_time(as.integer(Time))+labs(title = "Time: {frame_time} ms")

gganimate::animate(p, fps=10)
return(p)
}

perf_inset(perf)

}
