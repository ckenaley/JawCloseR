
i <- n


j[n,"Time"]<- j[n,"Time"]+TimeStep
j[n,"Iteration"]<- i
		#************* pressure **********
		j[n,"PressureT"]=-0.40373*(pi/2)*(a/Scale)^2*((JawWidth/2)/Scale)*((Gape-j[n-1,"JawAng"])/Gape)*Pressure
		
		#************* Acceleration and velocity **********
	
j[n,"AngAcc"]=(j[n-1,"AwTorque"]+j[n-1,"TendonTorque"]+j[n-1,"DragT"]+j[n-1,"PressureT"])/Itot

		j[n,"AngVel"]=j[n-1,"AngAcc"]*TimeStepSecs+j[n-1,"AngVel"];
		j[n,"dJawAng"]=j[n,"AngVel"]*TimeStepSecs; 
		j[n,"JawAng"]=j[n-1,"JawAng"]-j[n,"dJawAng"]
	
	open.angle <- -1*(start+j[n,"JawAng"]-SpecAng)

	#************* Drag torques ************

	j[n,"DragT"]=(-2/15)*1000*j[n-1,"AngVel"]^2*(a/Scale)^4*(JawWidth/2)/Scale	


inlever.p <- point.ang.r(joint.p,A2Li,open.angle)
JawTip<- point.ang.r(joint.p,JL,open.angle)

j[n,"Time"]=j[n-1,"Time"]+TimeStep;
	

#remove negative muscle forces
if(j[n-1,"A2Force"]<0) j[n-1,"A2Force"] <- 1e-15
if(j[n-1,"A3Force"]<0) j[n-1,"A3Force"] <- 1e-15
if(j[n-1,"AwForce"]<0) j[n-1,"AwForce"] <- 1e-15


###### cart coordinates at nth iteration (instantaneous)####

AwOr<- point.ang.r(joint.p,AwLi,open.angle)#this is dynamic
#which line of action inputs most torque 



####### RESULTANTS  #############

#****** find resultant of A2 and A3 forces
#need ang between A2 and A3 LOA 

A2A3.ang <- cos.ang(j[n-1,"A3ML"],j[n-1,"A2ML"],dist.2d(A2Or,A3Or))

#resultant and its angle
res <- sqrt(j[n-1,"A2Force"]^2+j[n-1,"A3Force"]^2+2*j[n-1,"A2Force"]*j[n-1,"A3Force"]*cos(A2A3.ang))


res.theta <- atan((j[n-1,"A2Force"]*sin(A2A3.ang))/(j[n-1,"A3Force"]+j[n-1,"A2Force"]*cos(A2A3.ang)))
if(res==0) res.theta=A2A3.ang

#res <- sqrt(j[n,"A2Force"]^2+j[n,"A3Force"]^2+2*j[n,"A2Force"]*j[n,"A3Force"]*cos(A2A3.ang))

#angle between line def by A2 and A3 origs and A3ML
A2A3.A3ang <- cos.ang(dist.2d(A2Or,A3Or),j[n-1,"A3ML"],j[n-1,"A2ML"]) 

#theoretical length of resultant forces from A2 and A3
res.ML <- sin(A2A3.A3ang)*(j[n-1,"A3ML"]/sin(pi-A2A3.A3ang-res.theta))

#using resultant lengths, compute resultants position
#distance btween A3 orig and resultant orig
A3.res <- cos.side(j[n-1,"A3ML"],res.ML,res.theta)
#cos.side <- function(l,r,ang){sqrt(l^2+ r^2-2*l*r* cos(ang))}
resOr <- point.ang.r(A3Or,
                     A3.res,
                     pi-cos.ang(dist.2d(A3Or,c(A2Or[1],A3Or[2])),dist.2d(A3Or,A2Or),dist.2d(A2Or,c(A2Or[1],A3Or[2]))) #horiz ang of A3or-A2or
)



tendon.theta <- cos.ang(A2Li,dist.2d(inlever.p,resOr),dist.2d(resOr,joint.p)) #tendon ang based on resultants

tendon.thetaNoAw <- tendon.theta

nexus <- point.ang.r(inlever.p,j[1,"TendonLength"],pi+open.angle-tendon.theta)



######muscle lengths
m.max <- ifelse(i==2,1,which.max(abs(c(j[n-1,"A2TendonTorque"],j[n-1,"A3TendonTorque"],j[n-1,"AwTendonTorque"]))))


### reposition according to sag and equilibrium imposed by Aw on resultant line.
nexus <- point.ang.r(inlever.p,j[1,"TendonLength"],pi+open.angle-tendon.theta)

j[n,"AwML"] <- dist.2d(nexus,AwOr)
j[n,"A2ML"] <-dist.2d(nexus,A2Or)
j[n,"A3ML"] <-dist.2d(nexus,A3Or)

### Add effect of Aw on tendon position, that is, sag imposed by its force
#if(j[n,"AwML"]>=min.ML*AwMLClosed & config==3){
AwTendonTheta <- cos.ang(j[n,"AwML"],j[1,"TendonLength"],dist.2d(AwOr,inlever.p))

AwTendonF <- sin(AwTendonTheta)*j[n-1,"AwForce"]

resTendonTheta <- pi-asin(AwTendonF/(AwTendonF+res))

res.theta.sag <- asin(((resTendonTheta)*j[1,"TendonLength"])/dist.2d(inlever.p,resOr))

TendonThetaAdj <- pi-res.theta.sag-resTendonTheta


#### recalculate tendon theta and nexus according to sag (but only for config==3, sling in place) 
 if(config==3) tendon.theta <- tendon.theta-TendonThetaAdj #tendon ang based on resultants


if(n>2 & tendon.theta>(pi-rad(10))) tendon.theta <- pi-rad(10) #include in MS, keeps tendon from being flat

nexus <- point.ang.r(inlever.p,j[1,"TendonLength"],pi+(open.angle-tendon.theta))

#Meckelian tendon theta ignoring Aw
noAw.TendonTheta <- ifelse(m.max==1,cos.ang(A2Li,dist.2d(inlever.p,A2Or),dist.2d(A2Or,joint.p)),cos.ang(A2Li,dist.2d(inlever.p,A3Or),dist.2d(A3Or,joint.p))) 
noAw.nexus <- point.ang.r(inlever.p,j[1,"TendonLength"], pi+open.angle-noAw.TendonTheta)

if(n==2) {tendon.theta <- noAw.TendonTheta 
nexus <- noAw.nexus
}

j[n,"AwML"] <- dist.2d(nexus,AwOr)
j[n,"A2ML"] <-dist.2d(nexus,A2Or)
j[n,"A3ML"] <-dist.2d(nexus,A3Or)

#}



# #Meckelian tendon theta ignoring Aw (for 1st it, too)
# noAw.TendonTheta <- ifelse(m.max==1,cos.ang(A2Li,dist.2d(inlever.p,A2Or),dist.2d(A2Or,joint.p)),tendon.theta) 
# 
# noAw.nexus <- point.ang.r(inlever.p,j[1,"TendonLength"], pi+open.angle-noAw.TendonTheta)
# 
# if(config==1|config==2|n==2){
# tendon.theta <- noAw.TendonTheta 
# nexus <- noAw.nexus
# }



#could constraint AwML to x% of AwFLopt (0.6 default, reasonable according to 
#T.J. Burkholder, R.L. Lieber. Journal of Experimental Biology 2001 204: 1529-1536)

#if(j[n,"AwML"]<min.ML*AwMLClosed){
  #  tendon.theta <- cos.ang(j[1,"TendonLength"],dist.2d(inlever.p,AwOr),min.ML*AwMLClosed)
   # nexus <- point.ang.r(inlever.p,j[1,"TendonLength"], pi+open.angle-tendon.theta)
    #j[n,"AwML"] <- min.ML*AwMLClosed
#}


j[n,"A2ML"] <-dist.2d(nexus,A2Or)
j[n,"A3ML"] <-dist.2d(nexus,A3Or)


j[n,"AwTheta"] <- cos.ang(dist.2d(inlever.p,AwOr),j[n,"AwML"],j[1,"TendonLength"])

if(n>22){range.n <- (n-21):(n-1)}else{range.n <- n-1}

# nexus is bound by lines of action of A2 and A3
j[n,"TendonLength"] <-j[1,"TendonLength"] 

j[n,"TendonTheta"]=cos.ang(A2Li,j[1,"TendonLength"],dist.2d(joint.p,nexus))

		#************* A3 force ********
		j[n,"A3dML"] =j[n-1,"A3ML"]-j[n,"A3ML"];TimeStepSecs
		j[n,"A3v"] =abs(j[n,"A3dML"]/A3MLOpen/TimeStepSecs);
		j[n,"A3FL"] =sqrt((sin(j[n-1,"A3Penn"])*j[n-1,"A3FL"])^2+((cos(j[n-1,"A3Penn"])*j[n-1,"A3FL"])-j[n,"A3dML"])^2);
		j[n,"A3Penn"] =asin((sin(j[n-1,"A3Penn"])*j[n-1,"A3FL"])/j[n,"A3FL"]);
		j[n,"A3Fmax"] =A3FmaxNoPenn*cos(j[n,"A3Penn"]);
		j[n,"A3Ffv"] =vl.ffv(V=mean(j[range.n,"A3v"]))
		#j[n,"A3Ffv"] =(Vmax-j[n-1,"A3v"])/(Vmax+j[n-1,"A3v"]*G);
		ifelse(j[n,"A3FL"] >A3FLOpt,j[n,"A3Ffl"]<- -6.25*(j[n,"A3FL"]/A3FLOpt)^2+12.5*(j[n,"A3FL"]/A3FLOpt)-5.25,j[n,"A3Ffl"]  <- 1)
		
		#ifelse(j[n,"A3FL"]/A3FLOpt >= 0.5 & j[n,"A3FL"]/A3FLOpt <= 1.56, j[n,"A3Ffl"]<- predict.lm(fl.poly,data.frame(p.fl=j[n,"A3FL"]/A3FLOpt)),j[n,"A3Ffl"]  <- 0)# from Porro et al. (2002)
		
		if(j[n,"Time"]<ActRiseTime){j[n,"A3Fact"]=0.5-0.5*cos(pi*j[n,"Time"]/ActRiseTime)} else{j[n,"A3Fact"]=1};
		
			ifelse(j[n,"A3ML"]>A3MLOpt,j[n,"A3Fpar"] <- A3FmaxClosed*exp(2*log(1+MaxFpar)*(j[n,"A3ML"]/A3MLOpt-1))-A3FmaxClosed,j[n,"A3Fpar"] <- 0)
		j[n,"A3Force"] =j[n,"A3Fmax"] *j[n,"A3Ffv"]*ifelse(j[n,"A3Ffl"]<0,0,j[n,"A3Ffl"])*j[n,"A3Fact"]+j[n,"A3Fpar"];


		#************* A2 force ********
		j[n,"A2dML"] =j[n-1,"A2ML"]-j[n,"A2ML"];
		j[n,"A2v"] =abs(j[n,"A2dML"]/A2MLOpen/TimeStepSecs);
		j[n,"A2FL"] =sqrt((sin(j[n-1,"A2Penn"])*j[n-1,"A2FL"])^2+((cos(j[n-1,"A2Penn"])*j[n-1,"A2FL"])-j[n,"A2dML"])^2);
		j[n,"A2Penn"] =asin((sin(j[n-1,"A2Penn"])*j[n-1,"A2FL"])/j[n,"A2FL"]);
		j[n,"A2Fmax"] =A2FmaxNoPenn*cos(j[n,"A2Penn"]);
		j[n,"A2Ffv"] =vl.ffv(V=mean(j[range.n,"A2v"]));
		#j[n,"A2Ffv"]=(Vmax-j[n-1,"A2v"])/(Vmax+j[n-1,"A2v"]*G);
		ifelse(j[n,"A2FL"] >A2FLOpt,j[n,"A2Ffl"]<- -6.25*(j[n,"A2FL"]/A2FLOpt)^2+12.5*(j[n,"A2FL"]/A2FLOpt)-5.25,j[n,"A2Ffl"]  <- 1)
		
		#ifelse(j[n,"A2FL"]/A2FLOpt >= 0.5 & j[n,"A2FL"]/A2FLOpt <= 1.56, j[n,"A2Ffl"]<- predict.lm(fl.poly,data.frame(p.fl=j[n,"A2FL"]/A2FLOpt)),j[n,"A2Ffl"]  <- 0)# from Porro et al. (2002)
		
		if(j[n,"Time"]<ActRiseTime){j[n,"A2Fact"]=0.5-0.5*cos(pi*j[n,"Time"]/ActRiseTime)} else{j[n,"A2Fact"]=1};
		ifelse(j[n,"A2ML"]>A2MLOpt,j[n,"A2Fpar"] <- A2FmaxClosed*exp(2*log(1+MaxFpar)*(j[n,"A2ML"]/A2MLOpt-1))-A2FmaxClosed,j[n,"A2Fpar"] <- 0)
			
		j[n,"A2Force"] =j[n,"A2Fmax"] *j[n,"A2Ffv"]*ifelse(j[n,"A2Ffl"]<0,0,j[n,"A2Ffl"])*j[n,"A2Fact"]+j[n,"A2Fpar"];


		
				#************* Aw force ********
		j[n,"AwdML"] =j[n-1,"AwML"]-j[n,"AwML"];
		j[n,"Awv"] =abs(j[n,"AwdML"]/AwMLOpen/TimeStepSecs);
		j[n,"AwFL"] =sqrt((sin(j[n-1,"AwPenn"])*j[n-1,"AwFL"])^2+((cos(j[n-1,"AwPenn"])*j[n-1,"AwFL"])-j[n,"AwdML"])^2);
		j[n,"AwPenn"] =asin((sin(j[n-1,"AwPenn"])*j[n-1,"AwFL"])/j[n,"AwFL"]);
		j[n,"AwFmax"] =AwFmaxNoPenn*cos(j[n,"AwPenn"]);
		j[n,"AwFfv"] =vl.ffv(V=mean(j[range.n,"Awv"]))
		#j[n,"AwFfv"] =abs(Vmax-j[n-1,"Awv"])/(Vmax+j[n-1,"Awv"]*G);
		ifelse(j[n,"AwFL"] >AwFLOpt,j[n,"AwFfl"]<- -6.25*(j[n,"AwFL"]/AwFLOpt)^2+12.5*(j[n,"AwFL"]/AwFLOpt)-5.25,j[n,"AwFfl"]  <- 1) #old length-tension relationship from Van Wassenbergh
		
		
		#ifelse(j[n,"AwFL"]/AwFLOpt >= 0.5 & j[n,"AwFL"]/AwFLOpt <= 1.56, j[n,"AwFfl"]<- predict.lm(fl.poly,data.frame(p.fl=j[n,"AwFL"]/AwFLOpt)),j[n,"AwFfl"]  <- 0)# from Porro et al. (2002)
		
	
		if(j[n,"Time"]<ActRiseTime){j[n,"AwFact"]=0.5-0.5*cos(pi*j[n,"Time"]/ActRiseTime)} else{j[n,"AwFact"]=1};
		ifelse(j[n,"AwML"]>AwMLOpt,j[n,"AwFpar"] <- AwFmaxClosed*exp(2*log(1+MaxFpar)*(j[n,"AwML"]/AwMLOpt-1))-AwFmaxClosed,j[n,"AwFpar"] <- 0)
		
		j[n,"AwForce"] =j[n,"AwFmax"] *j[n,"AwFfv"]*ifelse(j[n,"AwFfl"]<0,0,j[n,"AwFfl"])*j[n,"AwFact"]+j[n,"AwFpar"];


	#************** Torques 	#Scaled to N*m  *****************
	#Angles of A2 and A3 with Aw, both constrained to conform to A2-AW line by F condition, see above for d

	j[n,"A3AwTheta"] <-  cos.ang(j[n,"AwML"],j[n,"A3ML"],dist.2d(A3Or,AwOr))#ifelse(line.max==2,rad(180),cos.ang(j[n,"A3ML"],j[n,"AwML"],dist.2d(A3Or,AwOr)))
	
	j[n,"A2AwTheta"] <-cos.ang(j[n,"A2ML"],j[n,"AwML"],dist.2d(A2Or,AwOr))
	

#Angles of A2 and A3 with tendon
	j[n,"A3TendonTheta"] <-cos.ang(j[n,"TendonLength"],j[n,"A3ML"],dist.2d(inlever.p,A3Or))

	j[n,"A2TendonTheta"] <-cos.ang(j[n,"TendonLength"],j[n,"A2ML"],dist.2d(inlever.p,A2Or))

	j[n,"AwTendonTheta"] <- cos.ang(j[n,"TendonLength"],j[n,"AwML"],dist.2d(AwOr,inlever.p))


	#************** Aw 
#j[n,"AwTorque"]=(j[n,"AwForce"]+j[n,"A2Force"]*cos(j[n,"A2AwTheta"])*-1+j[n,"A3Force"]*-1*cos(j[n,"A3AwTheta"])*sin(j[n,"AwTheta"]))*AwLi/Scale
	
	#remove negative muscle forces
	if(j[n,"A2Force"]<0) j[n,"A2Force"] <- 1e-15
	if(j[n,"A3Force"]<0) j[n,"A3Force"] <- 1e-15
	if(j[n,"AwForce"]<0) j[n,"AwForce"] <- 1e-15
	
j[n,"AwTorque"]=min(c(
  j[n,"AwForce"]*abs(cos(j[n,"AwTheta"])), #Aw component
  j[n,"A2Force"]*abs(cos(j[n,"A2AwTheta"]))+j[n,"A3Force"]*abs(cos(j[n,"A3AwTheta"])) #facialis components
))*AwLi/Scale


#j[n,"AwTorque"] <- ifelse(line.max==4,0,j[n,"AwTorque"])

#(j[n,"AwForce"]+j[n,"A2Force"])*sin(j[n,"AwTheta"])*AwLi/Scale+(j[n,"A3Force"])*sin(cos.ang(AwLi,dist.2d(AwOr,A3Or),dist.2d(joint.p,A3Or)))*AwLi/Scale


	#************** Tendon input

j[n,"TendonTorque"]=((j[n,"A2Force"]*abs(cos(j[n,"A2TendonTheta"])+j[n,"A3Force"]*abs(cos(j[n,"A3TendonTheta"])))))*sin(j[n,"TendonTheta"])*A2Li/Scale

					#************** Torques about tendon*****************
	#************** A2
	#Scaled to N*m 
		j[n,"A3TendonTorque"]=j[n,"A3Force"]*j[1,"TendonLength"]/Scale*sin(j[n,"A3TendonTheta"]);#A3TendonTheta is outside?
	#************** A2
		j[n,"A2TendonTorque"]=j[n,"A2Force"]*j[1,"TendonLength"]/Scale*sin(j[n,"A2TendonTheta"]);#A2TendonTheta is outside?
		#************** Aw	
		j[n,"AwTendonTorque"]=j[n,"AwForce"]*j[1,"TendonLength"]/Scale*sin(j[n,"AwTendonTheta"]);#AwTendonTheta is outside?
			#****** Tendon: Torque applied by A2, A3, Aw at tendon insertion		

#************* performance
j[n,"AwFout"]= j[n,"AwTorque"]/(AwLo/Scale) #Aw bite force at tip of the jaw
j[n,"TendFout"]=j[n,"TendonTorque"]/(A2Lo/Scale) #A2 bite force at tip of the jaw
		
j[n,"JawVelTip"]=j[n,"AngVel"]*(A2Lo) #in cm/s
j[n,"StaticBite"]=j[n,"AwFout"]+j[n,"TendFout"] #Total bite force at tip of jaw
j[n,"EMA"]=j[n,"StaticBite"]/(j[n,"AwForce"]+j[n,"A2Force"]+j[n,"A3Force"]) #Total bite force at tip of jaw



######### p.var it coords

inset.dat <- c(n=j[n,"Time"],
               JawAng=deg(j[n,"JawAng"]),
               res.ML=res.ML,
               #A3TendTor=j[n,"A3TendonTorque"],
               #AwTendTor=j[n,"AwTendonTorque"],
               #A2TendTheta=deg(j[n,"A2TendonTheta"]),
               #A3TendTheta=deg(j[n,"A3TendonTheta"]),
               #AwTendTheta=deg(j[n,"AwTendonTheta"]),
               #A3AwTheta=deg(j[n,"A3AwTheta"]),
               A2AwTheta=deg(j[n,"A2AwTheta"]),
               AwTor=j[n,"AwTorque"],
               TendTor=j[n,"TendonTorque"],
               A2F=j[n,"A2Force"],A3F=j[n,"A3Force"],
               AwF=j[n,"AwForce"],AwFfv=j[n,"AwFfv"],
               A2Ffv=j[n,"A2Ffv"],A3Ffv=j[n,"A2Ffv"],
               A3Ffl=j[n,"AwFfl"],Awv=j[n,"Awv"],
               AwdML=j[n,"AwdML"],A2v=j[n,"A2v"],
               A2v=j[n,"A2v"],acc=j[n,"AngAcc"],
               drag=j[n,"DragT"],
               pres=j[n,"PressT"],
               inertia=Itot)
inset.names <- names(inset.dat)
n.dat <- data.frame(dat=inset.dat,name=inset.names,x=rep(-2,length(inset.dat)),y=1:length(inset.dat)*-0.25+2.5)

#get signs right
plot.ls <- list(joint.p=joint.p,A2Or=A2Or,A3Or=A3Or,JawTip=JawTip,inlever=inlever.p,AwIns=AwOr,nexus=nexus,resOr=resOr)
geom.df <- as.data.frame(do.call(rbind,plot.ls))

geom.out[[paste(n)]] <- data.frame(ms=j[n,"Time"],ang=j[n,"JawAng"],var=unlist(plot.ls))

geom.df$var <- rownames(geom.df)
colnames(geom.df) <- c("x","y","var")
geom.out[[paste(n)]] <- data.frame(ms=j[n,"Time"],ang=j[n,"JawAng"],geom.df)



if(progress==T){
  if(i %in% c(2,seq(print.every,MaxIts,print.every))){
    
    
    circs <- data.frame(t(sapply(seq(pi/10,pi*2,by = pi/10),function(x) point.ang.r(inlever.p,j[1,"TendonLength"],x))))
    colnames(circs) <- c("x","y")
    #p <- ggplot(data=geom.df,aes(x=x,y=y))+geom_point(aes(colour=var),size=3)+theme_classic(20)+ylim(-4,4)+line.2pt("joint.p","JawTip",geom.df)+theme(legend.position="none")+geom_text(data=n.dat,aes(x=x-1.25,y=y,label=paste0(name,": ",dat)),size=2,hjust=0)+line.2pt("nexus","AwIns",geom.df,col="darkgreen")+line.2pt("A3Or","nexus",geom.df,col="brown")+line.2pt("A2Or","nexus",geom.df,col="red")+line.2pt("inlever.p","nexus",geom.df,col="gray")+xlim(-4,4)+geom_point(aes(x=x,y=y),dat=circs,size=0.5)
    suppressWarnings( p <- ggplot(data=geom.df,aes(x=x,y=y))+geom_point(aes(colour=var),size=3)+theme_classic(20)+ylim(-4,4)+line.2pt("joint.p","JawTip",geom.df)+theme(legend.position="none")+geom_text(data=n.dat,aes(x=x-1.25,y=y,label=paste0(name,": ",dat)),size=2,hjust=0)+line.2pt("nexus","AwIns",geom.df,col="darkgreen")+line.2pt("resOr","nexus",geom.df,col="red")+line.2pt("mand.resOr","nexus",geom.df,col="gray")+line.2pt("inlever.p","nexus",geom.df,col="gray")+xlim(-4,4)+geom_point(aes(x=x,y=y),dat=circs,size=0.5)
    )
    
    if(out=="png"){png(filename=paste(graph.dir,"/geom_closed",i,".j.png",sep=""),width=700,height=700)
      print(p)
      dev.off()
    }else{
      if(file.exists(paste(graph.dir,"/geom_closed",i,".j.pdf",sep=""))) file.remove(paste(graph.dir,"/geom_closed",i,".j.pdf",sep=""))
      pdf(file=paste(graph.dir,"/geom_closed",i,".j.pdf",sep=""),width=7,height=7)
      suppressWarnings(print(p))
      dev.off()}
  }
    
}

setwd("~/Documents/Mikhaila")


##############