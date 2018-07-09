data.frame(
Iteration=0,
Time=Start,
AngAcc=0,
AngVel=0,
JawAng=Gape,
dJawAng=0,
dTendonTheta=0,
TendonTheta=TendonThetaOpen,
TendonLength=TendonLengthClosed,
TendonAngAcc=0,
TendonAngVel=0,
NormalDrag=0,
LooseDrag=0,
DragT=0,
PreyDrag=0,
PressureT=0,

###########  A3
A3ML=A3MLOpen,
A3dML=0,
A3v=0,
A3Penn=0,
A3FL=A3FLOpen, #Instantaneous FL
A3Fmax=0,
A3Ffv=1, #Force-velocity factor
A3Ffl=1, #Force-length factor
A3Fact=0, #Activation rise-time factor
A3Fpar=ifelse(A3MLOpen>A3MLOpt,A3Fmax*exp(2*log(1+MaxFpar)*(A3MLOpen/A3MLOpt-1))-A3Fmax,0), #Parallel elastic force
A3Force=0, #Force of the A3 division
A3Theta=A3MuscleThetaOpen, #A3 angle of instertion on lower jaw;
A3TendonTorque=0, #Torque imparted by A3 on tendon;
A3Torque=0, #Torque imparted by A3 on lower jaw;
A3Fout=0,
A3TendonTheta=A3TendonThetaOpen, #angle of A3 relative to tendon
A3AwTheta=A3AwThetaOpen, #angle of A3 relative to Aw

################ A2 vectors
A2ML=A2MLOpen,
A2dML=0,
A2v=0, #Muscle velocity in ml/s
A2Penn=A2PennOpen, #Instantaneous pennation angle
A2FL=A2FLOpen, #Instantaneous FL
A2Fmax=0, #Max force produced by A2
A2Ffv=1, #Force-velocity factor
A2Ffl=1, #Force-length factor
A2Fact=0, #Activation rise-time factor
A2Fpar=ifelse(A2MLOpen>A2MLOpt,A2Fmax*exp(2*log(1+MaxFpar)*(A2MLOpen/A2MLOpt-1))-A2Fmax,0), #Parallel elastic force
A2Force=0, #Force of the A2 division
A2Theta=A2MuscleThetaOpen, #A2 angle of instertion on lower jaw;
A2TendonTorque=0, #Torque imparted by A2 on tendon;
A2Torque=0, #Torque imparted by A2 on tendon;
A2Fout=0,
A2TendonTheta=A2TendonThetaOpen, #angle of A2 relative to tendon
A2AwTheta=A2AwThetaOpen, #angle of A3 relative to Aw

#Aw vectorsAw
AwML=AwMLOpen,
AwdML=0,
Awv=0, #Muscle velocity in ml/s
AwPenn=AwPennOpen, #Instantaneous pennation angle
AwFL=AwFLOpen, #Instantaneous FL
AwFmax=0, #Max force produced by Aw
AwFfv=1, #Force-velocity factor
AwFfl=1, #Force-length factor
AwFact=0,#Activation rise-time factor
AwFpar=ifelse(AwMLOpen>AwMLOpt,AwFmax*exp(2*log(1+MaxFpar)*(AwMLOpen/AwMLOpt-1))-AwFmax,0), #Parallel elastic force
AwForce=0,#Force of the Aw division
AwTheta=AwMuscleThetaOpen, #Aw angle of instertion on lower jaw;
AwTendonTorque=0, #Torque imparted by Aw on tendon;
AwTorque=0, #Torque imparted by Aw on lower jaw;
AwFout=0,
TendFout=0,
AwTendonTheta=AwTendonThetaOpen, #angle of Aw relative to tendon

#A2/3, Aw vectors
A2A3Fsum=0,


###Tendon
TendonTorque=0,
#vectors for static bite force and velocity
StaticBite=0,
JawVelTip=0
)
