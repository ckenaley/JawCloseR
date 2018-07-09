# JawCloseR

JawCloseR is an assembly of two R scripts that, with a few other text files, models lower-jaw adduction in fishes using a dynamic equilibrium model based on that of van Wassenbergh et al. (2005) and Kenaley (2012). The "jaw.model" script defines the function "run.jaw.prey.2" commences simulations. Each row of an input file (see example 

function(spec.n=3,dat=NULL,config=1,Loose=F,min.ML=0.6,OutPutVerb=T,prey=F,press=100,prey.per=NULL,prey.strike.ang=NULL,prey.pos="flat",MaxIts=1000,progress=F,print.every=10,inset=T,out="pdf")
