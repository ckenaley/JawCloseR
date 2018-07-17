bw.theme <- function(begin=0,stop,axis.p=c(.85, .25),x.lab,y.lab,leg.text=2){
	theme.bw <-theme_bw(20)+theme(legend.title=element_blank(),legend.text=element_text(size=leg.text),panel.grid=element_blank(),legend.position=axis.p,legend.background =element_blank(),legend.key = element_blank(),axis.title.x = element_text(vjust=-0.5),axis.title.y = element_text(vjust=0.25))

theme.opts <- list(xlab(x.lab),ylab(y.lab),scale_x_continuous(limits=c(begin,stop)))
return(list(theme.bw,theme.opts))
}

bw.theme.box <- function(axis.p=c(.85, .25),x.lab,y.lab,leg.text=2){
	theme.bw <-theme_bw(20)+theme(legend.title=element_blank(),legend.text=element_text(size=leg.text),panel.grid=element_blank(),legend.position=axis.p,legend.background =element_blank(),legend.key = element_blank(),axis.title.x = element_text(vjust=-0.5),axis.title.y = element_text(vjust=0.25))

theme.opts <- list(xlab(x.lab),ylab(y.lab))
return(list(theme.bw,theme.opts))
}

#free axis
bw.theme.free <- function(axis.p=c(.85, .25),x.lab,y.lab,leg.text=2){
	theme.bw <-theme_bw(20)+theme(legend.title=element_blank(),legend.text=element_text(size=leg.text),panel.grid=element_blank(),legend.position=axis.p,legend.background =element_blank(),legend.key = element_blank(),axis.title.x = element_text(vjust=-0.5),axis.title.y = element_text(vjust=0.25))

theme.opts <- list(xlab(x.lab),ylab(y.lab))
return(list(theme.bw,theme.opts))
}
