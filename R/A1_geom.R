if(F){theta=pi/40
open.angle <- -1*(x$start+j[n,"JawAng"]-x$SpecAng-theta)
post=A1TendonPostSeg
ant=x$A1TendonAntSeg
ant_p=x$A1AntSeg
inlev_l=A1Li
ins=A1Ins.open
inlev_p=point_ang_r(x$joint.p,x$A1Li,open.angle)

A1_close <- function(ant_p,inlev_p,ant,post){
TendonL <- dist_2d(ant_p,inlev_p)
dTheta <- cos_ang(post,TendonL,ant)
A1_ins_p <- point_ang_r(inlev_p,post,asin((ant_p[2]-inlev_p[2])/TendonL)+dTheta)

return(A1_ins_p)
}

theta=pi/40

d <- tibble(theta=theta*seq(0,9,0.5)) %>%
  mutate(open.angle=-1*(x$start+j[n,"JawAng"]-x$SpecAng)+theta)



JawTip<- point_ang_r(x$joint.p,x$JL,open.angle)

pt.l <- list()
for(i in d$open.angle){
  pt.i <- A1_close(ant_p=x$A1AntSeg,inlev_p=point_ang_r(x$joint.p,x$A1Li,i),ant=x$A1TendonAntSeg,post=x$A1TendonPostSeg)
  inp <- point_ang_r(x$joint.p,x$A1Li,i)
  jaw <- point_ang_r(x$joint.p,x$JL,i)
  A1Tins <- x$A1AntSeg
  pt.l[[paste(i)]] <- tibble(angle=i %>% deg,pt="A1ins",x=pt.i[1],y=pt.i[2]) %>%
    add_row(angle=i %>% deg,pt="A1Inlever",x=inp[1],y=inp[2]) %>%
    add_row(angle=i %>% deg,pt="JawTip",x=jaw[1],y=jaw[2]) %>%
  add_row(angle=i %>% deg,pt="joint",x=0,y=0)%>%
    add_row(angle=i %>% deg,pt="A1TendonAntIns",x= A1Tins[1],y= A1Tins[2])

}

pts <- do.call(rbind,pt.l)

for(a in d$open.angle){
p <- pts %>%
  filter(angle==a %>% deg) %>%
  ggplot(aes(x,y,col=pt))+geom_point()+coord_fixed() +ggtitle(label = a %>% deg)
print(p)
}
A1_close(ant_p=x$A1AntSeg,inlev_p=point_ang_r(x$joint.p,x$A1Li,open.angle),ant=x$A1TendonAntSeg,post=x$A1TendonPostSeg)


f=10
#### Sling load calculations

ant=x$A1TendonAntSeg
post=x$A1TendonPostSeg
tend=l+r-5

#confirmed with: https://www.northernstrands.com/sling-length-calculator.aspx
A1_sling <- function(post,ant,ant_p,inlev_p,f){

  span <- dist_2d(ant_p,inlev_p)
  theta <- cos_ang(post,span,ant)
  h <- sin(theta)*post #sin(theta)=o/hyp=h/post
  d1 <- h/tan(theta)#tan(theta)=o/a=h/d1
  d2 <- span-d1


  f_out <- (f*d2*post)/(h*(span))

return(f_out)
  }

#based on.confirmed w: https://www.omnicalculator.com/physics/tension

#T₁ = W / [cos(α) × sin(β) / cos(β) + sin(α)]

orig=x$A1Or
ant_p
ins_d <-pts %>% filter(pt=="A1ins")
ins=ins_d[,c("x","y")][2,]
ant


A1_tension <- function(orig_p,ins_p,inlev_p,ant_p,f){
  alpha <- abs(cos_ang(dist_2d(inlev_p,ins_p),dist_2d(orig_p,ins_p),dist_2d(inlev_p,orig_p)))
  if(alpha>pi/2) alpha <- alpha-pi/2

  beta <- cos_ang(dist_2d(orig_p,ins_p),dist_2d(ant_p,ins_p),dist_2d(orig_p,ant_p)) %>% abs
  if(beta>pi/2) beta <- beta-pi/2

  tension <- f/(cos(alpha)*sin(beta)/cos(beta)+sin(alpha))

  return(tension)
}
}
