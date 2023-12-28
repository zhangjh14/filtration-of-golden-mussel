
# calculate filtration rate -----------------------------------------------
cal_q<-function(L=seq(0.1,3,0.1),
                T=seq(1,30,1),
                ah=0.05,bh=2,
                aa=0.2,ba=1,
                ab=0.1,bb=1,
                ad=0.005,bd=2,
                al=0.04,bl=1,
                alpha=NULL,miu=NULL,
                V=0,H12=0){
  n<-length(L)
  m<-length(T)
  if (is.null(alpha)){
    alpha<-T^3*exp(-3/800*T*T)/1785.041
  }else{
    alpha<-rep(alpha,m)
  }
  if (is.null(miu)){
    miu<-(-0.0322*T+1.6927)*36
  }else{
    miu<-rep(miu,m)
  }
  
  if (n!=m){
    cat("为保证出图，L T 长度应一样")
  }else{
    Q<-matrix(0,ncol=n,nrow = m)
    H<-ah*L^bh
    g<-981*3600*3600
    
    sa<-aa*L^ba
    sb<-ab*L^bb
    sd<-ad*L^bd
    sl<-al*L^bl
    
    for (i in 1:n){
      for (j in 1:m){
        temp_a<-1.7*1.7/alpha[j]/alpha[j]/(3.14*sa[i]*sb[i])/(3.14*sa[i]*sb[i])
        temp_b<-256/3.14*sl[i]/(sd[i])^4*miu[j]
        temp_c<-c(-2*g*alpha[j]*H[i]-V*V*360000*360000+2*g*alpha[j]*H12)
        
        temp_delta<-temp_b*temp_b-4*temp_a*temp_c
        result<-(-temp_b+sqrt(temp_delta))/2/temp_a
        Q[i,j]<-result
      }
    }
  }
  Q
}


# calculation 2 -----------------------------------------------------------

cal_q1<-function(L=seq(0.1,3,0.1),
                T=seq(1,30,1),
                ah=0.05,bh=2,
                aA=0.02,bA=2,
                apump=10000000,bpump=-7,
                miu=NULL,aalpha=1/1785.041,balpha=-3/800,
                V=0,H12=0){

  
  alpha<-aalpha*T^3*exp(balpha*T*T)
  if (is.null(miu)){
    miu<-(-0.0322*T+1.6927)*36
  }else{
    miu<-rep(miu,m)
  }
  n<-length(L)
  Q<-rep(0,n)
  
  H<-ah*L^bh
  g<-981*3600*3600
  
  A<-aA*L^bA
  pump<-apump*L^bpump
  
  
  for (i in 1:n){
    
      temp_a<-1.7*1.7/alpha[i]/alpha[i]/(3.14*3.14*A[i]*A[i])
      temp_b<-256/3.14*pump[i]*miu[i]
      temp_c<-c(-2*g*alpha[i]*H[i]-V*V*360000*360000+2*g*alpha[i]*H12)
      
      temp_delta<-temp_b*temp_b-4*temp_a*temp_c
      result<-(-temp_b+sqrt(temp_delta))/2/temp_a
      Q[i]<-result
    
  }
  Q
}

# calculate 3 -------------------------------------------------------------

cal_q2<-function(L=seq(0.1,3,0.1),
                 T=seq(1,30,1),
                 ah=0.05,bh=2,
                 aA=0.02,bA=2,
                 apump=10000000,bpump=-7,
                 miu=NULL,aalpha=1/1785.041,balpha=-3/800,
                 V=0,H12=0){
  
  
  n<-length(L)
  m<-length(T)
  alpha<-T^3*exp(balpha*T*T)*aalpha
  
  
  if (is.null(miu)){
    miu<-(-0.0322*T+1.6927)*36
  }else{
    miu<-rep(miu,m)
  }
  
  if (n!=m){
    # cat(c("should same length"))
  }else{
    Q<-matrix(0,ncol=n,nrow = m)
    H<-ah*L^bh
    g<-981*3600*3600
    
    A<-aA*L^bA
    pump<-apump*L^bpump
    (m/s)
    for (i in 1:n){
      for (j in 1:m){
        temp_a<-1.7*1.7/alpha[j]/alpha[j]/(3.14*3.14*A[i]*A[i])
        temp_b<-256/3.14*pump[i]*miu[j]
        temp_c<-c(-2*g*alpha[j]*(H[i]-H12)-V*V*360000*360000)
        
        temp_delta<-temp_b*temp_b-4*temp_a*temp_c
        result<-(-temp_b+sqrt(temp_delta))/2/temp_a
        Q[i,j]<-result
      }
    }
  }
  Q
}
# sensitivity analysis ----------------------------------------------------
#1
sensitive<-function(x,flag,value="relative"){
  sensitive.data<-data.frame(x=x,
                             L20T05=0,
                             L20T10=0,
                             L20T15=0,
                             L20T20=0)
  w<-length(x)
  for(i in 1:w){
    if (flag=="ah"){
      Q<-cal_q(ah=sensitive.data$x[i])
    }else if(flag=="bh"){
      Q<-cal_q(bh=sensitive.data$x[i])
    }else if (flag=="aa"){
      Q<-cal_q(aa=sensitive.data$x[i])
    }else if (flag=="ba"){
      Q<-cal_q(ba=sensitive.data$x[i])
    }else if (flag=="ab"){
      Q<-cal_q(ab=sensitive.data$x[i])
    }else if (flag=="bb"){
      Q<-cal_q(bb=sensitive.data$x[i])
    }else if (flag=="ad"){
      Q<-cal_q(ad=sensitive.data$x[i])
    }else if (flag=="bd"){
      Q<-cal_q(bd=sensitive.data$x[i])
    }else if (flag=="al"){
      Q<-cal_q(al=sensitive.data$x[i])
    }else if (flag=="bl"){
      Q<-cal_q(bl=sensitive.data$x[i])
    }else if (flag=="alpha"){
      Q<-cal_q(alpha=sensitive.data$x[i])
    }else if (flag=="miu"){
      Q<-cal_q(miu=sensitive.data$x[i])
    }else if (flag=="V"){
      Q<-cal_q(V=sensitive.data$x[i])
    }else if (flag=="H12"){
      Q<-cal_q(H12=sensitive.data$x[i])
    }
    
    sensitive.data$L20T05[i]<-Q[20,5]
    sensitive.data$L20T10[i]<-Q[20,10]
    sensitive.data$L20T15[i]<-Q[20,15]
    sensitive.data$L20T20[i]<-Q[20,20]
    
  }
  if (value=="relative"){
    sensitive.data$L20T05<-(sensitive.data$L20T05-sensitive.data$L20T05[w/2])/sensitive.data$L20T05[w/2]
    sensitive.data$L20T10<-(sensitive.data$L20T10-sensitive.data$L20T10[w/2])/sensitive.data$L20T10[w/2]
    sensitive.data$L20T15<-(sensitive.data$L20T15-sensitive.data$L20T15[w/2])/sensitive.data$L20T15[w/2]
    sensitive.data$L20T20<-(sensitive.data$L20T20-sensitive.data$L20T20[w/2])/sensitive.data$L20T20[w/2]
  }
  sensitive.data
}

#2
sensitive2<-function(x,flag,value="relative"){
  sensitive.data<-data.frame(x=x,
                             L20T05=0,
                             L20T10=0,
                             L20T15=0,
                             L20T20=0)
  w<-length(x)
  for(i in 1:w){
    if (flag=="ah"){
      Q<-cal_q2(ah=sensitive.data$x[i])
    }else if(flag=="bh"){
      Q<-cal_q2(bh=sensitive.data$x[i])
    }else if (flag=="aA"){
      Q<-cal_q2(aA=sensitive.data$x[i])
    }else if (flag=="bA"){
      Q<-cal_q2(bA=sensitive.data$x[i])
    }else if (flag=="apump"){
      Q<-cal_q2(apump=sensitive.data$x[i])
    }else if (flag=="bpump"){
      Q<-cal_q2(bpump=sensitive.data$x[i])
    }else if (flag=="aalpha"){
      Q<-cal_q2(aalpha=sensitive.data$x[i])
    }else if (flag=="balpha"){
      Q<-cal_q2(balpha=sensitive.data$x[i])
    }else if (flag=="miu"){
      Q<-cal_q2(miu=sensitive.data$x[i])
    }else if (flag=="V"){
      Q<-cal_q2(V=sensitive.data$x[i])
    }else if (flag=="H12"){
      Q<-cal_q2(H12=sensitive.data$x[i])
    }
    
    sensitive.data$L20T05[i]<-Q[20,5]
    sensitive.data$L20T10[i]<-Q[20,10]
    sensitive.data$L20T15[i]<-Q[20,15]
    sensitive.data$L20T20[i]<-Q[20,20]
    
  }
  if (value=="relative"){
    sensitive.data$L20T05<-(sensitive.data$L20T05-sensitive.data$L20T05[w/2])/sensitive.data$L20T05[w/2]
    sensitive.data$L20T10<-(sensitive.data$L20T10-sensitive.data$L20T10[w/2])/sensitive.data$L20T10[w/2]
    sensitive.data$L20T15<-(sensitive.data$L20T15-sensitive.data$L20T15[w/2])/sensitive.data$L20T15[w/2]
    sensitive.data$L20T20<-(sensitive.data$L20T20-sensitive.data$L20T20[w/2])/sensitive.data$L20T20[w/2]
  }
  sensitive.data
}


# boxplot -----------------------------------------------------------------

myboxplot<-function(x,y,data,width=0.5){
  plotdata<-data.frame(X=data[[x]],y=data[[y]])
  colnames(plotdata)<-c("x","y")
  stats<-plotdata %>%
    group_by(x) %>%
    summarise(
      stats = list(boxplot.stats(y)$stats),

    )
  p1<-ggplot(data=plotdata)+
    geom_point(data=plotdata,aes(x=x,y=y))
  a<-width
  n<-length(unique(plotdata$x))
  for (i in 1:n){
    
    x<-stats$x[i]
    quantile<-stats$stats[[i]]
    p1<-p1+geom_segment(x=x,xend=x,y=quantile[1],yend=quantile[2])
    p1<-p1+geom_segment(x=x-a,xend=x+a,y=quantile[2],yend=quantile[2])
    p1<-p1+geom_segment(x=x-a,xend=x+a,y=quantile[3],yend=quantile[3])
    p1<-p1+geom_segment(x=x-a,xend=x+a,y=quantile[4],yend=quantile[4])
    p1<-p1+geom_segment(x=x,xend=x,y=quantile[4],yend=quantile[5])
    p1<-p1+geom_segment(x=x-a,xend=x-a,y=quantile[2],yend=quantile[4])
    p1<-p1+geom_segment(x=x+a,xend=x+a,y=quantile[2],yend=quantile[4])
  }
  p1

}
