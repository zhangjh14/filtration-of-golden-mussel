# get data ----------------------------------------------------------------

library(digitize)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggridges)
library(gghalves)

ReadAndCal("pestana fig/1.png")
digitize("pestana fig/10.png")

digitize("tokumon fig/gazulha.png")

digitize("sylvester fig/1.png")


# analysis ----------------------------------------------------------------
source("filtration model.r")
library(ggplot2)
library(tidyr)
library(leaps)
library(corrplot)
library(mgcv)
library(visreg)
library(pheatmap)
library(lme4)
library(vegan)
library(car)
library(effects)

mytheme<-function(){
    theme(axis.line = element_line(linetype = "solid"), 
          axis.ticks = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          panel.border = element_rect(fill = NA, 
                                      colour = "grey20"),
          panel.grid.major = element_line(colour = NA), 
          panel.grid.minor = element_line(colour = NA), 
          axis.title = element_text(size = 14, 
                                    face = "bold"), 
          axis.text = element_text(size = 12, 
                                   face = "bold", 
                                   colour = "black"), 
          axis.text.x = element_text(size = 12,
                                     vjust = -2), 
          axis.text.y = element_text(size = 12,
                                     hjust = -0.25), 
          plot.title = element_text(size = 15,
                                    face = "bold"), 
          panel.background = element_rect(fill = "white"), 
          legend.position = "right",
          plot.background = element_rect(colour = NA))
}

fil<-read.csv("fil.csv",header = TRUE)


fil%>% 
  subset(reference=='Tokumon et al. 2015')%>%
  ggplot(aes(x=length,y=rate,colour=T,shape=ft))+
  geom_point()

fil%>% 
  ggplot(aes(x=T,y=rate,colour=length,shape=reference))+
  geom_point()+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_text(size = 12,
                                   hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "right",
        plot.background = element_rect(colour = NA)) 

fil%>% 
  subset(reference=='Tokumon et al. 2015')%>%
  ggplot(aes(x=length,y=rate,colour=T,shape=ft))+
  geom_point()

fil%>%
  subset(reference=='Pestana, et al, 2009')%>%
  ggplot(aes(x=rate))+
  geom_histogram(binwidth = 10)+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_text(size = 12,
                                   hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "right",
        plot.background = element_rect(colour = NA))
 

m<-regsubsets(rate~T+time+ft+fc+sediment+V,data=fil,nbest=2)
plot(m,scale = "adjr2")
plot(m,scale = "bic")
fit_glm<-glm(rate~T+time+fc+V+length,family=poisson(link="log"),data=fil)
summary(fit_glm)

fit_gam<-gam(rate~s(T,k=3,by=ft)+s(time,k=3)+fc+V+length,data=fil)
summary(fit_gam)





# 实验 ----------------------------------------------------------------

# ##length ----------------------------------------------------------------
SL<-read.csv('length.csv')

ggplot(SL,aes(x=group,y=length))+
  geom_boxplot()+
  labs(x="Group",y="length (mm)")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_text(size = 12,
                                   hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "none",
        plot.background = element_rect(colour = NA)) 
ggsave("length.pdf",width=12,height=8,units="cm",dpi=600)


ggplot(SL,aes(x=group,y=weight))+
  geom_boxplot()+
  labs(x="Group",y="weight (g)")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_text(size = 12,
                                   hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "none",
        plot.background = element_rect(colour = NA)) 
ggsave("weight.pdf",width=12,height=8,units="cm",dpi=600)

SL %>%
  split(.$group) %>%
  sapply(function(x) c(mean=mean(x$length),sd=sd(x$length)))
SL %>%
  split(.$group) %>%
  sapply(function(x) c(mean=mean(x$weight),sd=sd(x$weight)))


#blank
blank<-read.csv('blank.csv')

blank_model<-lm(log(value)~time,data=subset(blank,group=='c'))
summary(blank_model)

ggplot(blank,aes(x=time,y=log(value),color=group))+
  geom_point()+
  geom_smooth(method = "lm")

test1<-read.csv('test1.csv')

ggplot(test1,aes(x=time,y=log(value),color=group))+
  geom_point()+
  geom_smooth(method = "lm")



t1<-c(0,10,20,30,50,70,90,120,180,240,300,360,600)
t2<-c(10,20,30,60,90,120,180,300,540)
t3<-c(0,10,20,30,60,90,120,180,300,540)

a<-read.csv('a.csv')
b<-read.csv('b.csv')
c<-read.csv('c.csv')

a1<-read.csv('a1.csv')
b1<-read.csv('b1.csv')
c1<-read.csv('c1.csv')

#ml/ind/h

afil1<-matrix(0,nrow=dim(a)[1],ncol=dim(a)[1])
afil2<-matrix(0,nrow=dim(a)[1],ncol=dim(a)[1])
afil3<-matrix(0,nrow=dim(a)[1],ncol=dim(a)[1])
for (i in 1:dim(a)[1]){
  for (j in 1:dim(a)[1]){
    if(i!=j){
      afil1[i,j]<-1000*(log(a$a1[i]/a$a1[j])/(a$time[j]-a$time[i])+0.003037)
      afil2[i,j]<-1000*(log(a$a2[i]/a$a2[j])/(a$time[j]-a$time[i])+0.003037)
      afil3[i,j]<-1000*(log(a$a3[i]/a$a3[j])/(a$time[j]-a$time[i])+0.003037)
    }
    
  }
}

pheatmap(afil3,cluster_rows = FALSE,
         cluster_cols = FALSE)

bfil1<-matrix(0,nrow=dim(b)[1],ncol=dim(b)[1])
bfil2<-matrix(0,nrow=dim(b)[1],ncol=dim(b)[1])
bfil3<-matrix(0,nrow=dim(b)[1],ncol=dim(b)[1])
for (i in 1:dim(b)[1]){
  for (j in 1:dim(b)[1]){
    if(i!=j){
      bfil1[i,j]<-1000*(log(b$b1[i]/b$b1[j])/(b$time[j]-b$time[i])+0.003037)
      bfil2[i,j]<-1000*(log(b$b2[i]/b$b2[j])/(b$time[j]-b$time[i])+0.003037)
      bfil3[i,j]<-1000*(log(b$b3[i]/b$b3[j])/(b$time[j]-b$time[i])+0.003037)
    }
  }
}

cfil1<-matrix(0,nrow=dim(c)[1],ncol=dim(c)[1])
cfil2<-matrix(0,nrow=dim(c)[1],ncol=dim(c)[1])
cfil3<-matrix(0,nrow=dim(c)[1],ncol=dim(c)[1])
for (i in 1:dim(c)[1]){
  for (j in 1:dim(c)[1]){
    if(i!=j){
      cfil1[i,j]<-1000*(log(c$c1[i]/c$c1[j])/(c$time[j]-c$time[i])+0.003037)
      cfil2[i,j]<-1000*(log(c$c2[i]/c$c2[j])/(c$time[j]-c$time[i])+0.003037)
      cfil3[i,j]<-1000*(log(c$c3[i]/c$c3[j])/(c$time[j]-c$time[i])+0.003037)
    }
    
  }
}


afil<-data.frame(a1=afil1[lower.tri(afil1,diag = FALSE)],
                 a2=afil2[lower.tri(afil2,diag = FALSE)],
                 a3=afil3[lower.tri(afil3,diag = FALSE)])

bfil<-data.frame(b1=bfil1[lower.tri(bfil1,diag = FALSE)],
                 b2=bfil2[lower.tri(bfil2,diag = FALSE)],
                 b3=bfil3[lower.tri(bfil3,diag = FALSE)])

cfil<-data.frame(c1=cfil1[lower.tri(cfil1,diag = FALSE)],
                 c2=cfil2[lower.tri(cfil2,diag = FALSE)],
                 c3=cfil3[lower.tri(cfil3,diag = FALSE)])

af<-pivot_longer(afil,a1:a3,names_to = "group")
bf<-pivot_longer(bfil,b1:b3,names_to = "group")
cf<-pivot_longer(cfil,c1:c3,names_to = "group")
fil<-rbind(af,bf,cf)
fil<-fil[which(fil$value>0),]
fil$length<-fil$value

fil$length[which(fil$group=="a1")]<-20.5
fil$length[which(fil$group=="a2")]<-14
fil$length[which(fil$group=="a3")]<-10
fil$length[which(fil$group=="b1")]<-18
fil$length[which(fil$group=="b2")]<-13
fil$length[which(fil$group=="b3")]<-9
fil$length[which(fil$group=="c1")]<-17
fil$length[which(fil$group=="c2")]<-16
fil$length[which(fil$group=="c3")]<-8

fil$value[which(fil$group=="a1")]<-fil$value[which(fil$group=="a1")]/10
fil$value[which(fil$group=="a2")]<-fil$value[which(fil$group=="a2")]/10
fil$value[which(fil$group=="a3")]<-fil$value[which(fil$group=="a3")]/15
fil$value[which(fil$group=="b1")]<-fil$value[which(fil$group=="b1")]/10
fil$value[which(fil$group=="b2")]<-fil$value[which(fil$group=="b2")]/9
fil$value[which(fil$group=="b3")]<-fil$value[which(fil$group=="b3")]/15
fil$value[which(fil$group=="c1")]<-fil$value[which(fil$group=="c1")]/10
fil$value[which(fil$group=="c2")]<-fil$value[which(fil$group=="c2")]/9
fil$value[which(fil$group=="c3")]<-fil$value[which(fil$group=="c3")]/13

fil$value[which(fil$group=="a1")]<-fil$value[which(fil$group=="a1")]/0.84
fil$value[which(fil$group=="a2")]<-fil$value[which(fil$group=="a2")]/0.41
fil$value[which(fil$group=="a3")]<-fil$value[which(fil$group=="a3")]/0.12
fil$value[which(fil$group=="b1")]<-fil$value[which(fil$group=="b1")]/0.58
fil$value[which(fil$group=="b2")]<-fil$value[which(fil$group=="b2")]/0.30
fil$value[which(fil$group=="b3")]<-fil$value[which(fil$group=="b3")]/0.09
fil$value[which(fil$group=="c1")]<-fil$value[which(fil$group=="c1")]/0.52
fil$value[which(fil$group=="c2")]<-fil$value[which(fil$group=="c2")]/0.24
fil$value[which(fil$group=="c3")]<-fil$value[which(fil$group=="c3")]/0.07
ggplot(fil1)+
  geom_boxplot(aes(x=as.factor(length),y=value*60),outlier.shape = 1)+
  coord_cartesian(ylim =  c(0, 180))

ggplot(fil1)+
  geom_point(aes(x=group,y=value*60))+
  coord_cartesian(ylim =  c(0, 180))

compare_means(value~group,data=fil,method = 't.test')

ka_model<-lm(log(blank)~time,data=a)
a1_model<-lm(log(a1)~time,data=a)
a2_model<-lm(log(a2)~time,data=a)
a3_model<-lm(log(a3)~time,data=a)

plot(x=a$time,y=log(a$blank))

kb_model<-lm(log(blank)~time,data=b)
b1_model<-lm(log(b1)~time,data=b)
b2_model<-lm(log(b2)~time,data=b)
b3_model<-lm(log(b3)~time,data=b)

1000*(kb_model$coefficients[2]-b1_model$coefficients[2])
1000*(kb_model$coefficients[2]-b2_model$coefficients[2])
1000*(kb_model$coefficients[2]-b3_model$coefficients[2])

kc_model<-lm(log(blank)~time,data=c)
c1_model<-lm(log(c1)~time,data=c)
c2_model<-lm(log(c2)~time,data=c)
c3_model<-lm(log(c3)~time,data=c)

1000*(kc_model$coefficients[2]-c1_model$coefficients[2])
1000*(kc_model$coefficients[2]-c2_model$coefficients[2])
1000*(kc_model$coefficients[2]-c3_model$coefficients[2])

                  
afil11<-matrix(0,nrow=dim(a1)[1],ncol=dim(a1)[1])
afil12<-matrix(0,nrow=dim(a1)[1],ncol=dim(a1)[1])
afil13<-matrix(0,nrow=dim(a1)[1],ncol=dim(a1)[1])
for (i in 1:dim(a1)[1]){
  for (j in 1:dim(a1)[1]){
    if(i!=j){
      afil11[i,j]<-1000*(log(a$a1[i]/a$a1[j])/(a$time[j]-a$time[i])+0.007423)
      afil12[i,j]<-1000*(log(a$a2[i]/a$a2[j])/(a$time[j]-a$time[i])+0.007423)
      afil13[i,j]<-1000*(log(a$a3[i]/a$a3[j])/(a$time[j]-a$time[i])+0.007423)
    }
    
  }
}                
pheatmap(afil11)

ka1_model<-lm(log(blank)~time,data=a1)
a11_model<-lm(log(a1)~time,data=a1)
a12_model<-lm(log(a2)~time,data=a1)
a13_model<-lm(log(a3)~time,data=a1)

1000*(ka1_model$coefficients[2]-a11_model$coefficients[2])
1000*(ka1_model$coefficients[2]-a12_model$coefficients[2])
1000*(ka1_model$coefficients[2]-a13_model$coefficients[2])


kb1_model<-lm(log(blank)~time,data=b1)
b11_model<-lm(log(b1)~time,data=b1)
b12_model<-lm(log(b2)~time,data=b1)
b13_model<-lm(log(b3)~time,data=b1)

1000*(kb1_model$coefficients[2]-b11_model$coefficients[2])
1000*(kb1_model$coefficients[2]-b12_model$coefficients[2])
1000*(kb1_model$coefficients[2]-b13_model$coefficients[2])

kc1_model<-lm(log(blank)~time,data=c1)
c11_model<-lm(log(c1)~time,data=c1)
c12_model<-lm(log(c2)~time,data=c1)
c13_model<-lm(log(c3)~time,data=c1)

1000*(kc1_model$coefficients[2]-c11_model$coefficients[2])
1000*(kc1_model$coefficients[2]-c12_model$coefficients[2])
1000*(kc1_model$coefficients[2]-c13_model$coefficients[2])

afil11<-matrix(0,nrow=dim(a1)[1],ncol=dim(a1)[1])
afil12<-matrix(0,nrow=dim(a1)[1],ncol=dim(a1)[1])
afil13<-matrix(0,nrow=dim(a1)[1],ncol=dim(a1)[1])
for (i in 1:dim(a1)[1]){
  for (j in 1:dim(a1)[1]){
    if(i!=j){
      afil11[i,j]<-1000/(a1$time[j]-a1$time[i])*(log(a1$a1[i]/a$a1[j])-log(a1$blank[i]/a1$blank[j]))
      afil12[i,j]<-1000/(a1$time[j]-a1$time[i])*(log(a1$a2[i]/a$a2[j])-log(a1$blank[i]/a1$blank[j]))
      afil13[i,j]<-1000/(a1$time[j]-a1$time[i])*(log(a1$a3[i]/a$a3[j])-log(a1$blank[i]/a1$blank[j]))
    }
    
  }
}


bfil11<-matrix(0,nrow=dim(b1)[1],ncol=dim(b1)[1])
bfil12<-matrix(0,nrow=dim(b1)[1],ncol=dim(b1)[1])
bfil13<-matrix(0,nrow=dim(b1)[1],ncol=dim(b1)[1])
for (i in 1:dim(b1)[1]){
  for (j in 1:dim(b1)[1]){
    if(i!=j){
      bfil11[i,j]<-1000*(log(b$b1[i]/b$b1[j])/(b$time[j]-b$time[i])+0.007423)
      bfil12[i,j]<-1000*(log(b$b2[i]/b$b2[j])/(b$time[j]-b$time[i])+0.007423)
      bfil13[i,j]<-1000*(log(b$b3[i]/b$b3[j])/(b$time[j]-b$time[i])+0.007423)
      
    }
    
  }
}

cfil11<-matrix(0,nrow=dim(c1)[1],ncol=dim(c1)[1])
cfil12<-matrix(0,nrow=dim(c1)[1],ncol=dim(c1)[1])
cfil13<-matrix(0,nrow=dim(c1)[1],ncol=dim(c1)[1])
for (i in 1:dim(c1)[1]){
  for (j in 1:dim(c1)[1]){
    if(i!=j){
      cfil11[i,j]<-1000*(log(c$c1[i]/c$c1[j])/(c$time[j]-c$time[i])+0.007423)
      cfil12[i,j]<-1000*(log(c$c2[i]/c$c2[j])/(c$time[j]-c$time[i])+0.007423)
      cfil13[i,j]<-1000*(log(c$c3[i]/c$c3[j])/(c$time[j]-c$time[i])+0.007423)
    }
    
  }
}


a1fil<-data.frame(a1=afil11[lower.tri(afil11,diag = FALSE)],
                 a2=afil12[lower.tri(afil12,diag = FALSE)],
                 a3=afil13[lower.tri(afil13,diag = FALSE)])

b1fil<-data.frame(b1=bfil11[lower.tri(bfil11,diag = FALSE)],
                 b2=bfil12[lower.tri(bfil12,diag = FALSE)],
                 b3=bfil13[lower.tri(bfil13,diag = FALSE)])

c1fil<-data.frame(c1=cfil11[lower.tri(cfil11,diag = FALSE)],
                 c2=cfil12[lower.tri(cfil12,diag = FALSE)],
                 c3=cfil13[lower.tri(cfil13,diag = FALSE)])

af1<-pivot_longer(afil,a1:a3,names_to = "group")
bf1<-pivot_longer(bfil,b1:b3,names_to = "group")
cf1<-pivot_longer(cfil,c1:c3,names_to = "group")
fil1<-rbind(af1,bf1,cf1)
fil1<-fil[which(fil1$value>0),]
fil1$length<-fil1$value

fil1$length[which(fil1$group=="a1")]<-20.5
fil1$length[which(fil1$group=="a2")]<-14
fil1$length[which(fil1$group=="a3")]<-10
fil1$length[which(fil1$group=="b1")]<-18
fil1$length[which(fil1$group=="b2")]<-13
fil1$length[which(fil1$group=="b3")]<-9
fil1$length[which(fil1$group=="c1")]<-17
fil1$length[which(fil1$group=="c2")]<-16
fil1$length[which(fil1$group=="c3")]<-8

fil1$value[which(fil1$group=="a1")]<-fil1$value[which(fil1$group=="a1")]/4
fil1$value[which(fil1$group=="a2")]<-fil1$value[which(fil1$group=="a2")]/8
fil1$value[which(fil1$group=="a3")]<-fil1$value[which(fil1$group=="a3")]/6
fil1$value[which(fil1$group=="b1")]<-fil1$value[which(fil1$group=="b1")]/7
fil1$value[which(fil1$group=="b2")]<-fil1$value[which(fil1$group=="b2")]/8
fil1$value[which(fil1$group=="b3")]<-fil1$value[which(fil1$group=="b3")]/7
fil1$value[which(fil1$group=="c1")]<-fil1$value[which(fil1$group=="c1")]/8
fil1$value[which(fil1$group=="c2")]<-fil1$value[which(fil1$group=="c2")]/5
fil1$value[which(fil1$group=="c3")]<-fil1$value[which(fil1$group=="c3")]/11

fil1$value[which(fil1$group=="a1")]<-fil1$value[which(fil1$group=="a1")]/0.84
fil1$value[which(fil1$group=="a2")]<-fil1$value[which(fil1$group=="a2")]/0.41
fil1$value[which(fil1$group=="a3")]<-fil1$value[which(fil1$group=="a3")]/0.12
fil1$value[which(fil1$group=="b1")]<-fil1$value[which(fil1$group=="b1")]/0.58
fil1$value[which(fil1$group=="b2")]<-fil1$value[which(fil1$group=="b2")]/0.30
fil1$value[which(fil1$group=="b3")]<-fil1$value[which(fil1$group=="b3")]/0.09
fil1$value[which(fil1$group=="c1")]<-fil1$value[which(fil1$group=="c1")]/0.52
fil1$value[which(fil1$group=="c2")]<-fil1$value[which(fil1$group=="c2")]/0.24
fil1$value[which(fil1$group=="c3")]<-fil1$value[which(fil1$group=="c3")]/0.07

fil1<-fil1[!is.na(fil1$group),]


dev.new()
ggplot(rbind(fil,fil1))+
  geom_boxplot(aes(x=as.factor(length),y=60*value),outlier.shape = NA)+
  coord_cartesian(ylim =  c(0, 3))

ggplot(fil1)+
  geom_boxplot(aes(x=as.factor(length),y=60*value),outlier.shape = NA)+
  coord_cartesian(ylim =  c(0,100))+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "black",
                                    linetype = "solid"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 12, 
                                  face = "bold"), 
        axis.text = element_text(size = 10, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 10,
                                   vjust = -2), 
        axis.text.y = element_text(size = 10,
                                   hjust = -0.25), 
        panel.background = element_rect(fill = "white"),
        legend.text = element_text(size = 12,
                                   face = "bold"), 
        legend.title = element_text(size = 14, 
                                    face = "bold"),
        legend.position = c(0.7,0.7), #比例数字
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)
  ) 

ggplot(fil)+
  geom_boxplot(aes(x=as.factor(length),y=60*value),outlier.shape = NA)+
  coord_cartesian(ylim =  c(0, 100))+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "black",
                                    linetype = "solid"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 12, 
                                  face = "bold"), 
        axis.text = element_text(size = 10, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 10,
                                   vjust = -2), 
        axis.text.y = element_text(size = 10,
                                   hjust = -0.25), 
        panel.background = element_rect(fill = "white"),
        legend.text = element_text(size = 12,
                                   face = "bold"), 
        legend.title = element_text(size = 14, 
                                    face = "bold"),
        legend.position = c(0.7,0.7), #比例数字
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)
  )

colnames<-c("group","length","weight","num","deltaT","ci","cf")

t<-data.frame(matrix(ncol=length(colnames),nrow=0))
colnames(t)<-colnames

for (i in 1:dim(a)[1]){
  for (j in i:dim(a)[1]){
    if(i!=j){
      new_row1<-data.frame(group="a1",
                           length=21,
                           weight=0.839,
                           num=10,
                           deltat=a$time[j]-a$time[i],
                           ci=a$a1[i],
                           cf=a$a1[j])
      new_row2<-data.frame(group="a2",
                           length=16.49,
                           weight=0.414,
                           num=10,
                           deltat=a$time[j]-a$time[i],
                           ci=a$a2[i],
                           cf=a$a2[j])
      new_row3<-data.frame(group="a3",
                           length=10.48,
                           weight=0.117,
                           num=15,
                           deltat=a$time[j]-a$time[i],
                           ci=a$a3[i],
                           cf=a$a3[j])
      t<-rbind(t,new_row1)
      t<-rbind(t,new_row2)
      t<-rbind(t,new_row3)
      
      
    }
    
  }
}

for (i in 1:dim(b)[1]){
  for (j in i:dim(b)[1]){
    if(i!=j){
      new_row1<-data.frame(group="b1",
                           length=18.57,
                           weight=0.582,
                           num=10,
                           deltat=b$time[j]-b$time[i],
                           ci=b$b1[i],
                           cf=b$b1[j])
      new_row2<-data.frame(group="b2",
                           length=14.389,
                           weight=0.302,
                           num=9,
                           deltat=b$time[j]-b$time[i],
                           ci=b$b2[i],
                           cf=b$b2[j])
      new_row3<-data.frame(group="b3",
                           length=9.62,
                           weight=0.0953,
                           num=15,
                           deltat=b$time[j]-b$time[i],
                           ci=b$b3[i],
                           cf=b$b3[j])
      t<-rbind(t,new_row1)
      t<-rbind(t,new_row2)
      t<-rbind(t,new_row3)
      
      
    }
    
  }
}

for (i in 1:dim(c)[1]){
  for (j in i:dim(c)[1]){
    if(i!=j){
      new_row1<-data.frame(group="c1",
                           length=17.669,
                           weight=0.519,
                           num=10,
                           deltat=c$time[j]-c$time[i],
                           ci=c$c1[i],
                           cf=c$c1[j])
      new_row2<-data.frame(group="c2",
                           length=13.469,
                           weight=0.238,
                           num=9,
                           deltat=c$time[j]-c$time[i],
                           ci=c$c2[i],
                           cf=c$c2[j])
      new_row3<-data.frame(group="c3",
                           length=8.42,
                           weight=0.068,
                           num=13,
                           deltat=c$time[j]-c$time[i],
                           ci=c$c3[i],
                           cf=c$c3[j])
      t<-rbind(t,new_row1)
      t<-rbind(t,new_row2)
      t<-rbind(t,new_row3)
      
      
    }
    
  }
}

t<-t[!is.na(t$group),]
write.csv(t,"t.csv")


t1<-data.frame(matrix(ncol=length(colnames),nrow=0))
colnames(t1)<-colnames

for (i in 1:dim(a1)[1]){
  for (j in i:dim(a1)[1]){
    if(i!=j){
      new_row1<-data.frame(group="a1",
                           length=21,
                           weight=0.839,
                           num=10,
                           deltat=a1$time[j]-a1$time[i],
                           ci=a1$a1[i],
                           cf=a1$a1[j])
      new_row2<-data.frame(group="a2",
                           length=16.49,
                           weight=0.414,
                           num=10,
                           deltat=a1$time[j]-a1$time[i],
                           ci=a1$a2[i],
                           cf=a1$a2[j])
      new_row3<-data.frame(group="a3",
                           length=10.48,
                           weight=0.117,
                           num=15,
                           deltat=a1$time[j]-a1$time[i],
                           ci=a1$a3[i],
                           cf=a1$a3[j])
      t1<-rbind(t1,new_row1)
      t1<-rbind(t1,new_row2)
      t1<-rbind(t1,new_row3)
      
      
    }
    
  }
}

for (i in 1:dim(b1)[1]){
  for (j in i:dim(b1)[1]){
    if(i!=j){
      new_row1<-data.frame(group="b1",
                           length=18.57,
                           weight=0.582,
                           num=10,
                           deltat=b1$time[j]-b1$time[i],
                           ci=b1$b1[i],
                           cf=b1$b1[j])
      new_row2<-data.frame(group="b2",
                           length=14.389,
                           weight=0.302,
                           num=9,
                           deltat=b1$time[j]-b1$time[i],
                           ci=b1$b2[i],
                           cf=b1$b2[j])
      new_row3<-data.frame(group="b3",
                           length=9.62,
                           weight=0.0953,
                           num=15,
                           deltat=b1$time[j]-b1$time[i],
                           ci=b1$b3[i],
                           cf=b1$b3[j])
      t1<-rbind(t1,new_row1)
      t1<-rbind(t1,new_row2)
      t1<-rbind(t1,new_row3)
      
      
    }
    
  }
}

for (i in 1:dim(c1)[1]){
  for (j in i:dim(c1)[1]){
    if(i!=j){
      new_row1<-data.frame(group="c1",
                           length=17.669,
                           weight=0.519,
                           num=10,
                           deltat=c1$time[j]-c1$time[i],
                           ci=c1$c1[i],
                           cf=c1$c1[j])
      new_row2<-data.frame(group="c2",
                           length=13.469,
                           weight=0.238,
                           num=9,
                           deltat=c1$time[j]-c1$time[i],
                           ci=c1$c2[i],
                           cf=c1$c2[j])
      new_row3<-data.frame(group="c3",
                           length=8.42,
                           weight=0.068,
                           num=13,
                           deltat=c1$time[j]-c1$time[i],
                           ci=c1$c3[i],
                           cf=c1$c3[j])
      t1<-rbind(t1,new_row1)
      t1<-rbind(t1,new_row2)
      t1<-rbind(t1,new_row3)
      
    }
    
  }
}

t1<-t1[!is.na(t1$group),]
write.csv(t1,"t1.csv")

t$filtration<-1000*(log(t$ci/t$cf)/t$deltat+0.003037)*60/t$num
t$test<-c("t")
t$fil1<-t$filtration/t$weight

t1$filtration<-1000*(log(t1$ci/t1$cf)/t1$deltat+0.007423)*60/t1$num
t1$test<-c("t1")
t1$fil1<-t1$filtration/t1$weight

t1 %>%
  rbind(t) %>%
  subset(filtration>10) %>%
  subset(deltat<=180) %>%
  ggplot()+
  geom_point(aes(x=deltat,y=filtration,color=test))

t1 %>%
  rbind(t) %>%
  subset(filtration>10) %>%
  subset(deltat==60) %>%
  ggboxplot(x="length",y="filtration",color = "test", palette = "jco",
            add = "jitter")+
  stat_compare_means(aes(group = test), label = "p.signif") ##label = "p.format"


###处理后的数据
t_60<-read.csv("t_60.csv")
t1_60<-read.csv("t1_60.csv")

t_60$filtration<-1000*(log(t_60$ci/t_60$cf)/t_60$deltat-0.003037)*60/t_60$num
t_60$test<-c("t")
t_60$fil1<-t_60$filtration/t_60$weight

t1_60$filtration<-1000*(log(t1_60$ci/t1_60$cf)/t1_60$deltat-0.007423)*60/t1_60$num
t1_60$test<-c("t1")
t1_60$fil1<-t1_60$filtration/t1_60$weight

t1_60 %>%
  rbind(t_60) %>%
  subset(filtration>10) %>%
  subset(deltat<=180) %>%
  ggplot()+
  geom_point(aes(x=deltat,y=filtration,color=test))+
  scale_color_manual(values = c("#d9b611","#0aa344"))+
  scale_x_continuous(breaks = seq(0,180, by = 60))+
  labs(x="Experiment time (min)",y="Filtration rate (ml ind.-1 h-1)")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_text(size = 12,
                                   hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "right",
        plot.background = element_rect(colour = NA)) 
ggsave("filtration-time.pdf",width=12,height=8,units="cm",dpi=600)

# food type:t1是藻，t是酵母 -----------------------------------------------------------------------
t1_60 %>%
  rbind(t_60) %>%
  subset(deltat==60) %>%
  ggboxplot(x="length",y="filtration",color = "test", palette = "jco",
            add = "jitter")+
  stat_compare_means(aes(group = length), label = "p.format",method="kruskal.test",paired=TRUE)+  ##label = "p.format"，"p.signif"
  scale_color_manual(values = c("#d9b611","#0aa344"))+
  labs(x="Group",y="Filtration rate (ml ind.-1 h-1)")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_text(size = 12,
                                   hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "right",
        plot.background = element_rect(colour = NA)) 




ggsave("filtration-food type.pdf",width=12,height=8,units="cm",dpi=600)

t1_60 %>%
  rbind(t_60) %>%
  subset(deltat==60) %>%
  ggboxplot(x="test",y="fil1",color = "test", palette = "jco",
            add = "jitter")+
  stat_compare_means(aes(group = test), label = "p.format",method="wilcox.test",paired=TRUE) 

max(subset(t_60,deltat==60)$fil1)

# wilcox.test kruskal.test 
group.difference<-compare_means(filtration~length,data=subset(t1_60,deltat==60),method = 'kruskal.test')

t1_60 %>%
  rbind(t_60) %>%
  subset(deltat==60) %>%
  subset(length==18.57) %>%
  compare_means(fil1~test,data=.,method = 'wilcox.test')

# shell length -----------------------------------------------------------------------

t1_60 %>%
  rbind(t_60) %>%
  subset(deltat==60) %>%
  # ggboxplot(x="length",y="filtration",color = "test", palette = "jco",
  #           add = "jitter",facet.by = "test")+
  ggboxplot(x="length",y="filtration",color = "test", palette = "jco",
            add = "jitter")+
  stat_compare_means(method = "anova")+  ##label = "p.format"，"p.signif"
  scale_color_manual(values = c("#00ACFE","#0aa344"))+
  labs(x="Group",y=NULL)+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_blank (), 
        # axis.text.y = element_text(size = 12,
        #                            hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "right",
        plot.background = element_rect(colour = NA)) ->p1
ggsave("filtration-length.pdf",width=24,height=8,units="cm",dpi=600)

t1_60 %>%
  rbind(t_60) %>%
  subset(deltat==60) %>%
  # ggboxplot(x="length",y="fil1",color = "test", palette = "jco",
  #           add = "jitter",facet.by = "test")+
  ggboxplot(x="length",y="fil1",color = "test", palette = "jco",
            add = "jitter")+
  stat_compare_means(method = "anova")+  ##label = "p.format"，"p.signif"
  scale_color_manual(values = c("#00ACFE","#0aa344"))+
  labs(x="Group",y=NULL)+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_blank (), 
        # axis.text.y = element_text(size = 12,
        #                            hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "right",
        plot.background = element_rect(colour = NA)) ->p2
ggsave("fil1-length.pdf",width=24,height=8,units="cm",dpi=600)

g<-grid.arrange(p1,p2,ncol=1)
ggsave("fil-length.pdf",g,width=14,height=8,units="cm",dpi=600)

t1_60 %>%
  rbind(t_60) %>%
  subset(deltat==60) %>%
  ggplot(aes(x=length,y=fil1,color=test))+
  geom_point()+
  labs(x="length (mm)",y="Filtration rate (ml ind.-1 h-1)")+
  scale_color_manual(values = c("#00ACFE","#0aa344"))+
  geom_smooth(method = "glm")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_text(size = 12,
                                   hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "right",
        plot.background = element_rect(colour = NA)) 
ggsave("fil1-length.pdf",width=12,height=8,units="cm",dpi=600)

cor(subset(t1_60,deltat==60)$length, subset(t1_60,deltat==60)$filtration, method = 'pearson')

cor.test(subset(t1_60,deltat==60)$length, subset(t1_60,deltat==60)$fil1, method = 'pearson')

cor.test(subset(t1_60,deltat==60)$length, subset(t1_60,deltat==60)$fil1, method = 'pearson')

t1_60 %>%
  rbind(t_60) %>%
  subset(deltat==60) %>%
  ggplot(aes(y=as.factor(length),x=filtration,fill=test,color=test))+
  geom_density_ridges(rel_min_height = 0.005, # 剪尾
                      scale = 1, # 山脊比例
                      alpha = 0.8, # 透明度
                      linetype = 1,# 脊线条类型
                      lwd = 0.2 )+
  scale_fill_manual(values = c("#00ACFE","#85BD41"))+
  scale_color_manual(values=c("#000080","#2F4F4F"))+
  mytheme()
ggsave("ridgelength.pdf",width=12,height=8,units="cm",dpi=600 )

t1_model<-nls(fil1~a*length^(b),data=subset(t1_60,deltat==60),start = list(a=53103,b=-2.27))
fil1_predict<-predict(t1_model)
RSS<-sum((subset(t1_60,deltat==60)$fil1-fil1_predict)^2)
TSS<-sum((subset(t1_60,deltat==60)$fil1-mean(subset(t1_60,deltat==60)$fil1))^2)
1-RSS/TSS


t_model<-nls(fil1~a*length^(b),data=subset(t_60,deltat==60),start = list(a=53103,b=-2.27))
summary(t_model)
fil1_predict<-predict(t_model)
RSS<-sum((subset(t_60,deltat==60)$fil1-fil1_predict)^2)
TSS<-sum((subset(t_60,deltat==60)$fil1-mean(subset(t_60,deltat==60)$fil1))^2)
1-RSS/TSS


# temperature -----------------------------------------------------------------------

temperature<-read.csv('temperature.csv')

temperature[!is.na(temperature$length),] %>%
  ggboxplot(x="group",y="length", palette = "jco",
            add = "jitter")+
  stat_compare_means(aes(group = group), label = "p.signif") 


temperature %>%
  split(.$group) %>%
  sapply(function(x) c(mean=mean(x$length,na.rm = TRUE),sd=sd(x$length,na.rm = TRUE)))
temperature %>%
  split(.$group) %>%
  sapply(function(x) c(mean=mean(x$weight,na.rm = TRUE),sd=sd(x$weight,na.rm = TRUE)))
temperature %>%
  split(.$group) %>%
  sapply(function(x) c(mean=mean(x$T),sd=sd(x$T)))

temperature$filtration<-300*(log(temperature$ci/temperature$cf)/temperature$deltat-temperature$settel)
temperature$fil1<-temperature$filtration/temperature$weight

temperature %>%
  ggboxplot(x="tem",y="filtration", palette = "jco",
            add = "jitter")+
  stat_compare_means(aes(group = tem), label = "p.signif")+ ##label = "p.format"
  labs(x="Temperature (℃)",y="Filtration rate (ml ind.-1 h-1)")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_text(size = 12,
                                   hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "right",
        plot.background = element_rect(colour = NA)) 
ggsave("filtration-T.pdf",width=12,height=8,units="cm",dpi=600)

temperature %>%
  ggboxplot(x="tem",y="fil1", palette = "jco",
            add = "jitter")+
  labs(x="Temperature (℃)",y="Filtration rate (ml ind.-1 h-1)")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_text(size = 12,
                                   hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "right",
        plot.background = element_rect(colour = NA)) 
ggsave("fil1-T.pdf",width=12,height=8,units="cm",dpi=600)



temperature %>%
  ggplot(aes(y=as.factor(tem),x=filtration,fill=as.factor(tem)))+
  geom_density_ridges(rel_min_height = 0.005, # 剪尾
                      scale = 1, # 山脊比例
                      alpha = 0.8, # 透明度
                      linetype = 1,# 脊线条类型
                      lwd = 0.2 )+
  scale_fill_manual(values = c("#FFE4C4","#FFDEAD","#F5DEB3","#D2B48C","#F4A460","#FF7F50","#FF4500"))+
  mytheme()
ggsave("ridgetem.pdf",width=12,height=8,units="cm",dpi=600 )


temperature %>%
  ggplot(aes(x=as.factor(tem),y=filtration)) + 
  ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(width = .05, alpha = .3)+
  labs(x="Temperature (℃)",y="Filtration rate (ml ind.-1 h-1)")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_text(size = 12,
                                   hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "right",
        plot.background = element_rect(colour = NA)) 
ggsave("fig3.pdf",width=10,height=6,units="cm",dpi=600 )

  

cor.test(subset(temperature,T>20)$T, subset(temperature,T>20)$filtration, method = 'pearson')

# velocity -----------------------------------------------------------------------

velocity<-read.csv("lentic.csv")

velocity$filtration<-1000*(log(velocity$ci/velocity$cf)/velocity$deltat-velocity$settel)/3
velocity$fil1<-velocity$filtration/velocity$weight

velocity %>%
  ggboxplot(x="group",y="fil1", palette = "jco",
            add = "jitter")+
  stat_compare_means(aes(group = group), label = "p.format")+ ##label = "p.format"
  labs(x="Group",y="Filtration rate (ml ind.-1 h-1)")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_text(size = 12,
                                   hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "right",
        plot.background = element_rect(colour = NA)) 
ggsave("fil1-v.pdf",width=12,height=8,units="cm",dpi=600)


# piv -----------------------------------------------------------------------
v<-read.csv("v.csv")
v$velocity<-v$velocity*8000 ##返回为像素速度

v$velocity<-v$velocity/8000*0.5 ##cm/s

v %>%
  subset(velocity<4000) %>%
  
  ggplot(aes(x=velocity))+
  geom_histogram()
sd(subset(v,velocity<4000)$velocity)*1200*800*3.14*3600/((2000)^3) ##1cm=2000pix

mean(v$velocity)*0.5*0.2*3.14*3600

v %>%
  subset(velocity<4000) %>%
  split(.$line) %>%
  sapply(function(x) c(mean=mean(x$velocity)*1200*800*3.14*3600/((2000)^3),
                       sd=sd(x$velocity))*1200*800*3.14*3600/((2000)^3)) -> v_mean
write.csv(v_mean,"v_mean.csv")

piv<-read.csv("v_mean.csv")

piv %>%
  ggboxplot(x="group",y="filtration", palette = "jco",
            add = "jitter")+
  stat_compare_means(aes(group = group), label = "p.format")+ ##label = "p.format"
  labs(x="Method",y="Filtration rate (ml ind.-1 h-1)")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_text(size = 12,
                                   hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "right",
        plot.background = element_rect(colour = NA)) 
ggsave("piv.pdf",width=12,height=8,units="cm",dpi=600)

mean(subset(piv,group=="direct")$filtration)

# summary -----------------------------------------------------------------------

colnames(t1_60)
colnames(temperature)
colnames(velocity)

s_t<-data.frame(flow=rep("lentic",length(t_60$group)),
                temperature=rep(25,length(t_60$group)),
                deltat=t_60$deltat/60,
                length=t_60$length,
                weight=t_60$weight,
                foodtype=rep("c",length(t_60$group)),
                filtration=t_60$filtration,
                fil1=t_60$fil1)

s_t1<-data.frame(flow=rep("lentic",length(t1_60$group)),
                temperature=rep(25,length(t1_60$group)),
                deltat=t1_60$deltat/60,
                length=t1_60$length,
                weight=t1_60$weight,
                foodtype=rep("s",length(t1_60$group)),
                filtration=t1_60$filtration,
                fil1=t1_60$fil1)

s_tem<-data.frame(flow=rep("lentic",length(temperature$group)),
                 temperature=temperature$tem,
                 deltat=temperature$deltat,
                 length=temperature$length,
                 weight=temperature$weight,
                 foodtype=rep("s",length(temperature$group)),
                 filtration=temperature$filtration,
                 fil1=temperature$fil1)

s_v<-data.frame(flow=velocity$group,
                  temperature=rep(18,length(velocity$group)),
                  deltat=velocity$deltat,
                  length=velocity$length,
                  weight=velocity$weight,
                  foodtype=rep("s",length(velocity$group)),
                  filtration=velocity$filtration,
                  fil1=velocity$fil1)

all_fil<-rbind(s_t,s_t1,s_tem,s_v)

all_fil %>%
  subset(deltat==1) %>%
  ggplot()+
  geom_point(aes(x=temperature,y=filtration))
  ggboxplot(x="temperature",y="filtration",color = "foodtype", palette = "jco",
            add = "jitter")
all_fil1<-subset(all_fil,fil1>0)

ggplot(all_fil1,aes(sample = log(fil1))) +
  geom_qq() + geom_qq_line()

varpart
# glm -----------------------------------------------------------------------

model_glm<-glm(filtration~temperature+flow+temperature+length+weight,
               family =Gamma(link = "inverse"),
               data=all_fil1) 
glm_step <- step(model_glm, direction = "both")
summary(model_glm) 
summary(glm_step) 

#伪r2
with(summary(model_glm), 1 - deviance/null.deviance)

newdata<-data.frame(foodtype=all_fil1$foodtype,
                    deltat=all_fil1$deltat,
                    flow=all_fil1$flow,
                    length=all_fil1$length,
                    weight=all_fil1$weight)
predicted <- predict(glm_step, newdata = newdata)
# Plot predicted vs. observed values
pre_ob<-data.frame(ob=all_fil1$filtration,
                   predict=glm_step$fit,
                   re=residuals(glm_step)
                   )
ggplot(pre_ob)+
  geom_point(aes(x=predict,y=ob))+
  geom_abline(slope=1,intercept = 0)+
  labs(x="predicted value (ml ind.-1 h-1)",y="observed value (ml ind.-1 h-1)")+
  scale_x_continuous(limits = c(0,500),expand = c(0,0))+
  scale_y_continuous(limits = c(0,500),expand = c(0,0))+
  coord_fixed(ratio=1)+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_text(size = 12,
                                   hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "right",
        plot.background = element_rect(colour = NA)) 
ggsave("pre-ob.pdf",width=8,height=8,units="cm",dpi=600)

ggplot(pre_ob)+
  geom_point(aes(x=predict,y=re))+
  geom_abline(slope=0,intercept = 0)+
  labs(x="predicted value (ml ind.-1 h-1)",y="predicted value (ml ind.-1 h-1)")+
  scale_x_continuous(limits = c(0,500),expand = c(0,0))+
  scale_y_continuous(limits = c(-4,4),expand = c(0,0))+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_text(size = 12,
                                   hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "right",
        plot.background = element_rect(colour = NA)) 
ggsave("pre-re.pdf",width=8,height=8,units="cm",dpi=600)  

residuals<-data.frame(residuals=residuals)
ggplot(residuals,aes(sample =log(residuals))) +
  geom_qq() +
  geom_qq_line()

effect<-allEffects(glm_step)
plot(allEffects(glm_step))

leaps <- regsubsets(filtration~foodtype+deltat+temperature+flow+temperature+length+weight,
                    data = all_fil1,
                    nbest=2)

plot(leaps, scale = "adjr2")

model_glm1<-glm(fil1~foodtype+log(deltat)+temperature+flow+length+weight,
                family =Gamma(link = "log"),,
                data=all_fil1) 
glm1_step <- step(model_glm1, direction = "both")
summary(model_glm1) 
summary(glm1_step) 


predicted <- predict(glm1_step, newdata = all_fil1)
# Plot predicted vs. observed values
plot(predicted, all_fil1$fil1)
abline(0,1)

residuals <- residuals(glm1_step)
plot(glm1_step$fit, residuals)
abline(h = 0)

plot(allEffects(model_glm))

df_resid <- augment(glm1_step, type.predict = "response")

# Plot the partial residuals
ggplot(df_resid, aes(x = .fitted, y = .resid)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red") +
  xlab(expression(hat(y))) +
  ylab("Partial Residuals") +
  ggtitle("Partial Residual Plot for GLM")
 
model_gam<-gam(filtration~foodtype+deltat+temperature+flow+length+weight,data=all_fil)
gam_step <- step(model_gam, direction = "both")
summary(model_gam) 
summary(gam_step) 

#文献
write.csv(all_fil1,"all_fil1.csv")
lit_fil1<-read.csv("fil1.csv",header = TRUE)

lit_fil1 %>%
  subset(filtration>0) %>%
  ggboxplot(x="temperature",y="filtration", palette = "jco",
          add = "jitter")+
    stat_compare_means(aes(group = foodtype), label = "p.format")

lit_fil1 %>%
  subset(filtration>0) %>%
  ggplot()+
  geom_point(aes(x=temperature,y=filtration))


# piv -----------------------------------------------------------------------

v<-read.csv("v.csv")
v$velocity<-v$velocity*8000 ##返回为像素速度

v$velocity<-v$velocity/8000*0.5 ##cm/s

v %>%
  subset(velocity<4000) %>%

  ggplot(aes(x=velocity))+
  geom_histogram()
sd(subset(v,velocity<4000)$velocity)*1200*800*3.14*3600/((2000)^3) ##1cm=2000pix

mean(v$velocity)*0.5*0.2*3.14*3600

v %>%
  subset(velocity<4000) %>%
  split(.$line) %>%
  sapply(function(x) c(mean=mean(x$velocity)*1200*800*3.14*3600/((2000)^3),
                       sd=sd(x$velocity))*1200*800*3.14*3600/((2000)^3)) -> v_mean




# modeling the pump -------------------------------------------------------

L<-seq(0.1,3,0.1)
T<-seq(1,30,1)

H<-L*L*1
g<-9810*3600*3600
alpha<-T^3*exp(-3/800*T*T)/1785.041
sa<-0.2*L
sb<-0.1*L
sd<-0.3*L
sl<-0.001*L
miu<-(-0.0322*T+1.6927)*36
V<-2

Q<-matrix(0,ncol=30,nrow = 30)

for (i in 1:30){
  for (j in 1:30){
    temp_a<-1.7*1.7/alpha[j]/alpha[j]/(3.14*sa[i]*sb[i])/(3.14*sa[i]*sb[i])
    temp_b<-256/3.14*sl[i]/(sd[i])^4*miu[j]
    temp_c<--2*g*alpha[j]*H[i]-V*V/2/g
    
    temp_delta<-temp_b*temp_b-4*temp_a*temp_c
    result<-(-temp_b+sqrt(temp_delta))/2/temp_a
    Q[i,j]<-result
  }
}



Q<-cal_q(L,T,
         ah=0.05,bh=2,
         aa=0.2,ba=1,
         ab=0.1,bb=1,
         ad=0.005,bd=2,
         al=0.04,bl=1,
         alpha=NULL,miu=NULL,
         V=0.1)
persp(L,T,Q, theta = 225, phi = 20,
      expand = 0.5,
      r=180,
      ltheta = 225,
      shade = 0.75,
      ticktype = "detailed")

# sensitive_ah<-data.frame(ah=seq(0.1,2,0.1),
#                          L20T05=0,
#                          L20T10=0,
#                          L20T15=0,
#                          L20T20=0)
# w<-20
# for (i in 1:w){
#   Q<-cal_q(L,T,ah=sensitive_ah$ah[i])
#   sensitive_ah$L20T05[i]<-Q[20,5]
#   sensitive_ah$L20T10[i]<-Q[20,10]
#   sensitive_ah$L20T15[i]<-Q[20,15]
#   sensitive_ah$L20T20[i]<-Q[20,20]
#   
# }
# sensitive_ah$L20T05<-(sensitive_ah$L20T05-sensitive_ah$L20T05[w/2])/sensitive_ah$L20T05[w/2]
# sensitive_ah$L20T10<-(sensitive_ah$L20T10-sensitive_ah$L20T10[w/2])/sensitive_ah$L20T10[w/2]
# sensitive_ah$L20T15<-(sensitive_ah$L20T15-sensitive_ah$L20T15[w/2])/sensitive_ah$L20T15[w/2]
# sensitive_ah$L20T20<-(sensitive_ah$L20T20-sensitive_ah$L20T20[w/2])/sensitive_ah$L20T20[w/2]


flag<-c("al")
sensitive(seq(0.1,10,0.1),flag=flag) %>%
  pivot_longer(L20T05:L20T20,names_to = "group") %>%
  ggplot()+
  geom_line(aes(x=x,y=value,color=group))+
  geom_hline(yintercept =0, color = "red") +
  labs(x=flag,y="change percent")+
  scale_x_continuous(limits = c(0,10),expand = c(0,0))+
  scale_y_continuous(limits = c(-1,1),expand = c(0,0))+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_text(size = 12,
                                   hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "right",
        plot.background = element_rect(colour = NA)) 

power<-seq(0.1,3,0.1)
flag<-c("bh")
choice<-c("relative")
sensitive_bh<-data.frame(power=power,
                         change=sensitive(power,flag=flag,value = choice)$L20T20,
                         group=rep(flag,length(power)))

flag<-c("ba")
sensitive_ba<-data.frame(power=power,
                         change=sensitive(power,flag=flag,value = choice)$L20T20,
                         group=rep(flag,length(power)))

flag<-c("bb")
sensitive_bb<-data.frame(power=power,
                         change=sensitive(power,flag=flag,value = choice)$L20T20,
                         group=rep(flag,length(power)))

flag<-c("bd")
sensitive_bd<-data.frame(power=power,
                         change=sensitive(power,flag=flag,value = choice)$L20T20,
                         group=rep(flag,length(power)))

flag<-c("bl")
sensitive_bl<-data.frame(power=power,
                         change=sensitive(power,flag=flag,value = choice)$L20T20,
                         group=rep(flag,length(power)))

sensitive_power<-rbind(sensitive_ba,sensitive_bb,sensitive_bd,sensitive_bl,sensitive_bh)

ggplot(sensitive_power)+
  geom_line(aes(x=power,y=change,color=group))+
  labs(x="power",y="change percent")+
  scale_x_continuous(limits = c(0,3),expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  mytheme()
  

coef<-seq(0,3.1,0.1)
flag<-c("V")
choice<-c("relative")
sensitive_bh<-data.frame(coef=coef,
                         change=sensitive(coef,flag=flag,value = choice)$L20T20,
                         group=rep(flag,length(coef)))

ggplot(sensitive_bh)+
  geom_line(aes(x=coef,y=change,color=group))+
  labs(x="coef",y="change percent")+
  scale_x_continuous(limits = c(0,3),expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  mytheme()

coef<-seq(25,60,1)
flag<-c("miu")
choice<-c("relative")
sensitive_aA<-data.frame(coef=coef,
                         change=sensitive(coef,flag=flag,value = choice)$L20T20,
                         group=rep(flag,length(coef)))

ggplot(sensitive_aA)+
  geom_line(aes(x=coef,y=c(change),color=group))+
  labs(x="coef",y="change percent")+
  scale_x_continuous(limits = c(25,60),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,1),expand = c(0,0))+
  geom_hline(yintercept = 0,linetype = "dashed") +
  mytheme()
ggsave(paste(flag,".pdf"),width=12,height=8,units="cm",dpi=600)


# parameter adjustment ----------------------------------------------------
#1
parameter<-data.frame(ah=0.05,bh=1,
             aA=0.02,bA=2,
             apump=10000000,bpump=-5,
             aalpha=1/1785.041,balpha=-3/800)
q<-cal_q1(L=temperature$length/10,T=temperature$T,
          ah=parameter$ah,bh=parameter$bh,
          aA=parameter$aA,bA=parameter$bA,
          apump =parameter$apump,bpump =parameter$bpump,
          aalpha = parameter$aalpha,balpha = parameter$balpha)

rmse<-sqrt(sum((q-temperature$filtration)^2)/length(q))
for (i in 1:9999){
  parameter1<-parameter
  # parameter1$ah<-parameter1$ah*(1+sample(c(-0.01,0.01),1))
  # parameter1$bh<-parameter1$bh*(1+sample(c(-0.01,0.01),1))
  # # parameter1$aA<-parameter1$aA*(1+sample(c(-0.01,0.01),1))
  # # parameter1$bA<-parameter1$bA*(1+sample(c(-0.01,0.01),1))
  # parameter1$apump<-parameter1$apump*(1+sample(c(-0.01,0.01),1))
  # parameter1$bpump<-parameter1$bpump*(1+sample(c(-0.01,0.01),1))
  parameter1$aalpha<-parameter1$aalpha*(1+sample(c(-0.01,0.01),1))
  parameter1$balpha<-parameter1$balpha*(1+sample(c(-0.01,0.01),1))
  
  q1<-cal_q1(L=temperature$length/10,T=temperature$T,
            ah=parameter1$ah,bh=parameter1$bh,
            aA=parameter1$aA,bA=parameter1$bA,
            apump =parameter1$apump,bpump =parameter1$bpump,
            aalpha = parameter1$aalpha,balpha = parameter1$balpha)
  rmse1<-sqrt(sum((q1-temperature$filtration)^2)/length(q1))
  
  if(rmse1<rmse){
    rmse<-rmse1
    parameter<-parameter1
  }
}
parametert<-parameter
plot(temperature$T,q1)
points(temperature$T,temperature$filtration,col="red")

plot(temperature$length,q11)
points(temperature$length,temperature$filtration,col="red")

predict.L<-rep(mean(temperature$length/10),30)
predict.T<-seq(1,30,1)
predict.q<-cal_q1(L=predict.L,T=predict.T,
                  ah=parameter$ah,bh=parameter$bh,
                  aA=parameter$aA,bA=parameter$bA,
                  apump =parameter$apump,bpump =parameter$bpump,
                  aalpha = parameter$aalpha,balpha = parameter$balpha)
predict.data<-data.frame(T=predict.T,filtration=predict.q)
p1<-myboxplot(x="tem",y="filtration",temperature)
p1<-p1+
  scale_x_continuous(limits = c(0,30),expand = c(0,0))+
  scale_y_continuous(limits = c(0,200),expand = c(0,0))+
  mytheme()
ggsave("tem-fil.pdf",width=12,height=8,units="cm",dpi=600)
p2<-ggplot(predict.data)+
  geom_line(aes(x=T,y=filtration))+
  scale_x_continuous(limits = c(0,30),expand = c(0,0))+
  scale_y_continuous(limits = c(0,200),expand = c(0,0))+
  mytheme()
ggsave("pre_tem-fil.pdf",width=12,height=8,units="cm",dpi=600)
#1.1
parameter<-data.frame(ah=0.05,bh=1,
                      aA=0.02,bA=2,
                      apump=1000000,bpump=-4,
                      aalpha=1/2000,balpha=-3/800)
q<-cal_q1(L=modeldata$L/10,T=modeldata$T,
          ah=parameter$ah,bh=parameter$bh,
          aA=parameter$aA,bA=parameter$bA,
          apump =parameter$apump,bpump =parameter$bpump,
          aalpha = parameter$aalpha,balpha = parameter$balpha)

rmse<-sqrt(sum((q-modeldata$q)^2)/length(q))
for (i in 1:9999){
  parameter1<-parameter
  parameter1$ah<-parameter1$ah*(1+sample(c(-0.01,0.01),1))
  parameter1$bh<-parameter1$bh*(1+sample(c(-0.01,0.01),1))
  # parameter1$aA<-parameter1$aA*(1+sample(c(-0.01,0.01),1))
  # parameter1$bA<-parameter1$bA*(1+sample(c(-0.01,0.01),1))
  parameter1$apump<-parameter1$apump*(1+sample(c(-0.01,0.01),1))
  parameter1$bpump<-parameter1$bpump*(1+sample(c(-0.01,0.01),1))
  parameter1$aalpha<-parameter1$aalpha*(1+sample(c(-0.01,0.01),1))
  parameter1$balpha<-parameter1$balpha*(1+sample(c(-0.01,0.01),1))
  
  q1<-cal_q1(L=modeldata$L/10,T=modeldata$T,
             ah=parameter1$ah,bh=parameter1$bh,
             aA=parameter1$aA,bA=parameter1$bA,
             apump =parameter1$apump,bpump =parameter1$bpump,
             aalpha = parameter1$aalpha,balpha = parameter1$balpha)
  rmse1<-sqrt(sum((q1-modeldata$q)^2)/length(q1))
  
  if(rmse1<rmse){
    rmse<-rmse1
    parameter<-parameter1
  }
}

plot(modeldata$T,q1)
points(modeldata$T,modeldata$q,col="red")

plot(modeldata$L,q1)
points(modeldata$L,modeldata$q,col="red")

#2
parameter<-data.frame(ah=0.05,bh=1,
                      aA=0.02,bA=2,
                      apump=10000000,bpump=-7,
                      aalpha=1/1785.041,balpha=-3/800)
q<-cal_q1(L=subset(t_60,deltat==60)$length/10,T=rep(20,length(subset(t_60,deltat==60)$length)),
          ah=parameter$ah,bh=parameter$bh,
          aA=parameter$aA,bA=parameter$bA,
          apump =parameter$apump,bpump =parameter$bpump,
          aalpha = parameter$aalpha,balpha = parameter$balpha)

rmse<-sqrt(sum((q-subset(t_60,deltat==60)$filtration)^2)/length(q))
rmses<-rep(0,999)
for (i in 1:999){
  parameter1<-parameter
  parameter1$ah<-parameter1$ah*(1+sample(c(-0.01,0.01),1))
  parameter1$bh<-parameter1$bh*(1+sample(c(-0.01,0.01),1))
  # parameter1$aA<-parameter1$aA*(1+sample(c(-0.01,0.01),1))
  # parameter1$bA<-parameter1$bA*(1+sample(c(-0.01,0.01),1))
  parameter1$apump<-parameter1$apump*(1+sample(c(-0.01,0.01),1))
  parameter1$bpump<-parameter1$bpump*(1+sample(c(-0.01,0.01),1))
  # parameter1$aalpha<-parameter1$aalpha*(1+sample(c(-0.01,0.01),1))
  # parameter1$balpha<-parameter1$balpha*(1+sample(c(-0.01,0.01),1))
  
  q1<-cal_q1(subset(t_60,deltat==60)$length/10,T=rep(20,length(subset(t_60,deltat==60)$length)),
             ah=parameter1$ah,bh=parameter1$bh,
             aA=parameter1$aA,bA=parameter1$bA,
             apump =parameter1$apump,bpump =parameter1$bpump,
             aalpha = parameter1$aalpha,balpha = parameter1$balpha)
  rmse1<-sqrt(sum((q1-subset(t_60,deltat==60)$filtration)^2)/length(q1))
  rmses[i]<-parameter$bh
  if(rmse1<rmse){
    rmse<-rmse1
    parameter<-parameter1
  }
}
parameterl1<-parameter
plot(subset(t_60,deltat==60)$length,q1)
points(subset(t_60,deltat==60)$length,subset(t_60,deltat==60)$filtration,col="red")

predict.L<-seq(0.1,2.5,0.1)
predict.T<-rep(20,25)
predict.q<-cal_q1(L=predict.L,T=predict.T,
                  ah=parameter$ah,bh=parameter$bh,
                  aA=parameter$aA,bA=parameter$bA,
                  apump =parameter$apump,bpump =parameter$bpump,
                  aalpha = parameter$aalpha,balpha = parameter$balpha)
predict.data<-data.frame(L=predict.L,filtration=predict.q)
p1<-myboxplot(x="length",y="filtration",subset(t_60,deltat==60),width=0.2)
p1<-p1+
  scale_x_continuous(limits = c(0,25),expand = c(0,0))+
  scale_y_continuous(limits = c(0,120),expand = c(0,0))+
  mytheme()
ggsave("l1-fil.pdf",width=12,height=8,units="cm",dpi=600)
p2<-ggplot(predict.data)+
  geom_line(aes(x=L*10,y=filtration))+
  scale_x_continuous(limits = c(0,25),expand = c(0,0))+
  scale_y_continuous(limits = c(0,120),expand = c(0,0))+
  mytheme()
ggsave("pre_l1-fil.pdf",width=12,height=8,units="cm",dpi=600)



q<-cal_q1(L=subset(t1_60,deltat==60)$length/10,T=rep(20,length(subset(t1_60,deltat==60)$length)),
          ah=parameter$ah,bh=parameter$bh,
          aA=parameter$aA,bA=parameter$bA,
          apump =parameter$apump,bpump =parameter$bpump,
          aalpha = parameter$aalpha,balpha = parameter$balpha)

rmse<-sqrt(sum((q-subset(t1_60,deltat==60)$filtration)^2)/length(q))
rmses<-rep(0,999)
for (i in 1:999){
  parameter1<-parameter
  parameter1$ah<-parameter1$ah*(1+sample(c(-0.01,0.01),1))
  parameter1$bh<-parameter1$bh*(1+sample(c(-0.01,0.01),1))
  # parameter1$aA<-parameter1$aA*(1+sample(c(-0.01,0.01),1))
  # parameter1$bA<-parameter1$bA*(1+sample(c(-0.01,0.01),1))
  parameter1$apump<-parameter1$apump*(1+sample(c(-0.01,0.01),1))
  parameter1$bpump<-parameter1$bpump*(1+sample(c(-0.01,0.01),1))
  parameter1$aalpha<-parameter1$aalpha*(1+sample(c(-0.01,0.01),1))
  parameter1$balpha<-parameter1$balpha*(1+sample(c(-0.01,0.01),1))
  
  q1<-cal_q1(subset(t1_60,deltat==60)$length/10,T=rep(20,length(subset(t1_60,deltat==60)$length)),
             ah=parameter1$ah,bh=parameter1$bh,
             aA=parameter1$aA,bA=parameter1$bA,
             apump =parameter1$apump,bpump =parameter1$bpump,
             aalpha = parameter1$aalpha,balpha = parameter1$balpha)
  rmse1<-sqrt(sum((q1-subset(t1_60,deltat==60)$filtration)^2)/length(q1))
  rmses[i]<-parameter$bh
  if(rmse1<rmse){
    rmse<-rmse1
    parameter<-parameter1
  }
}
parameterl2<-parameter
plot(subset(t1_60,deltat==60)$length,q1)
points(subset(t1_60,deltat==60)$length,subset(t1_60,deltat==60)$filtration,col="red")

predict.L<-seq(0.1,2.5,0.1)
predict.T<-rep(20,25)
predict.q<-cal_q1(L=predict.L,T=predict.T,
                  ah=parameter$ah,bh=parameter$bh,
                  aA=parameter$aA,bA=parameter$bA,
                  apump =parameter$apump,bpump =parameter$bpump,
                  aalpha = parameter$aalpha,balpha = parameter$balpha)
predict.data<-data.frame(L=predict.L,filtration=predict.q)
p1<-myboxplot(x="length",y="filtration",subset(t1_60,deltat==60),width=0.2)
p1<-p1+
  scale_x_continuous(limits = c(0,25),expand = c(0,0))+
  scale_y_continuous(limits = c(0,190),expand = c(0,0))+
  mytheme()
ggsave("l2-fil.pdf",width=12,height=8,units="cm",dpi=600)
p2<-ggplot(predict.data)+
  geom_line(aes(x=L*10,y=filtration))+
  scale_x_continuous(limits = c(0,25),expand = c(0,0))+
  scale_y_continuous(limits = c(0,190),expand = c(0,0))+
  mytheme()
ggsave("pre_l2-fil.pdf",width=12,height=8,units="cm",dpi=600)
#3
parameter<-data.frame(ah=0.05,bh=0.5,
                      aA=0.02,bA=1,
                      apump=10000000,bpump=-3,
                      aalpha=1/1785.041,balpha=-3/800)
q<-cal_q1(L=subset(t_60,deltat==60)$length/10,T=rep(20,length(subset(t_60,deltat==60)$length)),
          ah=parameter$ah,bh=parameter$bh,
          aA=parameter$aA,bA=parameter$bA,
          apump =parameter$apump,bpump =parameter$bpump,
          aalpha = parameter$aalpha,balpha = parameter$balpha)

rmse0<-sqrt(sum((q-subset(t_60,deltat==60)$filtration)^2)/length(q))
q1<-cal_q1(L=temperature$length/10,T=temperature$T,
          ah=parameter$ah,bh=parameter$bh,
          aA=parameter$aA,bA=parameter$bA,
          apump =parameter$apump,bpump =parameter$bpump,
          aalpha = parameter$aalpha,balpha = parameter$balpha)

rmse1<-sqrt(sum((q1-temperature$filtration)^2)/length(q1))
rmse<-rmse0+rmse1
for (i in 1:9999){
  parameter1<-parameter
  parameter1$ah<-parameter1$ah*(1+sample(c(-0.01,0.01),1))
  parameter1$bh<-parameter1$bh*(1+sample(c(-0.01,0.01),1))
  # parameter1$aA<-parameter1$aA*(1+sample(c(-0.01,0.01),1))
  # parameter1$bA<-parameter1$bA*(1+sample(c(-0.01,0.01),1))
  parameter1$apump<-parameter1$apump*(1+sample(c(-0.01,0.01),1))
  parameter1$bpump<-parameter1$bpump*(1+sample(c(-0.01,0.01),1))
  parameter1$aalpha<-parameter1$aalpha*(1+sample(c(-0.01,0.01),1))
  parameter1$balpha<-parameter1$balpha*(1+sample(c(-0.01,0.01),1))
  
  q01<-cal_q1(subset(t_60,deltat==60)$length/10,T=rep(20,length(subset(t_60,deltat==60)$length)),
             ah=parameter1$ah,bh=parameter1$bh,
             aA=parameter1$aA,bA=parameter1$bA,
             apump =parameter1$apump,bpump =parameter1$bpump,
             aalpha = parameter1$aalpha,balpha = parameter1$balpha)
  rmse01<-sqrt(sum((q01-subset(t_60,deltat==60)$filtration)^2)/length(q01))
  q11<-cal_q1(L=temperature$length/10,T=temperature$T,
            ah=parameter$ah,bh=parameter$bh,
            aA=parameter$aA,bA=parameter$bA,
            apump =parameter$apump,bpump =parameter$bpump,
            aalpha = parameter$aalpha,balpha = parameter$balpha)
  
  rmse11<-sqrt(sum((q11-temperature$filtration)^2)/length(q11))
  rmse00<-rmse01+rmse11
  if(rmse00<rmse){
    rmse<-rmse00
    parameter<-parameter1
  }
}
parameterl1t<-parameter
plot(subset(t_60,deltat==60)$length,q01)
points(subset(t_60,deltat==60)$length,subset(t_60,deltat==60)$filtration,col="red")


temperature %>%
  split(.$group) %>%
  sapply(function(x) c(T=mean(x$T),L=mean(x$length),q=mean(x$filtration))) %>%
  t() %>%
  data.frame()->modeldata
###per unit weight
parameter<-data.frame(ah=0.05,bh=1,
                      aA=0.02,bA=2,
                      apump=10000000,bpump=-5,
                      aalpha=1/1785.041,balpha=-3/800)
q<-cal_q1(L=temperature$length/10,T=temperature$T,
          ah=parameter$ah,bh=parameter$bh,
          aA=parameter$aA,bA=parameter$bA,
          apump =parameter$apump,bpump =parameter$bpump,
          aalpha = parameter$aalpha,balpha = parameter$balpha)

rmse<-sqrt(sum((q-temperature$filtration)^2)/length(q))
for (i in 1:9999){
  parameter1<-parameter
  # parameter1$ah<-parameter1$ah*(1+sample(c(-0.01,0.01),1))
  # parameter1$bh<-parameter1$bh*(1+sample(c(-0.01,0.01),1))
  # # parameter1$aA<-parameter1$aA*(1+sample(c(-0.01,0.01),1))
  # # parameter1$bA<-parameter1$bA*(1+sample(c(-0.01,0.01),1))
  # parameter1$apump<-parameter1$apump*(1+sample(c(-0.01,0.01),1))
  # parameter1$bpump<-parameter1$bpump*(1+sample(c(-0.01,0.01),1))
  parameter1$aalpha<-parameter1$aalpha*(1+sample(c(-0.01,0.01),1))
  parameter1$balpha<-parameter1$balpha*(1+sample(c(-0.01,0.01),1))
  
  q1<-cal_q1(L=temperature$length/10,T=temperature$T,
             ah=parameter1$ah,bh=parameter1$bh,
             aA=parameter1$aA,bA=parameter1$bA,
             apump =parameter1$apump,bpump =parameter1$bpump,
             aalpha = parameter1$aalpha,balpha = parameter1$balpha)
  rmse1<-sqrt(sum((q1-temperature$filtration)^2)/length(q1))
  
  if(rmse1<rmse){
    rmse<-rmse1
    parameter<-parameter1
  }
}
parametert<-parameter
plot(temperature$T,q1)
points(temperature$T,temperature$filtration,col="red")

plot(temperature$length,q11)
points(temperature$length,temperature$filtration,col="red")

predict.L<-rep(mean(temperature$length/10),30)
predict.T<-seq(1,30,1)
predict.q<-cal_q1(L=predict.L,T=predict.T,
                  ah=parameter$ah,bh=parameter$bh,
                  aA=parameter$aA,bA=parameter$bA,
                  apump =parameter$apump,bpump =parameter$bpump,
                  aalpha = parameter$aalpha,balpha = parameter$balpha)
predict.data<-data.frame(T=predict.T,filtration=predict.q)
p1<-myboxplot(x="tem",y="filtration",temperature)
p1<-p1+
  scale_x_continuous(limits = c(0,30),expand = c(0,0))+
  scale_y_continuous(limits = c(0,200),expand = c(0,0))+
  mytheme()
ggsave("tem-fil.pdf",width=12,height=8,units="cm",dpi=600)
p2<-ggplot(predict.data)+
  geom_line(aes(x=T,y=filtration))+
  scale_x_continuous(limits = c(0,30),expand = c(0,0))+
  scale_y_continuous(limits = c(0,200),expand = c(0,0))+
  mytheme()
ggsave("pre_tem-fil.pdf",width=12,height=8,units="cm",dpi=600)

#文献汇总的数据模拟

parameter<-data.frame(ah=0.05,bh=1,
                      aA=0.02,bA=2,
                      apump=10000000,bpump=-5,
                      aalpha=1/1785.041,balpha=-3/800)
q<-cal_q1(L=subset(lit_fil1,reference=="sylvester et al. 2006")$length/10,
          T=subset(lit_fil1,reference=="sylvester et al. 2006")$temperature,
          ah=parameter$ah,bh=parameter$bh,
          aA=parameter$aA,bA=parameter$bA,
          apump =parameter$apump,bpump =parameter$bpump,
          aalpha = parameter$aalpha,balpha = parameter$balpha)

rmse<-sqrt(sum((q-subset(lit_fil1,reference=="sylvester et al. 2006")$filtration)^2)/length(q))
for (i in 1:9999){
  parameter1<-parameter
  parameter1$ah<-parameter1$ah*(1+sample(c(-0.01,0.01),1))
  parameter1$bh<-parameter1$bh*(1+sample(c(-0.01,0.01),1))
  parameter1$aA<-parameter1$aA*(1+sample(c(-0.01,0.01),1))
  parameter1$bA<-parameter1$bA*(1+sample(c(-0.01,0.01),1))
  parameter1$apump<-parameter1$apump*(1+sample(c(-0.01,0.01),1))
  parameter1$bpump<-parameter1$bpump*(1+sample(c(-0.01,0.01),1))
  parameter1$aalpha<-parameter1$aalpha*(1+sample(c(-0.01,0.01),1))
  parameter1$balpha<-parameter1$balpha*(1+sample(c(-0.01,0.01),1))
  
  q1<-cal_q1(L=subset(lit_fil1,reference=="sylvester et al. 2006")$length/10,
             T=subset(lit_fil1,reference=="sylvester et al. 2006")$temperature,
             ah=parameter1$ah,bh=parameter1$bh,
             aA=parameter1$aA,bA=parameter1$bA,
             apump =parameter1$apump,bpump =parameter1$bpump,
             aalpha = parameter1$aalpha,balpha = parameter1$balpha)
  rmse1<-sqrt(sum((q1-subset(lit_fil1,reference=="sylvester et al. 2006")$filtration)^2)/length(q1))
  
  if(rmse1<rmse){
    rmse<-rmse1
    parameter<-parameter1
  }
}
# plot(subset(lit_fil1,reference="Tokumon et al. 2015")$temperature,q1)
# points(subset(lit_fil1,reference=="Tokumon et al. 2015")$temperature,subset(lit_fil1,reference=="Tokumon et al. 2015")$filtration,col="red")
# plot(subset(lit_fil1,reference=="Tokumon et al. 2015")$filtration,q1)
# lines(x=c(0,150),y=c(0,150))
# plot(subset(lit_fil1,reference==reference)$length,q11)
# points(subset(lit_fil1,reference==reference)$length,subset(lit_fil1,reference==reference)$filtration,col="red")
predict.L<-rep(mean(temperature$length/10),30)
predict.T<-seq(1,30,1)
predict.q<-cal_q1(L=predict.L,T=predict.T,
                  ah=parameter$ah,bh=parameter$bh,
                  aA=parameter$aA,bA=parameter$bA,
                  apump =parameter$apump,bpump =parameter$bpump,
                  aalpha = parameter$aalpha,balpha = parameter$balpha)
predict.data<-data.frame(T=predict.T,filtration=predict.q)
p1<-myboxplot(x="temperature",y="filtration",subset(lit_fil1,reference=="Pestana, et al, 2009"),width = 0.05)
p1<-p1+
  scale_x_continuous(limits = c(0,25),expand = c(0,0))+
  scale_y_continuous(limits = c(0,20),expand = c(0,0))+
  mytheme()
ggsave("tem-fil.pdf",width=12,height=8,units="cm",dpi=600)
p2<-ggplot(predict.data)+
  geom_line(aes(x=T,y=filtration))+
  scale_x_continuous(limits = c(0,30),expand = c(0,0))+
  scale_y_continuous(limits = c(0,200),expand = c(0,0))+
  mytheme()
ggsave("pre_tem-fil.pdf",width=12,height=8,units="cm",dpi=600)


# eduils ------------------------------------------------------------------

mefil<-read.csv("D:/Seafile/私人资料库/fil/mytilus/edulis/fil.csv")
parameter<-data.frame(ah=0.05,bh=1,
                      aA=0.02,bA=2,
                      apump=10000000,bpump=-5,
                      aalpha=1/1785.041,balpha=-3/800)
q<-cal_q1(L=subset(mefil,reference=="Luskow et al 2018")$length,
          T=subset(mefil,reference=="Luskow et al 2018")$temperature,
          ah=parameter$ah,bh=parameter$bh,
          aA=parameter$aA,bA=parameter$bA,
          apump =parameter$apump,bpump =parameter$bpump,
          aalpha = parameter$aalpha,balpha = parameter$balpha)

rmse<-sqrt(sum((q-subset(mefil,reference=="Luskow et al 2018")$filtration)^2)/length(q))
for (i in 1:9999){
  parameter1<-parameter
  parameter1$ah<-parameter1$ah*(1+sample(c(-0.01,0.01),1))
  parameter1$bh<-parameter1$bh*(1+sample(c(-0.01,0.01),1))
  parameter1$aA<-parameter1$aA*(1+sample(c(-0.01,0.01),1))
  parameter1$bA<-parameter1$bA*(1+sample(c(-0.01,0.01),1))
  parameter1$apump<-parameter1$apump*(1+sample(c(-0.01,0.01),1))
  parameter1$bpump<-parameter1$bpump*(1+sample(c(-0.01,0.01),1))
  parameter1$aalpha<-parameter1$aalpha*(1+sample(c(-0.01,0.01),1))
  parameter1$balpha<-parameter1$balpha*(1+sample(c(-0.01,0.01),1))
  
  q1<-cal_q1(L=subset(mefil,reference=="Luskow et al 2018")$length,
             T=subset(mefil,reference=="Luskow et al 2018")$temperature,
             ah=parameter1$ah,bh=parameter1$bh,
             aA=parameter1$aA,bA=parameter1$bA,
             apump =parameter1$apump,bpump =parameter1$bpump,
             aalpha = parameter1$aalpha,balpha = parameter1$balpha)
  rmse1<-sqrt(sum((q1-subset(mefil,reference=="Luskow et al 2018")$filtration)^2)/length(q1))
  
  if(rmse1<rmse){
    rmse<-rmse1
    parameter<-parameter1
  }
}
rmse
parameter
p1<-myboxplot(x="length",y="filtration",subset(mefil,reference=="riisgard et al 2014"),width = 0.05)
p1<-p1+
  scale_x_continuous(limits = c(0,25),expand = c(0,0))+
  scale_y_continuous(limits = c(0,20),expand = c(0,0))+
  mytheme()
ggsave("tem-fil.pdf",width=12,height=8,units="cm",dpi=600)
p2<-ggplot(predict.data)+
  geom_line(aes(x=T,y=filtration))+
  scale_x_continuous(limits = c(0,30),expand = c(0,0))+
  scale_y_continuous(limits = c(0,200),expand = c(0,0))+
  mytheme()

predict.L<-rep(mean(temperature$length/10),30)
predict.T<-seq(1,30,1)
predict.q<-cal_q1(L=predict.L,T=predict.T,
                  ah=parameter$ah,bh=parameter$bh,
                  aA=parameter$aA,bA=parameter$bA,
                  apump =parameter$apump,bpump =parameter$bpump,
                  aalpha = parameter$aalpha,balpha = parameter$balpha)
predict.data<-data.frame(T=predict.T,filtration=predict.q)

t1_model<-nls(filtration~sqrt(a*length^(b)+c*length^(d)),
              data=subset(mefil,reference=="riisgard et al 2014"),
              start = list(a=1000000,b=4,c=0,d=4))
fil1_predict<-predict(t1_model)
RSS<-sum((subset(t1_60,deltat==60)$fil1-fil1_predict)^2)
TSS<-sum((subset(t1_60,deltat==60)$fil1-mean(subset(t1_60,deltat==60)$fil1))^2)
1-RSS/TSS


#ui
q2<-cal_q2(L,T,
           ah=best.parameter$ah,bh=best.parameter$bh,
           aA=best.parameter$aA,bA=best.parameter$bA,
           apump =best.parameter$apump,bpump =best.parameter$bpump,
           aalpha = best.parameter$aalpha,balpha = best.parameter$balpha)
parameter<-parametert


q2<-cal_q2(L,T,
           ah=parameter$ah,bh=parameter$bh,
           aA=parameter$aA,bA=parameter$bA,
           apump =parameter$apump,bpump =parameter$bpump,
           aalpha = parameter$aalpha,balpha = parameter$balpha)

jet.colors <- colorRamp2(c(min(q2), max(q2)), c("#4b5cc4","#f47983"))
# Generate the desired number of colors from this palette
nbcol <- 100
color <- jet.colors(nbcol)
facetcol <- cut(q2, nbcol)

persp(L,T,q2, theta = 315, phi = 20,
      expand = 0.5,
      r=180,
      ltheta = 225,
      shade = 0,
      col = c("#f47983"),
      ticktype = "detailed")
library(shiny)

ui <- fluidPage(
  titlePanel("Hello Shiny!"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("ah", "ah:", min = 0, max = 10, value = 0.5),
      sliderInput("bh", "bh:", min = 0, max = 4, value = 2),
      sliderInput("aA", "aA:", min = 0, max = 4, value = 2),
      sliderInput("bA", "bA:", min = 0, max = 4, value = 2),
      sliderInput("apump", "apump:", min = 0, max = 1000000, value = 50),
      sliderInput("bpump", "bpump:", min = -10, max = 0, value = -5),
      sliderInput("aalpha", "aalpha:", min = 0, max = 1, value = 0.5),
      sliderInput("balpha", "balpha:", min = -1, max = 0, value = -0.5),
      sliderInput("theta", "theta:", min = 0, max = 360, value = 225),
      sliderInput("phi", "phi:", min = 0, max = 90, value = 20),
    ),
    mainPanel(
      plotOutput("myplot")
      )
  )
   
)
server<-function(input,output){
  output$myplot<-renderPlot({
    L<-seq(0.1,3,0.1)
    T<-seq(1,30,1)
    q2<-cal_q2(L,T,
               ah=input$ah,bh=input$bh,
               aA=input$aA,bA=input$bA,
               apump =input$apump,bpump =input$bpump,
               aalpha = input$aalpha,balpha = input$balpha)
    persp(L,T,q2, theta = input$theta, phi = input$phi,
          expand = 0.5,
          r=180,
          ltheta = 225,
          shade = 0.75,
          ticktype = "detailed")
    
    
  })
}
shinyApp(ui = ui, server = server)

ui <- fluidPage(
  titlePanel("Hello Shiny!"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("cc", "c:", min = 0, max = 1000, value = 50),
      sliderInput("dd", "d:", min = 0, max = 1000, value = 50),
      sliderInput("vv", "v:", min = 0, max = 100, value = 2)
    ),
    mainPanel(
      plotOutput("myplot")
    )
  )
  
)
server<-function(input,output){
  output$myplot<-renderPlot({
    a<-seq(0.01,1,0.01)
    q2<-a*sqrt(input$cc*input$vv*a*a+input$dd*a)
    data=data.frame(a=a,q2=q2)
    ggplot(data)+
      geom_line(aes(x=a,y=q2))+
      geom_line(aes(x=a,y=max(q2)*a))
  })
}
shinyApp(ui = ui, server = server)
  
#目前最好的
#        ah        bh   aA bA    apump     bpump     aalpha       balpha
# 1.066166 0.7405881 0.02  2 713335.4 -1.446426 0.01124989 0.01124989
best.parameter<-data.frame(ah=1.066166,bh=0.7405881,
                           aA=0.02,bA=2,
                           apump=713335.4,bpump=-1.446426,
                           aalpha=0.01124989,balpha=-0.01124989)


# effect of the temperature on mussel pump model --------------------------
digitize("D:/Seafile/私人资料库/fil/1.png")

fil.temp1<-read.csv("D:/Seafile/私人资料库/fil/1.csv")
fil.temp1$miu<-(-0.0322*fil.temp1$temperature+1.6927)*6 #unit 10-3 cm2/min

control<-data.frame(temperature=c(0,0,0),
                    filtration=c(0.1,0.1,0.1),
                    group=c("a","b","c"),
                    miu=c(10.1562,10.1562,10.1562))

# miu simplified model---------------------------------------------------------------------

rbind(fil.temp1,control) %>%
  subset(group=="a") %>%
  nls(miu~a*filtration+b/filtration+c,data=.,
      start = list(a=-0.2,b=1,c=1)) -> fil.temp1.modela
rbind(fil.temp1,control) %>%
  subset(group=="b") %>%
  nls(miu~a*filtration+b/filtration+c,data=.,
      start = list(a=-0.2,b=1,c=1)) -> fil.temp1.modelb
rbind(fil.temp1,control) %>%
  subset(group=="c") %>%
  nls(miu~a*filtration+b/filtration+c,data=.,
      start = list(a=-0.2,b=1,c=1)) -> fil.temp1.modelc

fil.temp1.modela<-nls(miu~a*filtration+b/filtration+C,data=subset(fil.temp1,group=="a"),
                      start = list(a=-0.2,b=1,C=1))
fil.temp1.modelb<-nls(miu~a*filtration+b/filtration+C,data=subset(fil.temp1,group=="b"),
                      start = list(a=-0.2,b=1,C=1))
fil.temp1.modelc<-nls(miu~a*filtration+b/filtration+C,data=subset(fil.temp1,group=="c"),
                      start = list(a=-0.2,b=1,C=1))

fil.temp1.modelal<- lm(miu~filtration,data=subset(fil.temp1,group=="a"))
fil.temp1.modelbl<- lm(miu~filtration,data=subset(fil.temp1,group=="b"))
fil.temp1.modelcl<- lm(miu~filtration,data=subset(fil.temp1,group=="c"))

predict.filtration<-data.frame(filtration=seq(1,60,1))
predict.miua<-predict(fil.temp1.modela,newdata=predict.filtration)
predict.miub<-predict(fil.temp1.modelb,newdata=predict.filtration)
predict.miuc<-predict(fil.temp1.modelc,newdata=predict.filtration)

predict.miual<-predict(fil.temp1.modelal,newdata=predict.filtration)
predict.miubl<-predict(fil.temp1.modelbl,newdata=predict.filtration)
predict.miucl<-predict(fil.temp1.modelcl,newdata=predict.filtration)

predict.a<-data.frame(miu=predict.miua,
                      filtration=predict.filtration,
                      group="I",
                      method="simplified")
predict.b<-data.frame(miu=predict.miub,
                      filtration=predict.filtration,
                      group="II",
                      method="simplified")
predict.c<-data.frame(miu=predict.miuc,
                      filtration=predict.filtration,
                      group="III",
                      method="simplified")

predict.al<-data.frame(miu=predict.miual,
                      filtration=predict.filtration,
                      group="I",
                      method="lm")
predict.bl<-data.frame(miu=predict.miubl,
                       filtration=predict.filtration,
                       group="II",
                       method="lm")
predict.cl<-data.frame(miu=predict.miucl,
                       filtration=predict.filtration,
                       group="III",
                       method="lm")
orginal.data<-data.frame(miu=fil.temp1$miu,
                         filtration=fil.temp1$filtration,
                         group=fil.temp1$group,
                         method="org")
predict.model<-rbind(predict.al,predict.bl,predict.cl,predict.a,predict.b,predict.c)


ggplot()+
  geom_point(data=orginal.data,aes(y=miu,x=filtration,color=group))+
  geom_line(data=predict.model,aes(y=miu,x=filtration,color=group,linetype=method))+
  scale_x_continuous(limits = c(-10,100))

predict.miua<-predict(fil.temp1.modela,newdata=subset(fil.temp1,group=="a")$filtration)
ssr<-sum((subset(fil.temp1,group=="a")$miu-predict.miua)^2)
sst<-3*var(subset(fil.temp1,group=="a")$miu)
1-ssr/sst

# T simplified model ------------------------------------------------------

rbind(fil.temp1,control) %>%
  subset(group=="a") %>%
  nls(temperature~a*filtration+b/filtration+c,data=.,
      start = list(a=0.5,b=1.5,c=-15)) -> fil.temp1.modela
rbind(fil.temp1,control) %>%
  subset(group=="b") %>%
  nls(temperature~a*filtration+b/filtration+c,data=.,
      start = list(a=0.5,b=1.5,c=-15)) -> fil.temp1.modelb
rbind(fil.temp1,control) %>%
  subset(group=="c") %>%
  nls(temperature~a*filtration+b/filtration+c,data=.,
      start = list(a=0.5,b=1.5,c=-15)) -> fil.temp1.modelc

fil.temp1.modela<-nls(temperature~a*filtration+b/filtration+C,data=subset(fil.temp1,group=="a"),
                      start = list(a=-0.2,b=1,C=1))
fil.temp1.modelb<-nls(temperature~a*filtration+b/filtration+C,data=subset(fil.temp1,group=="b"),
                      start = list(a=-0.2,b=1,C=1))
fil.temp1.modelc<-nls(temperature~a*filtration+b/filtration+C,data=subset(fil.temp1,group=="c"),
                      start = list(a=-0.2,b=1,C=1))

fil.temp1.modelal<- lm(temperature~filtration,data=subset(fil.temp1,group=="a"))
fil.temp1.modelbl<- lm(temperature~filtration,data=subset(fil.temp1,group=="b"))
fil.temp1.modelcl<- lm(temperature~filtration,data=subset(fil.temp1,group=="c"))

predict.filtration<-data.frame(filtration=seq(10,60,1))
predict.temperaturea<-predict(fil.temp1.modela,newdata=predict.filtration)
predict.temperatureb<-predict(fil.temp1.modelb,newdata=predict.filtration)
predict.temperaturec<-predict(fil.temp1.modelc,newdata=predict.filtration)

predict.temperatureal<-predict(fil.temp1.modelal,newdata=predict.filtration)
predict.temperaturebl<-predict(fil.temp1.modelbl,newdata=predict.filtration)
predict.temperaturecl<-predict(fil.temp1.modelcl,newdata=predict.filtration)

predict.a<-data.frame(temperature=predict.temperaturea,
                      filtration=predict.filtration,
                      group="I",
                      method="simplified")
predict.b<-data.frame(temperature=predict.temperatureb,
                      filtration=predict.filtration,
                      group="II",
                      method="simplified")
predict.c<-data.frame(temperature=predict.temperaturec,
                      filtration=predict.filtration,
                      group="III",
                      method="simplified")

predict.al<-data.frame(temperature=predict.temperatureal,
                       filtration=predict.filtration,
                       group="I",
                       method="lm")
predict.bl<-data.frame(temperature=predict.temperaturebl,
                       filtration=predict.filtration,
                       group="II",
                       method="lm")
predict.cl<-data.frame(temperature=predict.temperaturecl,
                       filtration=predict.filtration,
                       group="III",
                       method="lm")
orginal.data<-data.frame(temperature=fil.temp1$temperature,
                         filtration=fil.temp1$filtration,
                         group=fil.temp1$group,
                         method="org")
predict.model<-rbind(predict.al,predict.bl,predict.cl,predict.a,predict.b,predict.c)


ggplot()+
  geom_point(data=orginal.data,aes(y=temperature,x=filtration,color=group))+
  geom_line(data=predict.model,aes(y=temperature,x=filtration,color=group,linetype=method))+
  scale_y_continuous(limits = c(0,24))+
  mytheme()
ggsave("simplfied model.pdf",width=12,height=8,units="cm",dpi=600)

predict.temperaturec<-predict(fil.temp1.modelc,newdata=subset(fil.temp1,group=="c")$filtration)
ssr<-sum((subset(fil.temp1,group=="c")$temperature-predict.temperaturec)^2)
sst<-24*var(subset(fil.temp1,group=="c")$temperature)
1-ssr/sst
##kittnner
fil.temp1.modeld<-nls(temperature~a*filtration+b/filtration+C,data=subset(fil.temp1,group=="d"),
                      start = list(a=-0.2,b=1,C=1))
fil.temp1.modele<-nls(temperature~a*filtration+b/filtration+C,data=subset(fil.temp1,group=="e"),
                      start = list(a=-0.2,b=1,C=1))

fil.temp1.modeldl<- lm(temperature~filtration,data=subset(fil.temp1,group=="d"))
fil.temp1.modelel<- lm(temperature~filtration,data=subset(fil.temp1,group=="e"))

predict.filtration<-data.frame(filtration=seq(1,120,1))
predict.temperatured<-predict(fil.temp1.modeld,newdata=predict.filtration)
predict.temperaturee<-predict(fil.temp1.modele,newdata=predict.filtration)

predict.temperaturedl<-predict(fil.temp1.modeldl,newdata=predict.filtration)
predict.temperatureel<-predict(fil.temp1.modelel,newdata=predict.filtration)

predict.d<-data.frame(temperature=predict.temperatured,
                      filtration=predict.filtration,
                      group="IV",
                      method="simplified")
predict.e<-data.frame(temperature=predict.temperaturee,
                      filtration=predict.filtration,
                      group="V",
                      method="simplified")
predict.dl<-data.frame(temperature=predict.temperaturedl,
                       filtration=predict.filtration,
                       group="IV",
                       method="lm")
predict.el<-data.frame(temperature=predict.temperatureel,
                       filtration=predict.filtration,
                       group="V",
                       method="lm")
predict.model<-rbind(predict.dl,predict.el,predict.d,predict.e)

ggplot()+
  geom_point(data=subset(fil.temp1,group=="d"|group=="e"),
             aes(y=temperature,x=filtration,color=group))+
  geom_line(data=predict.model,aes(y=temperature,x=filtration,color=group,linetype=method))+
  scale_y_continuous(limits = c(5,24))+
  mytheme()
ggsave("simplfied model2.pdf",width=12,height=8,units="cm",dpi=600)

predict.temperaturee<-predict(fil.temp1.modele,newdata=subset(fil.temp1,group=="e")$filtration)
ssr<-sum((subset(fil.temp1,group=="e")$temperature-predict.temperaturee)^2)
sst<-16*var(subset(fil.temp1,group=="e")$temperature)
1-ssr/sst



# L simplified model ------------------------------------------------------
digitize("D:/Seafile/私人资料库/fil/2.png")
fil.length<-read.csv("D:/Seafile/私人资料库/fil/2.csv")

fil.length.modelall<-nls(filtration~a*length^b,data=fil.length,
                      start = list(a=0.00079,b=2.508))


predict.length<-data.frame(length=seq(0.5,35,0.5))
predict.filtrationall<-predict(fil.length.modelall,newdata=predict.length)
predict.filtration<-0.00010277*predict.length$length^4.3207*(sqrt(1+162580*predict.length$length^(-3.1059))-1)

predict.all<-data.frame(length=predict.length,
                      filtration=predict.filtrationall,
                      group="I",
                      method="all")
predict.sim<-data.frame(length=predict.length,
                        filtration=predict.filtration,
                        group="I",
                        method="simplified")

predict.lmodel<-rbind(predict.all,predict.sim)

ggplot()+
  geom_point(data=fil.length,
             aes(x=length,y=filtration))+
  geom_line(data=predict.lmodel,aes(x=length,y=filtration,linetype=method))+
  scale_x_continuous(limits = c(0,30))+
  mytheme()
ggsave("simplfied modell.pdf",width=12,height=8,units="cm",dpi=600)

predict.filtration<-predict(fil.length.modelall,newdata=fil.length$length)
predict.filtration<-0.00010277*fil.length$length^4.3207*(sqrt(1+162580*fil.length$length^(-3.1059))-1)
ssr<-sum((fil.length$filtration-predict.filtration)^2)
sst<-25*var(fil.length$filtration)
1-ssr/sst

# model = inline('a(1)*x.^a(2).*((1+a(3)*x.^a(4)).^0.5-1)','a','x');
# a0=[0.00079,2.508,173000,-1];
# a=nlinfit(x,y,model,a0)


ui <- fluidPage(
  titlePanel("Hello Shiny!"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("cc", "c:", min = -1, max = 1, value = 1),
      sliderInput("dd", "d:", min = 0, max = 1000, value = 50)
          ),
    mainPanel(
      plotOutput("myplot")
    )
  )
  
)
server<-function(input,output){
  output$myplot<-renderPlot({
    a<-seq(1,60,1)
    
    q2<-input$cc*(a)^(-1)+input$vv*a
    data=data.frame(a=a,q2=q2)
    ggplot(data)+
      geom_line(aes(x=a,y=q2))
      
  })
}
shinyApp(ui = ui, server = server)
