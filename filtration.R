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
# velocity -----------------------------------------------------------------------

velocity<-read.csv("lentic.csv")
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