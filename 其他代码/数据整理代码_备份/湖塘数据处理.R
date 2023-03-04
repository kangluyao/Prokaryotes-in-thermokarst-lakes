pond_table<-read.csv('E:/湖塘/数据处理/ponddata.csv',header=TRUE,stringsAsFactors = FALSE)
library(MASS)
doc_dat<-data.frame(t(sapply(pond_table[c('c_TOC_T0','c_TOC_T28')],function(x)
  (c(mean=mean(x),SD=sd(x))))))
se=doc_dat$SD/sqrt(102)
with(pond_table,t.test(c_TOC_T0,c_TOC_T28,paired=T))
sig<-c('a','b')
doc_dat<-data.frame(Time=rownames(doc_dat),doc_dat,se,sig,stringsAsFactors=F)

p1<-ggplot(pond_dat,aes(x=Time,y=mean))+
  geom_bar(stat = 'identity',fill=c('darkgreen','pink3'),color='black',width = 0.45)+
  geom_errorbar(aes(ymin=mean,ymax=mean+se),width=.2)+
  geom_text(aes(label=sig,y=(mean+se)*1.02),position = position_dodge(0.9),vjust = 0)+
  scale_x_discrete(limits=c('c_TOC_T0','c_TOC_T28'))+
  scale_y_continuous(expand = c(0,0),limits = c(0,18))+
  xlab('Time')+ylab('DOC (mg/L)')+
  theme_classic()

TN_dat<-data.frame(t(sapply(pond_table[c('c_TN_T0','c_TN_T28')],function(x)
  (c(mean=mean(x),SD=sd(x))))))
se=TN_dat$SD/sqrt(102)
with(pond_table,t.test(c_TN_T0,c_TN_T28,paired=T))
sig<-c('a','b')
TN_dat<-data.frame(Time=rownames(TN_dat),TN_dat,se,sig,stringsAsFactors=F)


p2<-ggplot(TN_dat,aes(x=Time,y=mean))+
  geom_bar(stat = 'identity',fill=c('darkgreen','pink3'),color='black',width = 0.45)+
  geom_errorbar(aes(ymin=mean,ymax=mean+se),width=.2)+
  geom_text(aes(label=sig,y=(mean+se)*1.02),position = position_dodge(0.9),vjust = 0)+
  scale_x_discrete(limits=c('c_TN_T0','c_TN_T28'))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.5))+
  xlab('Time')+ylab('TN (mg/L)')+
  theme_classic()


NH4_dat<-data.frame(t(sapply(pond_table[c('NH4_N_T0','NH4_N_T28')],function(x)
  (c(mean=mean(x),SD=sd(x))))))
se=NH4_dat$SD/sqrt(102)
with(pond_table,t.test(NH4_N_T0,NH4_N_T28,paired=T))
sig<-c('b','a')
NH4_dat<-data.frame(Time=rownames(NH4_dat),NH4_dat,se,sig,stringsAsFactors=F)


p3<-ggplot(NH4_dat,aes(x=Time,y=mean))+
  geom_bar(stat = 'identity',fill=c('darkgreen','pink3'),color='black',width = 0.45)+
  geom_errorbar(aes(ymin=mean,ymax=mean+se),width=.2)+
  geom_text(aes(label=sig,y=(mean+se)*1.02),position = position_dodge(0.9),vjust = 0)+
  scale_x_discrete(limits=c('NH4_N_T0','NH4_N_T28'))+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.18))+
  xlab('Time')+ylab('NH4_N (mg/L)')+
  theme_classic()



NO3_dat<-data.frame(t(sapply(pond_table[c('NO3_N_T0','NO3_N_T28')],function(x)
  (c(mean=mean(x),SD=sd(x))))))
se=NO3_dat$SD/sqrt(102)
with(pond_table,t.test(NO3_N_T0,NO3_N_T28,paired=T))
NO3_dat<-data.frame(Time=rownames(NO3_dat),NO3_dat,se,stringsAsFactors=F)


p4<-ggplot(NO3_dat,aes(x=Time,y=mean))+
  geom_bar(stat = 'identity',fill=c('darkgreen','pink3'),color='black',width = 0.45)+
  geom_errorbar(aes(ymin=mean,ymax=mean+se),width=.2)+
  scale_x_discrete(limits=c('NO3_N_T0','NO3_N_T28'))+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.18))+
  xlab('Time')+ylab('NO3_N (mg/L)')+
  theme_classic()


library(cowplot)
p<-plot_grid(p1,p2,p3,p4,
             labels = "AUTO", ncol = 2, nrow = 2, label_x = 0.015,
             label_y = 1.02,hjust = 0, label_size=17,align = "v")
pdf(file = 'E:/湖塘/数据处理/CN变化t检验.pdf',width=6.5,height= 6.5)
print(p)
dev.off()



vars <- colnames(pond_table)[c(5:21)]
model <- lapply(vars, function(x) {
  lm(substitute(log(ΔTOC)~i,list(i = as.name(x))), data = pond_table)})

r.squre<-round(as.vector(unlist(lapply(model, function(x) summary(x)$r.squared))),3)
p.value<-round(as.vector(unlist(lapply(model, function(x) anova(x)$'Pr(>F)'[1]))),3)
sig<-as.vector(unlist(lapply(p.value,formatPvalues)))
normality<-round(as.vector(unlist(lapply(model, function(x) shapiro.test(residuals(x))$p.value))),3)

summary_model<-data.frame(vars,r.squre,p.value,sig,normality)
summary_model


p5<-ggplot(pond_table,aes(x=c_TN_T0,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('TN (mg/L)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()


p6<-ggplot(pond_table,aes(x=SR_T0,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('SR')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()

p7<-ggplot(pond_table,aes(x=SUVA254_T0,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('SUVA254')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()

p8<-ggplot(pond_table,aes(x=Ca_T0,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('Ca (ppm)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()

p9<-ggplot(pond_table,aes(x=K_T0,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('K (ppm)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()


p10<-ggplot(pond_table,aes(x=Mg_T0,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('Mg (ppm)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()

p11<-ggplot(pond_table,aes(x=Na_T0,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('Na (ppm)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()

p12<-ggplot(pond_table,aes(x=P_T0,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('P (ppm)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()

p13<-ggplot(pond_table,aes(x=Al_T0,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('Al (ppb)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()


p14<-ggplot(pond_table,aes(x=Fe_T0,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('Fe (ppb)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()

p15<-ggplot(pond_table,aes(x=NH4_N_T0,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('NH4_N (mg/L)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()

p16<-ggplot(pond_table,aes(x=NO3_N_T0,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('NO3_N (mg/L)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()


p_ALL<-plot_grid(p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,
             labels = "AUTO", ncol = 4, nrow = 3, label_x = 0.015,
             label_y = 1.02,hjust = 0, label_size=17,align = "v")
pdf(file = 'E:/湖塘/数据处理/bdoc相关性检验.pdf',width=9.5,height= 6.5)
print(p_ALL)
dev.off()




##########差值相关性分析######
vars <- colnames(pond_table)[c(41:57)]
model <- lapply(vars, function(x) {
  lm(substitute(log(ΔTOC)~i,list(i = as.name(x))), data = pond_table)})
r.squre<-round(as.vector(unlist(lapply(model, function(x) summary(x)$r.squared))),3)
p.value<-round(as.vector(unlist(lapply(model, function(x) anova(x)$'Pr(>F)'[1]))),3)
sig<-as.vector(unlist(lapply(p.value,formatPvalues)))
normality<-round(as.vector(unlist(lapply(model, function(x) shapiro.test(residuals(x))$p.value))),3)

summary_model<-data.frame(vars,r.squre,p.value,sig,normality)
summary_model


p17<-ggplot(pond_table,aes(x=ΔTN,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('TN (mg/L)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()


p18<-ggplot(pond_table,aes(x=ΔSR,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('SR')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0),limits = c(-0.26,0.26))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()

p19<-ggplot(pond_table,aes(x=ΔSUVA254,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('SUVA254')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()

p20<-ggplot(pond_table,aes(x=ΔCa,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('Ca (ppm)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()

p21<-ggplot(pond_table,aes(x=ΔK,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('K (ppm)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()


p22<-ggplot(pond_table,aes(x=ΔMg,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('Mg (ppm)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()

p23<-ggplot(pond_table,aes(x=Δna,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('Na (ppm)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()

p24<-ggplot(pond_table,aes(x=ΔP,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('P (ppm)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()

p25<-ggplot(pond_table,aes(x=Δal,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('Al (ppb)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()


p26<-ggplot(pond_table,aes(x=ΔFe,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('Fe (ppb)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()

p27<-ggplot(pond_table,aes(x=ΔNH4.N,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('NH4_N (mg/L)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()

p28<-ggplot(pond_table,aes(x=ΔN03.N,y=log(ΔTOC)))+
  geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('N03_N (mg/L)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()


p_all<-plot_grid(p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,
                 labels = "AUTO", ncol = 4, nrow = 3, label_x = 0.015,
                 label_y = 1.02,hjust = 0, label_size=17,align = "v")
pdf(file = 'E:/湖塘/数据处理/bdoc差值相关性检验.pdf',width=9.5,height= 6.5)
print(p_all)
dev.off()



########序列数据分析#######
z3_table<-read.csv('E:/湖塘/数据处理/序列数据.csv',header=TRUE,stringsAsFactors = FALSE)
library(dplyr)
z3_bdoc<-data.frame(z3_table %>% group_by(Ponds) %>% summarise(ΔTOC=mean(ΔTOC)))
area<-c(130,41,104,84,146,69,31,8,5,2.8)
z3_bdoc<-cbind(z3_bdoc,area=area)
fit1<-lm(log(ΔTOC)~area,data=z3_bdoc)
summary(fit1)


p29<-ggplot(z3_bdoc,aes(x=area,y=log(ΔTOC)))+
  geom_point(shape=19, size=3,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('area (m2)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()
p29



library(dplyr)
pond_bdoc<-data.frame(pond_table %>% group_by(sites,Ponds) %>% summarise(ΔTOC=mean(ΔTOC)))
area<-c(10.5,12.5,36.7,2.1,68.5,29.8,46.2,35.6,130,41,104,84,146,69,31,8,5,
        2.8,1589,257,74,18,640,75,87,53,102,72,360,13,4.2,21,16,19)
pond_bdoc<-cbind(pond_bdoc,area=area)[-c(19,23),]
fit2<-lm(log(ΔTOC)~area,data=pond_bdoc)
summary(fit2)


p30<-ggplot(pond_bdoc,aes(x=area,y=log(ΔTOC)))+
  geom_point(shape=19, size=3,colour='tomato3',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  xlab('area (m2)')+ylab('Log BDOC (mg/L)')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.line=element_line(colour = 'black'))+
  theme_classic()
p30

####画图####
p_3<-plot_grid(p29,p30,
                 labels = "AUTO", ncol = 2, nrow = 1, label_x = 0.015,
                 label_y = 1.02,hjust = 0, label_size=17,align = "v")
pdf(file = 'E:/湖塘/数据处理/bdoc面积相关性检验.pdf',width=6.5,height= 3.3)
print(p_3)
dev.off()



####野外数据处理
pondfileddata<-read.csv('E:/湖塘/数据处理/pondfileddata.csv',header=TRUE,row.names = 1, stringsAsFactors = FALSE)
colnames(pondfileddata)
dat<-pondfileddata[-3]
melted <- melt(dat, id.vars=c("sites", "Ponds"))
results<-ddply(melted, c("sites", "variable"), summarise,
      mean = mean(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)))


library(lme4)
library(lmerTest)

lm.area<-lmer(Area~sites+(1|Ponds),pondfileddata)
library(multcomp)
summary(glht(lm.area, linfct = mcp(sites = "Tukey")), test = adjusted("holm"))

lm.PH<-lmer(PH~sites+(1|Ponds),pondfileddata)
library(multcomp)
summary(glht(lm.PH, linfct = mcp(sites = "Tukey")), test = adjusted("holm"))

lm.DO<-lmer(DO~sites+(1|Ponds),pondfileddata)
library(multcomp)
summary(glht(lm.DO, linfct = mcp(sites = "Tukey")), test = adjusted("holm"))

lm.EC<-lmer(EC~sites+(1|Ponds),pondfileddata)
library(multcomp)
summary(glht(lm.EC, linfct = mcp(sites = "Tukey")), test = adjusted("holm"))

lm.Salinity<-lmer(Salinity~sites+(1|Ponds),pondfileddata)
library(multcomp)
summary(glht(lm.Salinity, linfct = mcp(sites = "Tukey")), test = adjusted("holm"))

lm.Temperature<-lmer(Temperature~sites+(1|Ponds),pondfileddata)
library(multcomp)
summary(glht(lm.Temperature, linfct = mcp(sites = "Tukey")), test = adjusted("holm"))

lm.Depth<-lmer(Depth~sites+(1|Ponds),pondfileddata)
library(multcomp)
summary(glht(lm.Depth, linfct = mcp(sites = "Tukey")), test = adjusted("holm"))

lm.ALT<-lmer(ALT~sites+(1|Ponds),pondfileddata)
library(multcomp)
summary(glht(lm.ALT, linfct = mcp(sites = "Tukey")), test = adjusted("holm"))


lmm.results<-read.csv(file='E:/湖塘/数据处理/lmm_results.csv',header=TRUE,row.names = 1, stringsAsFactors = FALSE)


df<-list(subset(lmm.results,lmm.results$variable=='Area'),subset(lmm.results,lmm.results$variable=='PH'),
         subset(lmm.results,lmm.results$variable=='DO'),subset(lmm.results,lmm.results$variable=='EC'),
         subset(lmm.results,lmm.results$variable=='Salinity'),subset(lmm.results,lmm.results$variable=='Temperature'),
         subset(lmm.results,lmm.results$variable=='Depth'),subset(lmm.results,lmm.results$variable=='ALT'))

vars<-
plot.importance<-function(data){
  require(ggplot2)
  require(cowplot)
  plot<-function(x){
    ggplot(x,aes(x=sites,y= mean,fill=sites))+
      geom_bar(stat="identity")+
      geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),width=.2)+
      geom_text(aes(label=sig,y=(mean+sem)*1.02),position = position_dodge(0.9),vjust = 0)+
      #scale_y_continuous(expand = c(0,))+
      xlab('Sites')+
      theme_classic()+
      theme(axis.title.y =element_text(family = 'Times',size = 14),
            axis.text.y = element_text(family = 'Times',size = 14),
            axis.title.x = element_text(family = 'Times',size = 14),
            axis.text.x = element_text(family = 'Times',size = 14),
            legend.position='none')
  }
  images<-lapply(data,function(x) plot(x))
  cow.p<-plot_grid(images[[1]],images[[2]],images[[3]],images[[4]],images[[5]],
                   images[[6]],images[[7]],images[[8]],
                   labels = NULL, ncol = 4, nrow = 2, align = "v")
  print(cow.p)
}
plot.importance(df)




p<-ggplot(lmm.results,aes(x=sites,y=mean,fill=sites))+
  #geom_violin(trim=FALSE)+
  #geom_boxplot(width=0.1,fill="white")+
  geom_bar(stat = 'identity',width = 0.7)+
  geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),width=.2)+
  geom_text(aes(label=sig,y=(mean+sem)*1.02),position = position_dodge(0.9),vjust = 0)+
  facet_wrap(~variable,nrow = 2,scales = 'free_y')+
  #scale_y_continuous(expand = c(0,0))+
  xlab('Sites')+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=12),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 15),legend.position='none',element_blank())
p
dp<-ggplot(melted,aes(x=value,fill=variable))+
  geom_density( alpha = 0.5,size=1)+
  facet_wrap(~variable,nrow = 2,scales = 'free')+
  #scale_y_continuous(expand = c(0,0))+
  xlab('Sites')+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=12),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 15),legend.position='none',element_blank())
dp
