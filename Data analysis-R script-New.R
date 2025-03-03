##data analysis for manuscript-Xiaocui Wang et al.
#Packages used
library(readxl)
library(tidyverse)
library(lme4)
library(aods3)
library(car)
library(emmeans)
library(xlsx)
library(lmerTest)
library(optimx)


# Data analysis for Figure 1-Feeding --------------------------------------


#import feeding data
Fed<-Figure_1_Feeding
#full model
#include starvationTime as random effects
#lg transform the feeding data and center lg data by adding its mean
Fed1<-mutate(Fed,Yeast=factor(Yeast),
             Population=factor(Population),
             Date=factor(Date),
             StarvationTime=factor(StarvationTime),
             YeastBatch=factor(YeastBatch),
             YeastDay=factor(YeastDay),
             lgYeastConsumption=log(YeastConsumption+0.06))
m1<-lmer(lgYeastConsumption~Yeast*Population+(1|Date)+(1|YeastBatch)+(1|YeastDay)+(1|StarvationTime),data=Fed1)
summary(m1)
Anova(m1)
drop1(m1)
ranova(m1)
plot(m1)
qqnorm(residuals(m1))
#Conclusion: starvationtime has significant effects on feeding. 
#Therefore, we need to include it as fixed effects
#Second model to include starvationtime as fixed effects 
#and then change starvation time as numeric data
#center starvation time by substracting mean
Fed2<-mutate(Fed,Yeast=factor(Yeast),
             Population=factor(Population),
             Date=factor(Date),
             StarvationTime=as.numeric(StarvationTime),
             STC=StarvationTime-mean(StarvationTime),
             YeastBatch=factor(YeastBatch),
             YeastDay=factor(YeastDay),
             lgYeastConsumption=log(YeastConsumption+0.06))
#roughly check the normality of lg transformed feeding data
hist(Fed2$lgYeastConsumption)
#run the second model
m2<-lmer(lgYeastConsumption~Yeast*Population+STC+
          +(1|YeastBatch/Date)+(1|YeastDay/Date),
         data=Fed2,
         control = lmerControl(
           optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))
summary(m2)
Anova(m2)
ranova(m2)
drop1(m2)
step(m2)
plot(m2)
qqnorm(residuals(m2))
leveneTest(YeastConsumption~Yeast*Population,data=Fed2)
#multiple comparison to compare population differences, 
#yeast differences and fly-yeast interaction
emmeans(m2, pairwise ~ Yeast | Population, adjust = "tukey",type="response")
emmeans(m2, pairwise ~ Population | Yeast, adjust = "tukey",type="response")
emmeans(m2, pairwise ~ Yeast : Population, adjust = "tukey",type="response")


# Data analysis for Figure 2 -Mating --------------------------------------


##Data analysis for Figure 2 -Mating
#import mating data
Mat<-Figure_2_Mating
##Data analysis for Figure 2a -First mating latency
Mat$VirginalMatingLatency=as.numeric(Mat$VirginalMatingLatency)
Mat1<-subset(Mat,VirginalMatingLatency>0)
Mat2<-mutate(Mat1,lgVML=log(VirginalMatingLatency),
             Yeast=factor(Yeast),
             Population=factor(Population),
             Date=factor(Date),
             DishNumber=factor(DishNumber),
             YeastBatch=factor(YeastBatch),
             YeastDay=factor(YeastDay))
m3<-lmer(lgVML~Yeast*Population+(1|YeastBatch/Date)+(1|YeastDay/Date)+
           (1|DishNumber),data=Mat2)
summary(m3)
Anova(m3)
plot(m3)
qqnorm(residuals(m3))
hist(Mat2$lgVML)
leveneTest(lgVML~Yeast*Population,data=Mat2)
#Testing random part of model
m32<-update(m3,~.-(1|YeastBatch/Date))
anova(m3,m32)
m33<-update(m3,~.-(1|YeastDay/Date))
anova(m3,m33)
m34<-update(m3,~.-(1|DishNumber))
anova(m3,m34)

#multiple comparison
emmeans(m3, pairwise ~ Yeast | Population, adjust = "tukey",type="response")
emmeans(m3, pairwise ~ Population | Yeast, adjust = "tukey",type="response")
emmeans(m3, pairwise ~ Yeast:Population, adjust = "tukey",type="response")

##Data analysis for Figure 2b -second-mating probability
Mat$Remating=as.numeric(Mat$Remating)
Mat3<-subset(Mat,Remating>=0)
Mat4<-mutate(Mat3,Yeast=factor(Yeast),
             Population=factor(Population),
             Date=factor(Date),
             DishNumber=factor(DishNumber),
             YeastBatch=factor(YeastBatch),
             YeastDay=factor(YeastDay))
m4<-glmer(Remating~Yeast*Population+(1|DishNumber)+(1|YeastBatch/Date)+(1|YeastDay/Date),
          family=binomial,
          control=glmerControl(calc.derivs=F),
          data=Mat4)
summary(m4)
Anova(m4)
plot(m4)
qqnorm(residuals(m4))
leveneTest(Remating~Yeast*Population,data=Mat4)
#Testing random part of model
m41<-update(m4,~.-(1|Date))
anova(m4,m41)
m42<-update(m4,~.-(1|YeastBatch/Date))
anova(m4,m42)
m43<-update(m4,~.-(1|YeastDay/Date))
anova(m4,m43)
m44<-update(m4,~.-(1|DishNumber))
anova(m4,m44)

#Multiple comparison
emmeans(m4, pairwise ~ Yeast, adjust = "tukey",type="response")
emmeans(m4, pairwise ~ Population, adjust = "tukey",type="response")


##Data analysis for Figure 2c -No of total matings
m5<-glmer(MatingFrequency~Yeast*Population+(1|DishNumber)+(1|YeastBatch/Date)+
            (1|YeastDay/Date),family=poisson,data=Mat)
summary(m5)
gof(m5)
Anova(m5)
plot(m5,type=c("p","smooth"))
qqnorm(residuals(m5))
leveneTest(MatingFrequency~Yeast*Population,data=Mat)

#test random effects
m51<-update(m5,~.-(1|Date))
anova(m5,m51)
m52<-update(m5,~.-(1|YeastBatch/Date))
anova(m5,m52)
m53<-update(m5,~.-(1|YeastDay/Date))
anova(m5,m53)
m54<-update(m5,~.-(1|DishNumber))
anova(m5,m54)

#multiple comparison
emmeans(m5, pairwise ~ Yeast , adjust = "tukey",type="response")
emmeans(m5, pairwise ~ Population , adjust = "tukey",type="response")


# Data analysis for Figure 3-Egg-laying -----------------------------------


##Data analysis for Figure 3-Egg-laying
#import egg-laying data
Ovi<-Figure_3_Egg_laying
colnames(Ovi)[7]<-"Egg"
#transform data for analysis
Ovi1<- mutate(Ovi,Egg=as.numeric(Egg),
              Yeast=factor(Yeast),
              Population=factor(Population),
              Date=factor(Date),
              DishNumber=factor(DishNumber),
              YeastBatch=factor(YeastBatch),
              YeastDay=factor(YeastDay))
#Run the GLMM model
m6<-glmer.nb(Egg~Yeast*Population+(1|DishNumber)+(1|YeastBatch/Date)+(1|YeastDay/Date),
             data=Ovi1,
             control=glmerControl(calc.derivs=F))
summary(m6)
Anova(m6)
plot(m6)
qqnorm(residuals(m6))
#Testing random part of model
m61<-update(m6,~.-(1|Date))
anova(m6,m61)
m62<-update(m6,~.-(1|YeastBatch/Date))
anova(m6,m62)
m63<-update(m6,~.-(1|YeastDay/Date))
anova(m6,m63)
m64<-update(m6,~.-(1|DishNumber))
anova(m6,m64)

#multiple comparison
emmeans(m6, pairwise ~ Yeast | Population, adjust = "tukey",type="response")
emmeans(m6, pairwise ~ Population| Yeast, adjust = "tukey",type="response")
emmeans(m6, pairwise ~ Yeast:Population, adjust = "tukey",type="response")

# Data analysis for Figure 4-Egg development ------------------------------


##Data analysis for Figure 4-Egg development
#import development data
Dev<-Figure_4_Egg_development
colnames(Dev)[8]<-"Adult"
colnames(Dev)[9]<-"Adu50"

##Data analysis for Figure 4a-Egg-to-adult time
#analysis for adult time
Dev1<-subset(Dev,Adu50>=0)
m7<-lmer(Adu50~Yeast*Population+Location+(1|Date)+(1|Time),data=Dev1)
summary(m7)
Anova(m7)
plot(m7)
qqnorm(residuals(m7))
#Testing random part of model
ranova(m7)
#multiple comparison
emmeans(m7, pairwise ~ Location,adjust = "tukey",type = "response")
emmeans(m7, pairwise ~ Yeast|Population,adjust = "tukey",type = "response")
emmeans(m7, pairwise ~ Population|Yeast, adjust = "tukey",type = "response")
emmeans(m7, pairwise ~ Yeast:Population, adjust = "tukey",type = "response")

##Data analysis for Figure 4b-Egg-to-adult survival
Dev2<-mutate(Dev,Yeast=factor(Yeast),
             Population=factor(Population),
             Time=factor(Time),
             Date=factor(Date),
             DishNumber=factor(DishNumber),
             Location=factor(Location))
y<-cbind(Dev2$Adult,10-Dev2$Adult)
m8<-glmer(y~Yeast*Population+Date+(1|Time)+
            (1|DishNumber)+(1|Location),
          family=binomial,data=Dev2,
          glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(m8)
Anova(m8)
plot(m8)
qqnorm(residuals(m8))
#test random effects
m82<-update(m8,~.-(1|Time))
anova(m8,m82)
m83<-update(m8,~.-(1|DishNumber))
anova(m8,m83)
m84<-update(m8,~.-(1|Location))
anova(m8,m84)

#multiple comparison
emmeans(m8, pairwise ~ Date,adjust = "tukey",type = "response")
emmeans(m8, pairwise ~ Yeast|Population,adjust = "tukey",type = "response")
emmeans(m8, pairwise ~ Population|Yeast,adjust = "tukey",type = "response")
emmeans(m8, pairwise ~ Yeast:Population,adjust = "tukey",type = "response")


# Data analysis for Figure S1 ---------------------------------------------


##Data analysis for Figure S1a-Egg-to-pupae time
Dev3<-subset(Dev,Pup50>=0)
m9<-lmer(Pup50~Yeast*Population+(1|Location)+(1|Date)+(1|Time),data=Dev3)
summary(m9)
Anova(m9)
plot(m9)
qqnorm(residuals(m9))
#Testing random part of model
ranova(m9)
#multiple comparison
emmeans(m9, pairwise ~ Yeast|Population,adjust = "tukey",type = "response")
emmeans(m9, pairwise ~ Population|Yeast, adjust = "tukey",type = "response")
emmeans(m9, pairwise ~ Yeast:Population, adjust = "tukey",type = "response")


##Data analysis for Figure S1b-Egg-to-pupae survival
y<-cbind(Dev2$Pupae,10-Dev2$Pupae)
m10<-glmer(y~Yeast*Population+Date+(1|Time)+
            (1|DishNumber)+(1|Location),
          family=binomial,data=Dev2,
          glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(m10)
Anova(m10)
plot(m10)
qqnorm(residuals(m10))
#test random effects
m102<-update(m10,~.-(1|Time))
anova(m10,m102)
m103<-update(m10,~.-(1|DishNumber))
anova(m10,m103)
m104<-update(m10,~.-(1|Location))
anova(m10,m104)

#multiple comparison
emmeans(m10, pairwise ~ Date,adjust = "tukey",type = "response")
emmeans(m10, pairwise ~ Yeast|Population,adjust = "tukey",type = "response")
emmeans(m10, pairwise ~ Population|Yeast,adjust = "tukey",type = "response")
emmeans(m10, pairwise ~ Yeast:Population,adjust = "tukey",type = "response")

