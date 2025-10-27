# Year- and species-specific analyses for each species

library(glmmTMB)
library(ggplot2)
require(plyr) # for ddply
require(ggeffects)
require(DHARMa)
require(car)
require(gridExtra)
require(dplyr)

rm(list=ls())

setwd()

#### PEAK RECRUITMENT ANALYSES ####
#read raw recruitment data

peak2 = read.csv("Peak_recruitment_plants.csv", sep=";", dec=",", header = T)
peak2

str(peak2)

table(peak2$forest, peak2$species)
hist(peak2$absrec2)

#include olre
peak2$olre = seq.int(nrow(peak2))

# test for ZI
m0 = glmmTMB(cbind(absrec2, denom2)~1+(1|plot), data=peak2, family=binomial)
m1 = glmmTMB(cbind(absrec2, denom2)~(1|species)+(1|plot), data=peak2, family=binomial)
m2 = glmmTMB(cbind(absrec2, denom2)~forest+(1|species)+(1|plot), data=peak2, family=binomial)
m3 = glmmTMB(cbind(absrec2, denom2)~forest+(1|species)+(1|plot), data=peak2, family=betabinomial)
# zi-formula for forest seems the best

m4 = glmmTMB(cbind(absrec2, denom2)~forest+(1|species)+(1|plot), data=peak2, family=betabinomial, ziformula=~forest)
m5 = glmmTMB(cbind(absrec2, denom2)~forest+(1|olre)+(1|species)+(1|plot), data=peak2, family=binomial, ziformula=~forest)

anova(m0,m1,m2,m3,m4,m5)
# olre > betabinomial

m6 = glmmTMB(cbind(absrec2, denom2)~forest*veg1+(1|olre)+(1|species)+(1|plot), data=peak2, family=binomial, ziformula=~forest)
m7 = glmmTMB(cbind(absrec2, denom2)~forest*veg1*z.cc+(1|olre)+(1|species)+(1|plot), data=peak2, family=binomial, ziformula=~forest)
m8 = glmmTMB(cbind(absrec2, denom2)~forest*veg1*z.cc+(1|species)+(1|plot), data=peak2, family=betabinomial, ziformula=~forest)
anova(m7,m8)
# ziformula forest -> robust
# olre > betabinomial -> robust

# do species respond differently to certain factors? #keep it maximal
m7 = glmmTMB(cbind(absrec2, denom2)~forest*veg1*z.cc+(1|olre)+(1|species)+(1|plot), data=peak2, family=binomial)
m8 = glmmTMB(cbind(absrec2, denom2)~forest*veg1*z.cc+(1|olre)+(1+forest|species)+(1|plot), data=peak2, family=binomial)

anova(m7, m8)
# random effect forest*species > ziformula forest
# further model assumptions not needed.
# Note: We also observed the model outcome and residual fit for the above-described models
# m8 seemed to fit the data best.
# continue with m8

plot(simulateResiduals(fittedModel = m8, plot = F))
# residuals are fine
Anova(m8)

# check model robustness against removal of z.cc
m9 = glmmTMB(cbind(absrec2, denom2)~forest*veg1+(1|olre)+(1|species)+(1|plot), data=peak2, family=binomial) 
Anova(m9)
plot(ggpredict(m9, terms=c("veg1[all]")))
plot(ggpredict(m9, terms=c("forest", "species"), type="random"))
# results are robust and qualitatively similar

#### Plotting ####
res3 = ggpredict(m8, terms=c("forest"))
res6 = ggpredict(m8, terms=c("veg1[all]", "forest"))

# fixed
# for species-specific lines in the graph
mfixed = glmmTMB(cbind(absrec2, denom2)~forest*veg1*z.cc*species+(1|olre)+(1|plot), data=peak2, family=binomial)
plot(simulateResiduals(fittedModel = mfixed, plot = F))
Anova(mfixed)
# qualitatively similar results across models with species as fixed and random factor

res72 = ggpredict(mfixed, terms=c("veg1[all]", "species", "forest"))

# take the predicted mean of both forest types
res74 = res72 %>%
  group_by(group, x) %>%
  summarise(predicted = mean(predicted, na.rm=T),
            conf.high = mean(conf.high, na.rm=T),
            conf.low = mean(conf.low, na.rm=T))

# take the predicted mean of both forest types
res6.1 = res6 %>%
  group_by(x) %>%
  summarise(predicted = mean(predicted, na.rm=T),
            conf.high = mean(conf.high, na.rm=T),
            conf.low = mean(conf.low, na.rm=T))

# scale colour gradient
cols8 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols5 = c("#999999", "#E69F00", "#F0E442", "#D55E00", "#CC79A7")

pforest2 = ggplot(res3, aes(x, predicted, group))+
  geom_point(res3, mapping=aes(x), size=3)+
  geom_linerange(res3, mapping=aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high), size=1)+
  scale_y_continuous(limits=c(0,0.08), name=c("Recruitment rate"), breaks=c(seq(0,1,0.02)))+
  scale_x_discrete(name = c("Forest type"), labels=c("Degraded", "Intact"))+
  theme_classic()#+theme(legend.position = "none")

fun_veg = function (x) (x - mean(spring2$bk1/3, na.rm=T))/ sd(spring2$bk1/3, na.rm=T)
fun_backveg = function (x) x *sd(spring2$bk1/3, na.rm=T) + mean(spring2$bk1/3, na.rm=T) 

# Note: The main effect is shown for species as a random factor.
# The species-dependent effect is with species as a fixed factor.
pveg = ggplot()+
  # geom_ribbon(res6.1, mapping=aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high), fill="lightgrey", alpha=0.5)+
  geom_smooth(res6.1, mapping=aes(x=x, y=predicted), method = "lm", formula = y~x+I(x^2), se=F, colour="black")+
  geom_line(res74, mapping=aes(x=x, y=predicted, colour=group), stat="smooth", method="lm", formula=y~x+I(x^2), 
            alpha=0.5, size=1)+
  scale_y_continuous(limits=c(0,0.45), name=c("Recruitment rate"), breaks=c(seq(0,1,0.05)))+
  scale_x_continuous(name = c("Ground vegetation [%]"), 
                     breaks=c(-1.7150416, -0.6785534, 0.3579347, 1.3944229, 2.4309110),
                     labels=c(10, 20, 30, 40, 50),
                     limits=c(-2.233, 2.431))+
  scale_colour_manual(values=cols8, name="Plant species",
                      breaks=c("Eueu", "Fral", "Prpa", "Rhca", "Risp", "Rini", "Soau", "Viop"),
                      labels=c("E. europaeus", "F. alnus", "P. padus", "R. cathartica", "R. spicatum",
                               "R. nigrum", "S. aucuparia", "V. opulus"))+
  theme_classic()+theme(legend.position = c(0.95, 0.95),
                        legend.justification=c("right", "top"),
                        legend.text=element_text(face="italic"),
                        legend.key.size = unit(1, "line"))+
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8))+
  guides(color = guide_legend(override.aes = list(size = 1)))+
  guides(shape = guide_legend(override.aes = list(size = 1)))

#
pforest2
pveg

grid.arrange(pforest2, pveg, nrow=1)

#### EARLY SURVIVAL ####

early0 = read.csv("Early_survival_plants.csv", header=T, sep=";", dec=",")

# Note:
table(early0$species, early0$surv)

# too few Prpa
# too few Rini surviving
# in addition, Rhca & Sani need to be removed, because the seedlings were spatially clustered in only two locations

early = early0[early0$species != "Prpa" & early0$species != "Sani" & early0$species != "Rhca" &
                 early0$species != "Rini",] 

# Note:
# Because Survival is binomial, using plots as random variable would lead to an overfitting of the
# model outcome at best, at worst in convergence issues.
# Thus, we only used the sites as random effects in the models

n7 = glmmTMB(surv~forest*z.veg*z.cc*species+(1|sites), early, family=binomial)
plot(simulateResiduals(fittedModel = n7, plot = F))
testOutliers(n7, type=c("bootstrap")) 
# seems fine
Anova(n7)

# Plotting model output of n7
nres1 = ggpredict(n7, terms=c("z.cc[all]", "species", "forest"))
nres1$group <- factor(nres1$group, levels=c("Eueu","Fral","Risp","Soau","Viop"))

# pool forest types
nres1.1 = nres1 %>%
  group_by(x, group) %>%
  summarise(predicted = mean(predicted, na.rm=T),
            conf.high = mean(conf.high, na.rm=T),
            conf.low = mean(conf.low, na.rm=T))

# pool species
nres1.2 = nres1.1 %>%
  group_by(x) %>%
  summarise(predicted = mean(predicted, na.rm=T),
            conf.high = mean(conf.high, na.rm=T),
            conf.low = mean(conf.low, na.rm=T))

untransform_canopy =  function (x) x*sd(canopy$cc, na.rm=T)+mean(canopy$cc, na.rm=T)

psurv = ggplot(nres1, aes(x, predicted))+
  # geom_ribbon(nres1.2, mapping=aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high), fill="lightgrey", alpha=0.5)+
  geom_line(nres1.2, mapping=aes(x=x, y=predicted), stat="smooth", method="glm", 
            method.args = list(family = "binomial"), size=1)+
  geom_line(nres1, mapping=aes(x=x, y=predicted, colour=group), stat="smooth", method="glm", 
            method.args = list(family = "binomial"), alpha=0.5, size=1)+
  scale_x_continuous(name="Canopy cover [%]",
                     breaks=c(-2.2946203, -1.4483703, -0.6021203,  0.2441297,  1.0903797,  1.9366297),
                     labels=c(70,75,80,85,90,95),
                     limits=c(-2.2946203,1.9366297))+
  geom_jitter(early, mapping=aes(x=z.cc, y=surv, colour=species), height=0.04,width=0.02, alpha=0.4)+
  scale_y_continuous(name="Early survival")+
  scale_colour_manual(values=cols5, name="Plant species",
                      breaks=c("Eueu" ,"Fral", "Risp", "Soau", "Viop"),
                      labels=c("E. europaeus", "F. alnus", "R. spicatum","S. aucuparia", "V. opulus"))+
  theme_classic()+theme(legend.position = "none",
                        legend.justification=c("right", "top"),
                        legend.text=element_text(face="italic"))

p1 = grid.arrange(pforest2, pveg, psurv, nrow=1)
p1
#ggsave("Supplementary_figure.png", p1, dpi=450, units=c("mm"), width=300, height=120)

#end#