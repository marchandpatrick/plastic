#####################################################################################################
#####     Factors influencing the detection of beach plastic debris	            #################
#####################################################################################################
## cite as: Lavers, JL; Oppel, S; Bond, AL (2016). Factors influencing the detection of beach plastic debris. Marine Environmental Research XXX
## written by steffen.oppel@rspb.org.uk on 12 June 2015
## EXCLUDED black fragments on 22 Aug due to obvious violation of critical assumptions (no false positives)


############### LOAD THE REQUIRED PACKAGES #########################################################

library(unmarked)
library(Hmisc)


#####################################################################################################################################################
#############    LOAD RAW DATA FROM GITHUB      ##############################################################################
#####################################################################################################################################################
setInternet2(TRUE)
con <- url("https://github.com/steffenoppel/plastic/blob/master/plastic.RData?raw=true")
load(con)
close(con)


### DATA 'dat2': data.frame that lists the number of fragments of a certain type and colour recorded by each observer in each quadrat
head(dat2)

### DATA 'validat2': data.frame that lists the true number of fragments of a certain type and colour recorded by sediment extraction in each quadrat
head(validat2)

### DATA 'SITEDAT': data.frame that lists the environmental covariates for each quadrat
head(SITEDAT)




#######################################################################################################
#############    SPECIFY PLASTIC 'SPECIES' WE WANT TO ANALYSE      ####################################
#######################################################################################################

types<-c('Fragment','Nurdles')
colours<-c('White','Blue','Black')






#####################################################################################################################################################
#############    ANALYSIS: EXPLORE MULTIPLE EXPLANATORY COVARIATES FOR DETECTION AND ESTIMATE TRUE ABUNDANCE      ###################################
#####################################################################################################################################################

output<-data.frame()
detout<-data.frame()
AICtables<-data.frame()


for (t in 1:length(types)){
typedat<-subset(dat2, Type==types[t])
valitypedat<-subset(validat2, Type==types[t])

for (col in 1:length(colours)){
coldat<-subset(typedat, Colour==colours[col])
valicoldat<-subset(valitypedat, Colour==colours[col])

if(dim(coldat)[1]>0){



##### PREPARE OBSERVATION COVARIATE DATA ####
OBSDAT<-coldat[order(coldat$Quadrat, coldat$Observer),]
OBSDAT$experience<-ifelse(OBSDAT$Observer %in% c('JL','AF'),1,0)
OBSDAT$Observer<-as.factor(OBSDAT$Observer)
head(OBSDAT)
str(OBSDAT)
dim(OBSDAT)


##### PREPARE COUNT DATA OF PLASTIC DEBRIS ####
obs.y<-cast(coldat, Quadrat~Observer, value='Number')
obs.y<-obs.y[order(obs.y$Quadrat),]
nobs<-dim(obs.y)[2]			### number of observers = repeat visits
obs.y<-as.matrix(obs.y[,2:nobs])

##### PREPARE SITE COVARIATE DATA ####
## black fragments were counted in trials B and C, the rest in A and B

if(t==1 & col==2){
SITEDAT2<-subset(SITEDAT, trial %in% c("B","C"))}else{
SITEDAT2<-subset(SITEDAT, trial %in% c("A","B"))}
SITEDAT2<-SITEDAT2[order(SITEDAT2$Quadrat),]
SITEDAT2$Quadrat<-as.factor(SITEDAT2$Quadrat)
SITEDAT2$a<-seq(1:dim(SITEDAT2)[1])

##### COMBINE RESPONSE AND OBSERVATION COVARIATES TO UNMARKED FRAME AND STANDARDIZE NUMERIC COVARIATES #######

UFP<-unmarkedFramePCount(y=obs.y, siteCovs=SITEDAT2, obsCovs=OBSDAT[,c(4,5,6,9)])
siteCovs(UFP)[,2:3] <- scale(siteCovs(UFP)[,2:3])
obsCovs(UFP)[,2:3] <- scale(obsCovs(UFP)[,2:3])
summary(UFP)

### model 6 does not work for black fragments, which were not counted in trial A, so sun is invariant
if(dim(OBSDAT)[1]==65){
obsCovs(UFP)[,5] <- rnorm(65,0)		### introduce random number for sun covariate
}


###################################################################################################
######## ANALYSIS OF DATA #########################################################################
###################################################################################################

### set K (upper limit of integration) to the maximum count for the sieved data

Kset<-max((max(valicoldat$Number)+1),(max(coldat$Number)+1))			### take the maximum value of either sieved data OR observations, whichever is larger


########### FIT PLAUSIBLE BINOMIAL MIXTURE MODELS WITH VARIOUS DETECTION COVARIATES ############
### first formula is for detection, second formula is for abundance  ####

m0 <- pcount( ~ 1 ~ as.factor(a), data= UFP,  mixture = "P", K=Kset)
m1 <- pcount( ~ Observer ~ as.factor(a), data= UFP,  mixture = "P", K=Kset)
m2 <- pcount( ~ Observer+BIOL ~ as.factor(a), data= UFP,  mixture = "P", K=Kset)
m3 <- pcount( ~ Observer+Order ~ as.factor(a), data= UFP,  mixture = "P", K=Kset)
m4 <- pcount( ~ Observer+STONES ~ as.factor(a), data= UFP,  mixture = "P", K=Kset)
m5 <- pcount( ~ STONES+BIOL ~ as.factor(a), data= UFP,  mixture = "P", K=Kset)
m6 <- pcount( ~  Observer+Sun ~ as.factor(a), data= UFP,  mixture = "P", K=Kset)
m7 <- pcount( ~  experience ~ as.factor(a), data= UFP,  mixture = "P", K=Kset)
m8 <- pcount( ~  experience+BIOL ~ as.factor(a), data= UFP,  mixture = "P", K=Kset)
m9 <- pcount( ~  experience+Order ~ as.factor(a), data= UFP,  mixture = "P", K=Kset)
m10 <- pcount( ~  experience+STONES ~ as.factor(a), data= UFP,  mixture = "P", K=Kset)
m11 <- pcount( ~  experience+Sun ~ as.factor(a), data= UFP,  mixture = "P", K=Kset)

fl <- fitList(null=m0, 'Observer'=m1, 'obs+biol.debris'=m2, 'obs+fatigue'=m3, 'obs+stones'=m4,'biol.debris+stones'=m5,
			'obs+visibility'=m6, 'experience'=m7, 'experience+biol.debris'=m8, 'experience+fatigue'=m9, 'experience+stones'=m10, 'experience+visibility'=m11)


### SAVE MODEL SELECTION OUTPUT ###############
ms <- modSel(fl, nullmod="null")
aicout<-ms@Full
aicout$Type<-types[t]
aicout$Colour<-colours[col]
AICtables<-rbind(AICtables,aicout[,c(1,91:100)])		### NEED TO CHANGE THIS IF LAMBDA FORMULA CHANGES!!


### COLLATE BEST UNBIASED PREDICTED ABUNDANCE FROM TOP MODEL ###
out<-predict(fl,type='state')											## discarded once state formula reduced to ~1


names(out)[3:4]<-c('lcl','ucl')
out$Type<-types[t]
out$Colour<-colours[col]
if(dim(out)[1]==13){
out$Quadrat<-SITEDAT$Quadrat[SITEDAT$trial=="B"]
out$REAL<-valicoldat$Number[valicoldat$trial=="B"]
}else{
out$Quadrat<-SITEDAT2$Quadrat
out$REAL<-valicoldat$Number
}
output<-rbind(output,out)


### COLLATE MODEL AVERAGED PREDICTED DETECTION PROBABILITY ###
detprobs<-predict(fl,type='det', newdat=data.frame(experience=c(1,0,0,1,0), Observer=c('JL', 'SO', 'LM', 'AF', 'AD'), BIOL=5, STONES=5, Order=1, Sun=1))
detprobs$experience<-c(1,0,0,1,0)
detprobs$Observer<-c('JL', 'SO', 'LM', 'AF', 'AD')
detprobs$Type<-types[t]
detprobs$Colour<-colours[col]
detout<-rbind(detout,detprobs)


}			### close if loop that excludes non-existent combinations

}			### close loop over colours

}			### close loop over fragments		


dim(output)
		





###################################################################################################
######## OUTPUT TABLES ##################################################################
###################################################################################################

### predicted abundance of plastic particles per quadrat
head(output)

### model selection table for each type and colour of plastic
head(AICtables)

### predicted detection probability of plastic particles for each observer
head(detout)


############# FIGURE ######################

output<-read.table("predicted_abundances_QUAD_experience_final.csv", header=T, sep=",")
output$ucl[is.infinite(output$ucl)] <- 0 
pd<-aggregate(Predicted~Type+Colour,output, FUN=sum)
pd$lcl<-aggregate(lcl~Type+Colour,output, FUN=sum)[,3]
pd$ucl<-aggregate(ucl~Type+Colour,output, FUN=sum)[,3]
pd$REAL<-aggregate(REAL~Type+Colour,output, FUN=sum)[,3]

par(mar=c(6,6,0,0))
errbar(c(1:4), pd$Predicted[c(3,4,5,2)], pd$lcl[c(3,4,5,2)], pd$ucl[c(3,4,5,2)],xlim=c(0,5), ylim=c(0,800), axes=F,ylab="Total number of plastic particles", xlab="",cex.lab=1.7,
pch=16, cex=2, mgp=c(4,1,0))
axis(1, at=c(0:5), labels=c("","Blue \n Fragments","White \n Fragments","White \n Nurdles","Black \n Nurdles",""), cex.axis=1.5, lwd=1, lwd.ticks=0, pos=0, padj=1)
axis(2, at=seq(0,800,200), labels=T, cex.axis=1.5, las=1)
points(pd$REAL[c(3,4,5,2)], pch=4, col='red', cex=2.5)


