setwd("C:/Users/aziz/Desktop/Canada DA")
require(mlbench)
require(car)
require(MASS)
require(plyr)
require(mice)
data("BreastCancer")
BC<-BreastCancer
BC.NA<-na.omit(BC)
BC.Ben<-subset(BC.NA,BC.NA$Class=="benign")
BC.Mal<-subset(BC.NA,BC.NA$Class=="malignant")
summary(BC.Mal$Normal.nucleoli)
summary(BC.Ben$Normal.nucleoli)
#require(foreign)
###write.csv(BC.NA, "C:/Users/aziz/Desktop/Canada DA/BC.fix.csv")
BC.fix<-read.csv("C:/Users/aziz/Desktop/Canada DA/BC.fix.csv")
summary(BC.fix)
#######################################################################
Flagging the missing data
#######################################################################
flag.Clt<-is.na(BC.fix$Cl.thickness)
flag.Clt<-as.numeric(flag.Clt)
flag.Clt[flag.Clt==0]<-5
flag.Clt[flag.Clt==1]<-4
flag.Clt<-flag.Clt-4

flag.csz<-is.na(BC.fix$Cell.size)
flag.csz<-as.numeric(flag.csz)
flag.csz[flag.csz==0]<-5
flag.csz[flag.csz==1]<-4
flag.csz<-flag.csz-4

flag.csh<-is.na(BC.fix$Cell.shape)
flag.csh<-as.numeric(flag.csh)
flag.csh[flag.csh==0]<-5
flag.csh[flag.csh==1]<-4
flag.csh<-flag.csh-4

flag.Mad<-is.na(BC.fix$Marg.adhesion)
flag.Mad<-as.numeric(flag.Mad)
flag.Mad[flag.Mad==0]<-5
flag.Mad[flag.Mad==1]<-4
flag.Mad<-flag.Mad-4

flag.Ecs<-is.na(BC.fix$Epith.c.size)
flag.Ecs<-as.numeric(flag.Ecs)
flag.Ecs[flag.Ecs==0]<-5
flag.Ecs[flag.Ecs==1]<-4
flag.Ecs<-flag.Ecs-4

flag.Bn<-is.na(BC.fix$Bare.nuclei)
flag.Bn<-as.numeric(flag.Bn)
flag.Bn[flag.Bn==0]<-5
flag.Bn[flag.Bn==1]<-4
flag.Bn<-flag.Bn-4

flag.Bl<-is.na(BC.fix$Bl.cromatin)
flag.Bl<-as.numeric(flag.Bl)
flag.Bl[flag.Bl==0]<-5
flag.Bl[flag.Bl==1]<-4
flag.Bl<-flag.Bl-4

flag.Nn<-is.na(BC.fix$Normal.nucleoli)
flag.Nn<-as.numeric(flag.Nn)
flag.Nn[flag.Nn==0]<-5
flag.Nn[flag.Nn==1]<-4
flag.Nn<-flag.Nn-4

flag.Mit<-is.na(BC.fix$Mitoses)
flag.Mit<-as.numeric(flag.Mit)
flag.Mit[flag.Mit==0]<-5
flag.Mit[flag.Mit==1]<-4
flag.Mit<-flag.Mit-4
# 1 missing 0 not missing


BC.fix<-cbind(BC.fix,flag.Clt,flag.csz,flag.csh,flag.Mad,flag.Ecs,flag.Bn,flag.Bl,flag.Nn,flag.Mit)
BC.fix<-BC.fix[c(1:3,13,4,14,5,15,6,16,7,17,8,18,9,19,10,20,11,21,12)]
View(BC.fix)
####################################
Casewise Deletion
######################################
BC.C<-na.omit(n.PID)
summary(BC.C)
BC.CBen<-subset(BC.C,BC.C$class=="benign")
summary(BC.CBen)
BC.CMal<-subset(BC.C,BC.C$class=="malignant")
summary(BC.CMal)
########################################################
Pairwise
#########################################################
cor(BC.C,use="pairwise.complete.obs") 
summary(N.pid,use="pairwise.complete.obs") 
mean(N.pid,na.rm=TRUE)
#########################################################
Sub group Mean imputation
Sub group being Class
#########################################################
impute.mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
BC.fix2 <- ddply(BC.fix, ~Class, transform, Cl.thickness = impute.mean(Cl.thickness),
     Cell.size = impute.mean(Cell.size),Cell.shape=impute.mean(Cell.shape),
Marg.adhesion=impute.mean(Marg.adhesion),Epith.c.size=impute.mean(Epith.c.size),
Bare.nuclei=impute.mean(Bare.nuclei),Bl.cromatin=impute.mean(Bl.cromatin),
Normal.nucleoli=impute.mean(Normal.nucleoli),Mitoses=impute.mean(Mitoses))
BC.fix2<-BC.fix2[order(BC.fix2$X), ]
View(BC.fix2)
###########################################################
Mean Imputations
#################################################################
BC.Mean<-BC.fix
BC.Mean$Cl.thickness<-as.numeric(BC.Mean$Cl.thickness)
BC.Mean$Cl.thickness[BC.Mean$Cl.thickness==NA]<-mean(BC.Mean$Cl.thickness,na.rm=T)

BC.Mean$Cl.thickness[is.na(BC.Mean$Cl.thickness)]<-4.436
BC.Mean$Cell.size[is.na(BC.Mean$Cell.size)]<-3.138
BC.Mean$Cell.shape[is.na(BC.Mean$Cell.shape)]<-3.223
BC.Mean$Marg.adhesion[is.na(BC.Mean$Marg.adhesion)]<-2.832
BC.Mean$Epith.c.size[is.na(BC.Mean$Epith.c.size)]<-3.235
BC.Mean$Bare.nuclei[is.na(BC.Mean$Bare.nuclei)]<-3.551
BC.Mean$Bl.cromatin[is.na(BC.Mean$Bl.cromatin)]<-3.428
BC.Mean$Normal.nucleoli[is.na(BC.Mean$Normal.nucleoli)]<-2.872
BC.Mean$Mitoses[is.na(BC.Mean$Mitoses)]<-1.608
View(BC.Mean)
###########################################################
Case Mean Imputations !!!!!!!!!!!!Need to work on this!!!!!!!!!!!!!! 
#################################################################
BC.row<-BC.NA
BC.row[is.na(BC.row)]<-rowMeans((BC.row))
###########################################################
Medium Imputations
#################################################################
impute <- function(x, fun) {
 missing <- is.na(x)
  replace(x, missing, fun(x[!missing]))
}
BC.fix3 <- ddply(BC.fix, ~Class, transform, Cl.thickness = impute(Cl.thickness,median),
     Cell.size = impute(Cell.size,median),Cell.shape=impute(Cell.shape,median),
Marg.adhesion=impute(Marg.adhesion,median),Epith.c.size=impute(Epith.c.size,median),
Bare.nuclei=impute(Bare.nuclei,median),Bl.cromatin=impute(Bl.cromatin,median),
Normal.nucleoli=impute(Normal.nucleoli,median),Mitoses=impute(Mitoses,median))
BC.fix3<-BC.fix3[order(BC.fix3$X), ]
View(BC.fix3)
##################################################################
Min Imputation
##################################################################
BC.fix4 <- ddply(BC.fix, ~Class, transform, Cl.thickness = impute(Cl.thickness,min),
     Cell.size = impute(Cell.size,min),Cell.shape=impute(Cell.shape,min),
Marg.adhesion=impute(Marg.adhesion,min),Epith.c.size=impute(Epith.c.size,min),
Bare.nuclei=impute(Bare.nuclei,min),Bl.cromatin=impute(Bl.cromatin,min),
Normal.nucleoli=impute(Normal.nucleoli,min),Mitoses=impute(Mitoses,min))
BC.fix4<-BC.fix4[order(BC.fix4$X), ]
View(BC.fix4)
##################################################################
Max Imputation
##################################################################
BC.fix5<- ddply(BC.fix, ~Class, transform, Cl.thickness = impute(Cl.thickness,max),
     Cell.size = impute(Cell.size,max),Cell.shape=impute(Cell.shape,max),
Marg.adhesion=impute(Marg.adhesion,max),Epith.c.size=impute(Epith.c.size,max),
Bare.nuclei=impute(Bare.nuclei,max),Bl.cromatin=impute(Bl.cromatin,max),
Normal.nucleoli=impute(Normal.nucleoli,max),Mitoses=impute(Mitoses,max))
BC.fix5<-BC.fix5[order(BC.fix5$X), ]
View(BC.fix5)

##################################################################
Regression Imputation
##################################################################
BC.fixed<-BC.fix[,-c(1,2,4,6,8,10,12,14,16,18,20,21)]

fit<-lm(Cl.thickness~.,data=BC.fixed)
COEF<-fit$coefficients
MIS<-BC.fixed[is.na(Cl.thickness),]
REPLACE<-ceiling(COEF[1]+(rowSums(COEF[-1]*MIS[,-1],na.rm=T)))
BC.fixed$Cl.thickness[which(is.na(BC.fixed$Cl.thickness))]<-REPLACE

LM.impute<-function(Cl,BC){
	fit<-lm(BC[,Cl]~.,data=BC)
	COEF<-fit$coefficients
	MIS<-BC[which(is.na(BC[,Cl])),]
	REPLACE<-ceiling(COEF[1]+(rowSums(COEF[-1]*MIS[,-1],na.rm=T)))
	BC[which(is.na(BC[,Cl])),Cl]<-REPLACE
	return(BC[,Cl])
}

reg.imp<-cbind(
	LM.impute(1,BC.fixed),
	LM.impute(2,BC.fixed),
	LM.impute(3,BC.fixed),
	LM.impute(4,BC.fixed),
	LM.impute(5,BC.fixed),
	LM.impute(6,BC.fixed),
	LM.impute(7,BC.fixed),
	LM.impute(8,BC.fixed),
	LM.impute(9,BC.fixed))
colnames(reg.imp)<-names(BC.fixed)
regg.imp<-cbind(BC.fix[,2],reg.imp,BC.NA$Class)
colnames(regg.imp)<-names(BC.NA)
regg.imp[,11][regg.imp[,11]==1]<-'benign'
regg.imp[,11][regg.imp[,11]==2]<-'malignant'
View(regg.imp)
##############################################################################








