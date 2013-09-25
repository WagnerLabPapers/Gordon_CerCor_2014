
# given a model, finds the significance of fixed effect specified with idx, using a log likelihood ratio test.
LRT <- function(theModel, idx=1){
	modFull = theModel
	
	modelTerms_h<-terms(theModel)
	modelTerms = attr(modelTerms_h, "term.labels")
	
	updateStr<-paste(".~.-",as.character(modelTerms[idx]))
	
	modNested = update(theModel, updateStr) 

  llFull = logLik(modFull)
	llNested = logLik(modNested)

	chsqval<-2*(llFull[1] - llNested[1])
	dfDiff<- attr(llFull,"df") - attr(llNested,"df") 

  res = list()
  res$pVal<-1-pchisq(chsqval,dfDiff)
	res$chsqval<-chsqval
  return(res)
}

#read in relevant libraries
library(lme4)

#read in raw data
dEE<-read.csv("AG1SubsequentMemoryMVPAAnalysisFS.csv") # for analyses involving encoding data only
dER<-read.csv("AG1_EEEvAndEREv_FS_withAnatHipp_EERSA.csv") # for analyses involving both encoding and retrieval data

#remove RTs that are too fast.
dEEClean<-dEE[dEE$RT>.3,]
dERClean<-dER[dER$RT>.3,]

#turn subsMem into 0s and 1s
dEEClean$subsMem <- dEEClean$subsMem-1

# add vectors of classifier correct status and of absolute encoding activity values
dEEClean$classifierCor<-as.numeric(dEEClean$EncAct>0)
dEEClean$absEncAct<-abs(dEEClean$EncAct)

#clear the res variable
remove(list = ls(pattern =  "res*"))

#combine the hippocampal values:
for (m in 1:27) {
dERClean$eHippBilat[dERClean$subs==m]<-(scale(dERClean$eHippLeft[dERClean$subs==m]) + scale(dERClean$eHippRight[dERClean$subs==m]))/2
dERClean$rHippBilat[dERClean$subs==m]<-(scale(dERClean$rHippLeft[dERClean$subs==m]) + scale(dERClean$rHippRight[dERClean$subs==m]))/2
}

################################################################
# ENC Analyses
################################################################

# Predict classifier accuracy with encoding strength
res.EE.classAccByEEStrength.mod<-glmer(classifierCor~absEncAct+class+subsMem+(absEncAct-1 | subs) + (1 | subs), data=dEEClean, family="binomial")

# Predict Subsequent Memory with Encoding Strength
res.EE.SMByEEStrength.mod<-glmer(subsMem~EncAct+class+ (EncAct-1 | subs) + (1 | subs), data=dEEClean, family="binomial")

# Predict RT with encoding strength (SM correct trials)
res.EE.RTByEEStrengthSMCor.mod<-lmer(RTs ~EncAct+class+(EncAct-1 | subs) + (1 | subs), data=dEEClean[dEEClean$subsMem==1,], REML=FALSE)

################################################################
# ENC/RET Analyses
################################################################

# predict ER (output off classifier trained on encoding, tested on retrieval) from EE (output of classifier trained and tested on Encoding)
res.ER.ERFromEE.mod<-lmer(ERAct ~EEAct+subsMem+class+(EEAct-1 | subs) + (1 | subs), data=dERClean, REML=FALSE)

# predict subsequent memory from ER
res.ER.SMFromER.mod<-glmer(subsMem ~ERActUnsigned+class+(ERActUnsigned-1 | subs) + (1 | subs), data=dERClean, family = "binomial")

# predict RT from ER - correct trials
res.ER.RTFromERCor.mod<-lmer(RT ~ERActUnsigned+class + (ERActUnsigned - 1 | subs) + (1 | subs), data=dERClean[dERClean$subsMem==1,], REML=FALSE)

################################################################
# Bilateral Hippocampal Analyses
################################################################

# predict Bilat Ret Hipp from EE Act
res.hipp.BilatHippAtRetFromEE.mod <-lmer(rHippBilat ~ EEActUnsigned + subsMem +class+(EEActUnsigned - 1 | subs) + (1 | subs), data=dERClean, REML=FALSE)

# predict Bilat Enc Hipp from EE Act
res.hipp.BilatHippAtEncFromEE.mod <-lmer(eHippBilat ~ EEActUnsigned + subsMem +class+(EEActUnsigned - 1 | subs) + (1 | subs), data=dERClean, REML=FALSE)

# predict Bilat Ret Hipp from Bilat Enc Hipp
res.hipp.BilatHippAtRetFromBilatHippAtEnc.mod <- lmer(rHippBilat ~ eHippBilat  +class+ subsMem + (eHippBilat - 1 | subs) + (1 | subs), data=dERClean, REML=FALSE)

# predict ER Act from Bilat Ret Hipp
res.hipp.ERFromBilatHippAtRet.mod <-lmer(ERActUnsigned ~ rHippBilat + subsMem +class+(rHippBilat - 1 | subs) + (1 | subs), data=dERClean, REML=FALSE)

# predict ER Act from Bilat Enc Hipp
res.hipp.ERFromBilatHippAtEnc.mod <-lmer(ERActUnsigned ~ eHippBilat + subsMem +class+(eHippBilat - 1 | subs) + (1 | subs), data=dERClean, REML=FALSE)

# predict SM from Bilat Enc Hipp
res.hipp.SMFromBilatHippAtEnc.mod <-glmer(subsMem ~ eHippBilat +class+(eHippBilat - 1 | subs) + (1 | subs), data=dERClean, family = "binomial", REML=FALSE)

# predict SM from Bilat Ret Hipp
res.hipp.SMFromBilatHippAtRet.mod <-glmer(subsMem ~ rHippBilat +class+(rHippBilat - 1 | subs) + (1 | subs), data=dERClean, family = "binomial", REML=FALSE)

# predict ER from Bilat Enc Hipp
res.hipp.ERHippFromEEHipp_Bilat.mod <-lmer(rHippBilat ~ eHippBilat + class + subsMem + (eHippBilat - 1 | subs) + (1 | subs), data=dERClean, REML=FALSE)

# predict Correct RT from Bilat Enc Hipp
res.hipp.RTFromBilatHippAtEnc.mod <-lmer(RT ~ eHippBilat +class+(eHippBilat - 1 | subs) + (1 | subs), data=dERClean[dERClean$subsMem==1,], REML=FALSE)

# predict Correct RT from Bilat Ret Hipp
res.hipp.RTFromBilatHippAtRet.mod <-lmer(RT ~ rHippBilat +class+(rHippBilat - 1 | subs) + (1 | subs), data=dERClean[dERClean$subsMem==1,], REML=FALSE)


################################################################
# Class Interactions
################################################################
res.EE.SMByEEStrength.ClassInteraction.mod<-glmer(subsMem~EncAct*class+ (EncAct:class-1 | subs) + (1 | subs), data=dEEClean, family="binomial")
res.EE.RTByEEStrengthSMCor.ClassInteraction.mod<-lmer(RTs ~EncAct*class+(EncAct:class-1 | subs) + (1 | subs), data=dEEClean[dEEClean$subsMem==1,], REML=FALSE)
res.ER.ERFromEE.ClassInteraction.mod<-lmer(ERAct ~EEAct*class+subsMem+(EEAct:class-1 | subs) + (1 | subs), data=dERClean, REML=FALSE)
res.ER.SMFromER.ClassInteraction.mod<-glmer(subsMem ~ERActUnsigned*class+(ERActUnsigned:class-1 | subs) + (1 | subs), data=dERClean, family = "binomial")
res.ER.RTFromERCor.ClassInteraction.mod<-lmer(RT ~ERActUnsigned*class + (ERActUnsigned:class - 1 | subs) + (1 | subs), data=dERClean[dERClean$subsMem==1,], REML=FALSE)
res.hipp.BilatHippAtRetFromEE.ClassInteraction.mod <-lmer(rHippBilat ~ EEActUnsigned*class + subsMem + (EEActUnsigned:class - 1 | subs) + (1 | subs), data=dERClean, REML=FALSE)
res.hipp.BilatHippAtEncFromEE.ClassInteraction.mod <-lmer(eHippBilat ~ EEActUnsigned*class + subsMem +(EEActUnsigned:class - 1 | subs) + (1 | subs), data=dERClean, REML=FALSE)
res.hipp.BilatHippAtRetFromBilatHippAtEnc.ClassInteraction.mod <- lmer(rHippBilat ~ eHippBilat*class + subsMem + (eHippBilat:class - 1 | subs) + (1 | subs), data=dERClean, REML=FALSE)
res.hipp.ERFromBilatHippAtRet.ClassInteraction.mod <-lmer(ERActUnsigned ~ rHippBilat*class + subsMem +(rHippBilat:class - 1 | subs) + (1 | subs), data=dERClean, REML=FALSE)
res.hipp.ERFromBilatHippAtEnc.ClassInteraction.mod <-lmer(ERActUnsigned ~ eHippBilat*class + subsMem +(eHippBilat:class - 1 | subs) + (1 | subs), data=dERClean, REML=FALSE)
res.hipp.SMFromBilatHippAtEnc.ClassInteraction.mod <-glmer(subsMem ~ eHippBilat*class +(eHippBilat:class - 1 | subs) + (1 | subs), data=dERClean, family = "binomial", REML=FALSE)
res.hipp.SMFromBilatHippAtRet.ClassInteraction.mod <-glmer(subsMem ~ rHippBilat*class +(rHippBilat:class - 1 | subs) + (1 | subs), data=dERClean, family = "binomial", REML=FALSE)
res.hipp.RTFromBilatHippAtEnc.ClassInteraction.mod <-lmer(RT ~ eHippBilat*class+(eHippBilat:class - 1 | subs) + (1 | subs), data=dERClean[dERClean$subsMem==1,], REML=FALSE)
res.hipp.RTFromBilatHippAtRet.ClassInteraction.mod <-lmer(RT ~ rHippBilat*class+(rHippBilat:class - 1 | subs) + (1 | subs), data=dERClean[dERClean$subsMem==1,], REML=FALSE)

################################################################
# Effects controlling for RSA measure of typicality
################################################################
res.EE.SMByEEStrength.RSAResids.mod<-glmer(subsMem~EEUnsRSAResid+class+ (EEUnsRSAResid-1 | subs) + (1 | subs), data=dRSAClean, family="binomial")
res.EE.RTByEEStrengthSMCor.RSAResids.mod<-lmer(RT ~EEUnsRSAResid+class+(EEUnsRSAResid-1 | subs) + (1 | subs), data=dRSAClean[dRSAClean$subsMem==1,], REML=FALSE)
res.ER.ERFromEE.RSAResids.mod<-lmer(ERUnsRSAResid ~EEUnsRSAResid+class+subsMem+(EEUnsRSAResid-1 | subs) + (1 | subs), data=dRSAClean, REML=FALSE)
res.ER.SMFromER.RSAResids.mod<-glmer(subsMem ~ERUnsRSAResid+class+(ERUnsRSAResid-1 | subs) + (1 | subs), data=dRSAClean, family = "binomial")
res.ER.RTFromERCor.RSAResids.mod<-lmer(RT ~ERUnsRSAResid+class + (ERUnsRSAResid - 1 | subs) + (1 | subs), data=dRSAClean[dERClean$subsMem==1,], REML=FALSE)
res.hipp.BilatHippAtRetFromEE.RSAResids.des  = "predict Bilat Ret Hipp from EE Act"
res.hipp.BilatHippAtRetFromEE.RSAResids.mod <-lmer(rHippBilat ~ EEUnsRSAResid+class + subsMem + (EEUnsRSAResid - 1 | subs) + (1 | subs), data=dRSAClean, REML=FALSE)
res.hipp.BilatHippAtEncFromEE.RSAResids.des = "predict Bilat Enc Hipp from EE Act"
res.hipp.BilatHippAtEncFromEE.RSAResids.mod <-lmer(eHippBilat ~ EEUnsRSAResid+class + subsMem +(EEUnsRSAResid - 1 | subs) + (1 | subs), data=dRSAClean, REML=FALSE)
res.hipp.BilatHippAtRetFromBilatHippAtEnc.RSAResids.des = "predict Bilat Ret Hipp from Bilat Enc Hipp"
res.hipp.BilatHippAtRetFromBilatHippAtEnc.RSAResids.mod <- lmer(rHippBilat ~ eHippBilat+class + subsMem + (eHippBilat - 1 | subs) + (1 | subs), data=dRSAClean, REML=FALSE)
res.hipp.ERFromBilatHippAtRet.RSAResids.des  = "predict ER Act from Bilat Ret Hipp"
res.hipp.ERFromBilatHippAtRet.RSAResids.mod <-lmer(ERUnsRSAResid ~ rHippBilat+class + subsMem +(rHippBilat - 1 | subs) + (1 | subs), data=dRSAClean, REML=FALSE)
res.hipp.ERFromBilatHippAtEnc.RSAResids.des  = "predict ER Act from Bilat Enc Hipp"
res.hipp.ERFromBilatHippAtEnc.RSAResids.mod <-lmer(ERUnsRSAResid ~ eHippBilat+class + subsMem +(eHippBilat - 1 | subs) + (1 | subs), data=dRSAClean, REML=FALSE)

################################################################
# Effects Driven By RSA mesaure of Typicality
################################################################
res.EE.SMByEEStrength.Typicality.mod<-glmer(subsMem~EERSA_CorClass+class+ (EERSA_CorClass-1 | subs) + (1 | subs), data=dRSAClean, family="binomial")
res.EE.RTByEEStrengthSMCor.Typicality.mod<-lmer(RT ~EERSA_CorClass+class+(EERSA_CorClass-1 | subs) + (1 | subs), data=dRSAClean[dRSAClean$subsMem==1,], REML=FALSE)
res.ER.ERFromEE.Typicality.mod<-lmer(ERRSA_CorClass ~EERSA_CorClass+class+subsMem+(EERSA_CorClass-1 | subs) + (1 | subs), data=dRSAClean, REML=FALSE)
res.ER.SMFromER.Typicality.mod<-glmer(subsMem ~ERRSA_CorClass+class+(ERRSA_CorClass-1 | subs) + (1 | subs), data=dRSAClean, family = "binomial")
res.ER.RTFromERCor.Typicality.mod<-lmer(RT ~ERRSA_CorClass+class + (ERRSA_CorClass - 1 | subs) + (1 | subs), data=dRSAClean[dRSAClean$subsMem==1,], REML=FALSE)
res.hipp.BilatHippAtRetFromEE.Typicality.des  = "predict Bilat Ret Hipp from EE Act"
res.hipp.BilatHippAtRetFromEE.Typicality.mod <-lmer(rHippBilat ~ EERSA_CorClass+class + subsMem + (EERSA_CorClass - 1 | subs) + (1 | subs), data=dRSAClean, REML=FALSE)
res.hipp.BilatHippAtEncFromEE.Typicality.des = "predict Bilat Enc Hipp from EE Act"
res.hipp.BilatHippAtEncFromEE.Typicality.mod <-lmer(eHippBilat ~ EERSA_CorClass+class + subsMem +(EERSA_CorClass - 1 | subs) + (1 | subs), data=dRSAClean, REML=FALSE)
res.hipp.ERFromBilatHippAtRet.Typicality.des  = "predict ER Act from Bilat Ret Hipp"
res.hipp.ERFromBilatHippAtRet.Typicality.mod <-lmer(ERRSA_CorClass ~ rHippBilat+class + subsMem +(rHippBilat - 1 | subs) + (1 | subs), data=dRSAClean, REML=FALSE)
res.hipp.ERFromBilatHippAtEnc.Typicality.des  = "predict ER Act from Bilat Enc Hipp"
res.hipp.ERFromBilatHippAtEnc.Typicality.mod <-lmer(ERRSA_CorClass ~ eHippBilat+class + subsMem +(eHippBilat - 1 | subs) + (1 | subs), data=dRSAClean, REML=FALSE)

################################################################
# Effects with hippocampus unstandardized
################################################################

res.hipp.BilatHippAtRetFromEE.Unstandardized.mod <-lmer(rHippBilatUnstd ~ EEActUnsigned + class + subsMem + (EEActUnsigned - 1 | subs) + (1 | subs), data=dRSAClean, REML=FALSE)
res.hipp.BilatHippAtEncFromEE.Unstandardized.mod <-lmer(eHippBilatUnstd ~ EEActUnsigned + class + subsMem +(EEActUnsigned - 1 | subs) + (1 | subs), data=dRSAClean, REML=FALSE)
res.hipp.BilatHippAtRetFromBilatHippAtEnc.Unstandardized.mod <- lmer(rHippBilatUnstd ~ eHippBilatUnstd + class + subsMem + (eHippBilatUnstd - 1 | subs) + (1 | subs), data=dRSAClean, REML=FALSE)
res.hipp.ERFromBilatHippAtRet.Unstandardized.mod <-lmer(ERActUnsigned ~ rHippBilatUnstd + class + subsMem +(rHippBilatUnstd - 1 | subs) + (1 | subs), data=dRSAClean, REML=FALSE)
res.hipp.ERFromBilatHippAtEnc.Unstandardized.mod <-lmer(ERActUnsigned ~ eHippBilatUnstd + class + subsMem +(eHippBilatUnstd - 1 | subs) + (1 | subs), data=dRSAClean, REML=FALSE)
res.hipp.SMFromBilatHippAtEnc.Unstandardized.mod <-glmer(subsMem ~ eHippBilatUnstd + class +(eHippBilatUnstd - 1 | subs) + (1 | subs), data=dRSAClean, family = "binomial", REML=FALSE)
res.hipp.SMFromBilatHippAtRet.Unstandardized.mod <-glmer(subsMem ~ rHippBilatUnstd + class +(rHippBilatUnstd - 1 | subs) + (1 | subs), data=dRSAClean, family = "binomial", REML=FALSE)
res.hipp.RTFromBilatHippAtEnc.Unstandardized.mod <-lmer(RT ~ eHippBilatUnstd + class+(eHippBilatUnstd - 1 | subs) + (1 | subs), data=dRSAClean[dRSAClean$subsMem==1,], REML=FALSE)
res.hipp.RTFromBilatHippAtRet.Unstandardized.mod <-lmer(RT ~ rHippBilatUnstd + class+(rHippBilatUnstd - 1 | subs) + (1 | subs), data=dRSAClean[dRSAClean$subsMem==1,], REML=FALSE)

################################################################
# pVals - linear models
################################################################
allModels<-ls(pattern = glob2rx("res.*.mod"))

for (m in 1:length(allModels)) {
	thisModName<-substr(allModels[m],1,nchar(allModels[m]))
	assign(paste(thisModName, ".pVal", sep = ""), LRT(eval(as.name(thisModName)),1))		
}

################################################################
# Mediation analysis
################################################################

# standardize all continuous variables.  
dCleanStandardized = dERClean
dCleanStandardized$rGenHipp<-scale(dCleanStandardized$rHippBilat)
dCleanStandardized$ERActUnsigned<-scale(dCleanStandardized$ERActUnsigned)
dCleanStandardized$EEActUnsigned<-scale(dCleanStandardized$EEActUnsigned)
dCleanStandardized$RT<-scale(dCleanStandardized$RT)
dCleanStandardized$eGenHipp<-scale(dCleanStandardized$eHippBilat)

# base model, from which we can get the magnitude of direct effects
eq.1<-glmer(subsMem ~ rGenHipp + ERActUnsigned + EEActUnsigned + class + (rGenHipp - 1 | subs) + (ERActUnsigned - 1 | subs) +(EEActUnsigned - 1 | subs) + (1 | subs), data=dCleanStandardized, family = "binomial", REML=FALSE)
eq.2<-lmer(ERActUnsigned ~ EEActUnsigned + rGenHipp + class  + (EEActUnsigned - 1 | subs)  + (rGenHipp - 1 | subs) + (1 | subs), data=dCleanStandardized, REML=FALSE)
eq.3<-lmer(rGenHipp ~ eGenHipp + class  + subsMem + (eGenHipp - 1 | subs) + (1 | subs) , data=dCleanStandardized, REML=FALSE)
eq.4<-lmer(eGenHipp ~ EEActUnsigned + class  + subsMem + (EEActUnsigned - 1 | subs) + (1 | subs) , data=dCleanStandardized, REML=FALSE)

# direct effects
DE1<-fixef(eq.1)[2] 
DE2<-fixef(eq.1)[3] 
DE3<-fixef(eq.1)[4] 
DE4<-fixef(eq.2)[2] 
DE5<-fixef(eq.2)[3] 
DE6<-fixef(eq.3)[2] 
DE7<-fixef(eq.4)[2]

# use bootstrapping to find the significance of indirect effects.
bootDist1<-0
bootDist2<-0
bootDist3<-0
bootDist4<-0
bootDist5<-0
bootDist6<-0
bootDist7<-0

for (i in 1:10000) {
  
  # randomly sample from data, with replacement
  o1<-dCleanStandardized[sample(1:dim(dCleanStandardized)[1], size=dim(dCleanStandardized)[1], replace=T),]
   
  # all trials
  beq.1<-glmer(subsMem ~ rGenHipp + ERActUnsigned + EEActUnsigned  + class + (rGenHipp - 1 | subs) + (ERActUnsigned - 1 | subs) +(EEActUnsigned - 1 | subs) + (1 | subs), data=o1, family = "binomial", REML=FALSE)
  beq.2<-lmer(ERActUnsigned ~ EEActUnsigned + rGenHipp  +  class + subsMem + (EEActUnsigned - 1 | subs)  + (rGenHipp - 1 | subs) + (1 | subs), data=o1, REML=FALSE)
  beq.3<-lmer(rGenHipp ~ eGenHipp  + class + subsMem + (eGenHipp - 1 | subs) + (1 | subs) ,data=o1, REML=FALSE)
  beq.4<-lmer(eGenHipp ~ EEActUnsigned + class + subsMem + (EEActUnsigned - 1 | subs) + (1 | subs) , data=o1, REML=FALSE)
  
  bootDist1[i]<-fixef(beq.1)[2] 
  bootDist2[i]<-fixef(beq.1)[3] 
  bootDist3[i]<-fixef(beq.1)[4] 
  bootDist4[i]<-fixef(beq.2)[2] 
  bootDist5[i]<-fixef(beq.2)[3] 
  bootDist6[i]<-fixef(beq.3)[2] 
  bootDist7[i]<-fixef(beq.4)[2]
}

#distribution for indirect effects 
distIE1 = bootDist5 * bootDist2
distIE2 = bootDist4 * bootDist2
