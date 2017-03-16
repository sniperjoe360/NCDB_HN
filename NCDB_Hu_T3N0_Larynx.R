# NCDB_HN.R: Data analysis for the NCDB Larynx Project with Dr. Hu and Givi
read<-read.csv("H:/NCDB Files/Larynx.csv",strip.white=TRUE)

# filter for only T2-3 NX M0 cases
dat <- read[(read$TNM_CLIN_T == 3 & (read$TNM_CLIN_N == 0 | read$TNM_CLIN_N == 1)) & read$TNM_CLIN_M == 0,]

###########################################
# Here are some relevant variables (see complete list seperately)
# [1] "PUF_CASE_ID"                 "FACILITY_TYPE_CD"           
# [3] "AGE"                         "DX_DEFSURG_STARTED_DAYS"                     
# [5] "RACE"                        "INSURANCE_STATUS"           
# [7] "YEAR_OF_DIAGNOSIS"           "HISTOLOGY"                  
# [9] "GRADE"                       "TNM_CLIN_T"                 
# [11] "TNM_CLIN_N"                  "TNM_CLIN_M"                 
# [13] "TNM_CLIN_STAGE_GROUP"        "TNM_PATH_T"                 
# [15] "TNM_PATH_N"                  "TNM_PATH_M"                 
# [17] "TNM_PATH_STAGE_GROUP"        "TNM_EDITION_NUMBER"         
# [19] "RX_SUMM_RADIATION"           "RAD_LOCATION_OF_RX"         
# [21] "RAD_TREAT_VOL"               "RAD_REGIONAL_RX_MODALITY"   
# [23] "RAD_REGIONAL_DOSE_CGY"       "RAD_BOOST_RX_MODALITY"      
# [25] "RAD_BOOST_DOSE_CGY"          "RAD_NUM_TREAT_VOL"          
# [27] "RX_SUMM_SURGRAD_SEQ"         "PUF_VITAL_STATUS"           
# [29] "DX_LASTCONTACT_DEATH_MONTHS" "RX_SUMM_IMMUNOTHERAPY"      
# [31] "RX_SUMM_SYSTEMIC_SUR_SEQ"    "DX_SYSTEMIC_STARTED_DAYS"   
# [33] "DX_CHEMO_STARTED_DAYS"       "CDCC_TOTAL"                 
# [35] "CROWFLY"                     "DX_RAD_STARTED_DAYS"        
#########################################
# install dependencies
install.packages('tableone')
# install.packages('survival')
# install.packages("ReporteRs")
# install.packages("magrittr")
library(tableone)
library(survival)
library(ReporteRs)
library(magrittr)

#########################################
# split data into comparative groups
C_Date <- dat$DX_CHEMO_STARTED_DAYS
R_Date <- dat$DX_RAD_STARTED_DAYS 
S_Date <- dat$DX_SURG_STARTED_DAYS

R_Leng <- dat$RAD_ELAPSED_RX_DAYS

CS_Seq <- dat$RX_SUMM_SYSTEMIC_SUR_SEQ 
RS_Seq <- dat$RX_SUMM_SURGRAD_SEQ

# 1A. Concurrent ChemoRT no Laryngectomy
# 2A. Induction/Sequential ChemoRT no Laryngectomy
# 3A. RT Only no Laryngectomy
# 4A. No therapy


CCRT <- dat[(C_Date < (R_Date + 7)) & (C_Date > (R_Date - 7)) & is.na(S_Date), ]
CCRT <- CCRT[!is.na(CCRT[1]),]
SCRT <- dat[(C_Date >= (R_Date + 7)) | (C_Date <= (R_Date - 7)) & is.na(S_Date), ]
SCRT <- SCRT[!is.na(SCRT[1]),]
RT   <- dat[is.na(C_Date) & !is.na(R_Date) & is.na(S_Date) ,]
None <- dat[is.na(C_Date) & is.na(R_Date) & is.na(S_Date),]

# 1B. Concurrent ChemoRT with Laryngectomy
# 2B. Induction/Sequential ChemoRT with Laryngectomy
# 3B. RT Only with Laryngectomy

CCRTL <- dat[(C_Date < (R_Date + 7)) & (C_Date > (R_Date - 7)) & !is.na(S_Date), ]
CCRTL <- CCRTL[!is.na(CCRTL[1]),]
SCRTL <- dat[(C_Date >= (R_Date + 7)) | (C_Date <= (R_Date - 7)) & !is.na(S_Date), ]
SCRTL <- SCRTL[!is.na(SCRTL[1]),]
RTL   <- dat[is.na(C_Date) & !is.na(R_Date) & !is.na(S_Date) ,]
L     <- dat[is.na(C_Date) & is.na(R_Date) & !is.na(S_Date),]

# final coding of the treatment groups
CCRT$RX <- 'CCRT'
CCRTL$RX<- 'CCRTL'
SCRT$RX <- 'SCRT'
SCRTL$RX<- 'SCRTL'
RT$RX   <- 'RT'
RTL$RX  <- 'RTL'
L$RX    <- 'L'

# final dataframe
larynx <- rbind(CCRT,CCRTL)
larynx <- rbind(larynx,SCRT)
larynx <- rbind(larynx,SCRTL)
larynx <- rbind(larynx,RT)
larynx <- rbind(larynx,RTL)
larynx <- rbind(larynx,L)

#recode for IMRT
larynx$IMRT <- "No IMRT"
larynx$IMRT[larynx$RAD_REGIONAL_RX_MODALITY == 31] = "IMRT"

# descriptive statistics
listVars<-(c("AGE", "SEX", "RACE","YEAR_OF_DIAGNOSIS","INSURANCE_STATUS","FACILITY_TYPE_CD","TNM_CLIN_T", "TNM_CLIN_N",'RAD_REGIONAL_RX_MODALITY','CROWFLY'))
catVars<-(c( "SEX", "TNM_CLIN_T",'TNM_CLIN_N','RAD_REGIONAL_RX_MODALITY','FACILITY_TYPE_CD','INSURANCE_STATUS','RACE','YEAR_OF_DIAGNOSIS')) 
table1 <- CreateTableOne(vars = listVars, data = larynx, factorVars = catVars, strata = 'RX')

docx( ) %>% 
     addFlexTable(table1 %>%
     FlexTable(add.rownames = TRUE ) %>%
     setZebraStyle( odd = "#DDDDDD", even = "#FFFFFF" ) ) %>%
     writeDoc(file = "H:/NCDB Files/larynx_table1_facility.docx")
	 
#########################################
# survival analysis for entire cohort and by treatment
attach(larynx)

#find the follow up and censor time
followup_time = larynx["DX_LASTCONTACT_DEATH_MONTHS"][[1]]
censor  = larynx["PUF_VITAL_STATUS"][[1]]
censor = 1-censor 

x11()

# create survival object
allpts <- survfit(Surv(followup_time, censor)~ 1)
plot(allpts, xlab="Time", ylab="Survival Probability (%)",lwd=2, main="Kaplan Meier Plot for Survival of All T3N0-1M0 Larynx")



x11()


# stratify by some variable
subset <- survfit(Surv(followup_time, censor)~RX)
plot(subset, xlab="Months", ylab="Survival Probability (%)", main='Overall Survival by Therapy for T3N0-1',lwd=2,col=c(1:length(names(subset$strata))),lty=1,cex.axis = 2,cex.lab = 1.5,cex.main = 2)
legend(80, 1.0, c("CCRT","CCRTL","L","RT",'RTL',"SCRT","SCRTL"),lwd=4,col=c(1:length(names(subset$strata))),lty=1)

# Cox Model
coxmodel <- coxph(Surv(followup_time, censor)~SEX + AGE + CDCC_TOTAL + as.factor(RX)+ as.factor(IMRT))
coxmodel

fit = summary(coxmodel)
hazard.ratio = fit$conf.int[,1]
lower95 = fit$conf.int[,3]
upper95 = fit$conf.int[,4]
pval   = fit$coefficients[,5]

finalCox<-cbind(hazard.ratio,lower95,upper95,pval)
write.csv(round(finalCox,2),file="H:/NCDB Files/T3N0-1_Larynx_IMRT.csv")


detach(larynx)

#########################################
# survival analysis for CCRT w/ and w/o IMRT
CCRT_ <- larynx[larynx$RX == 'CCRT' | larynx$RX == 'CCRTL' ,]
attach(CCRT_)

#find the follow up and censor time
followup_time = DX_LASTCONTACT_DEATH_MONTHS
censor  = PUF_VITAL_STATUS
censor = 1-censor 

x11()

# create survival object
allpts <- survfit(Surv(followup_time, censor)~ 1)
plot(allpts, xlab="Time", ylab="Survival Probability (%)",lwd=2, main="Kaplan Meier Plot for Survival of All T3N0-1 with CCRT")

x11()

# stratify by IMRT
subset <- survfit(Surv(followup_time, censor)~IMRT)
plot(subset, xlab="Months", ylab="Survival Probability (%)", main='Overall Survival of All CCRT +/- L T3N0-1 Patients \n Stratified by IMRT T3N0-1',lwd=2,col=c(1:length(names(subset$strata))),lty=1,cex.axis = 2,cex.lab = 1.5,cex.main = 2)
legend(80, 1.0, c("IMRT","Other"),lwd=4,col=c(1:length(names(subset$strata))),lty=1)
logrank <- survdiff(Surv(followup_time, censor)~IMRT)
pval <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
text(30,0.07,paste("logrank ","p =",pval,sep = ' '))

# Cox Model
coxmodel <- coxph(Surv(followup_time, censor)~SEX + AGE + CDCC_TOTAL + as.factor(IMRT))
coxmodel
detach(CCRT_)


fit = summary(coxmodel)
hazard.ratio = fit$conf.int[,1]
lower95 = fit$conf.int[,3]
upper95 = fit$conf.int[,4]
pval   = fit$coefficients[,5]

finalCox<-cbind(hazard.ratio,lower95,upper95,pval)
write.csv(round(finalCox,2),file="H:/NCDB Files/T3N0-1_Larynx_CCR IMRT.csv")



#########################################
# survival analysis for CCRT w/ and w/o IMRT
RT_ <- larynx[larynx$RX == 'RT',]
attach(RT_)

#find the follow up and censor time
followup_time = DX_LASTCONTACT_DEATH_MONTHS
censor  = PUF_VITAL_STATUS
censor = 1-censor 

x11()

# create survival object
allpts <- survfit(Surv(followup_time, censor)~ 1)
plot(allpts, xlab="Time", ylab="Survival Probability (%)",lwd=2, main="Kaplan Meier Plot for Survival of All T3N0-1 with RT Alone")

x11()

# stratify by IMRT
subset <- survfit(Surv(followup_time, censor)~IMRT)
plot(subset, xlab="Months", ylab="Survival Probability (%)", main='Overall Survival of All Definitive RT treated T3N0-1 Patients \n Stratified by IMRT',lwd=2,col=c(1:length(names(subset$strata))),lty=1,cex.axis = 2,cex.lab = 1.5,cex.main = 2)
legend(80, 1.0, c("IMRT","Other"),lwd=4,col=c(1:length(names(subset$strata))),lty=1)
logrank <- survdiff(Surv(followup_time, censor)~IMRT)
pval <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
text(30,0.07,paste("logrank ","p =",pval,sep = ' '))

# Cox Model
coxmodel <- coxph(Surv(followup_time, censor)~SEX + AGE + CDCC_TOTAL + as.factor(IMRT))
coxmodel
detach(RT_)


fit = summary(coxmodel)
hazard.ratio = fit$conf.int[,1]
lower95 = fit$conf.int[,3]
upper95 = fit$conf.int[,4]
pval   = fit$coefficients[,5]

finalCox<-cbind(hazard.ratio,lower95,upper95,pval)
write.csv(round(finalCox,2),file="H:/NCDB Files/T3N0-1_Larynx_IMRT_RT.csv")



#########################################
# survival analysis for CCRT vs. Laryngectomy
L_ <- larynx[larynx$RX == 'L' | larynx$RX == 'L' ,]
CCRT_ <- larynx[larynx$RX == 'CCRT' | larynx$RX == 'CCRT' ,]
L_CCRT_<-rbind(L_,CCRT_)
attach(L_CCRT_)



#find the follow up and censor time
followup_time = DX_LASTCONTACT_DEATH_MONTHS
censor  = PUF_VITAL_STATUS
censor = 1-censor 

x11()

# create survival object
allpts <- survfit(Surv(followup_time, censor)~ 1)
plot(allpts, xlab="Time", ylab="Survival Probability (%)",lwd=2, main="Kaplan Meier Plot for Survival of All T3N0-1 Patients with Laryngectomy or CCRT")

x11()

# stratify by Rx
subset <- survfit(Surv(followup_time, censor)~RX)
plot(subset, xlab="Months", ylab="Survival Probability (%)", main='Overall Survival of All T3N0-1 Patients \n Stratified by Concurrent Chemoradiation vs. Laryngectomy',lwd=2,col=c(1:length(names(subset$strata))),lty=1,cex.axis = 2,cex.lab = 1.5,cex.main = 2)
legend(80, 1.0, c("Concurrent CRT","Laryngectomy"),lwd=4,col=c(1:length(names(subset$strata))),lty=1)
logrank <- survdiff(Surv(followup_time, censor)~RX)
pval <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
text(30,0.07,paste("logrank","p =",pval,sep = ' '))

# Cox Model
coxmodel <- coxph(Surv(followup_time, censor)~SEX + AGE + CDCC_TOTAL + as.factor(RX))
coxmodel

detach(L_CCRT_)


fit = summary(coxmodel)
hazard.ratio = fit$conf.int[,1]
lower95 = fit$conf.int[,3]
upper95 = fit$conf.int[,4]
pval   = fit$coefficients[,5]

finalCox<-cbind(hazard.ratio,lower95,upper95,pval)
write.csv(round(finalCox,2),file="H:/NCDB Files/T3N0-1_Larynx_IMRT_CCRTvsL.csv")

detach(L_CCRT_)

# descriptive statistics
listVars<-(c("AGE", "SEX", "RACE","YEAR_OF_DIAGNOSIS",'CDCC_TOTAL'))
catVars<-(c( "SEX", "TNM_CLIN_T",'TNM_CLIN_N','YEAR_OF_DIAGNOSIS','CDCC_TOTAL')) 
table1 <- CreateTableOne(vars = listVars, data = L_CCRT_, factorVars = catVars, strata = 'RX')

#########################################
# propensity score matched cohort analysis

# KM survival analysis

# Cox Model

