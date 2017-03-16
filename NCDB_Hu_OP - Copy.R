# NCDB_HN.R: Data analysis for the NCDB Larynx Project with Dr. Hu and Givi
read<-read.csv("H:/NCDB Files/Larynx.csv",strip.white=TRUE)

# filter for only T2-3 NX M0 cases
dat <- read[(read$TNM_CLIN_T == 2 | read$TNM_CLIN_T  == 3) & read$TNM_CLIN_M == 0,]

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
# install.packages('tableone')
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

x11()

# create survival object
allpts <- survfit(Surv(followup_time, censor)~ 1)
plot(allpts, xlab="Time", ylab="Survival Probability (%)",lwd=2, main="Kaplan Meier Plot for Survival of All T2-3NXM0 Larynx")

x11()


# stratify by some variable
subset <- survfit(Surv(followup_time, censor)~RX)
plot(subset, xlab="Months", ylab="Survival Probability (%)", main='Overall Survival by Therapy',lwd=2,col=c(1:length(names(subset$strata))),lty=1,cex.axis = 2,cex.lab = 1.5,cex.main = 2)
legend(80, 1.0, c("CCRT","CCRTL","L","RT",'RTL',"SCRT","SCRTL"),lwd=4,col=c(1:length(names(subset$strata))),lty=1)

# Cox Model
coxmodel <- coxph(Surv(followup_time, censor)~SEX + AGE + CDCC_TOTAL + as.factor(RX)+ as.factor(IMRT)+ as.factor(droplevels(TNM_CLIN_T))+ as.factor(droplevels(TNM_CLIN_N)))
coxmodel
detach(larynx)

#########################################
# survival analysis for CCRT w/ and w/o IMRT
CCRT_ <- larynx[larynx$RX == 'CCRT' | larynx$RX == 'CCRTL' ,]
attach(CCRT_)

#find the follow up and censor time
followup_time = DX_LASTCONTACT_DEATH_MONTHS
censor  = PUF_VITAL_STATUS

x11()

# create survival object
allpts <- survfit(Surv(followup_time, censor)~ 1)
plot(allpts, xlab="Time", ylab="Survival Probability (%)",lwd=2, main="Kaplan Meier Plot for Survival of All Patients with CCRT")

x11()

# stratify by IMRT
subset <- survfit(Surv(followup_time, censor)~IMRT)
plot(subset, xlab="Months", ylab="Survival Probability (%)", main='Overall Survival of All CCRT +/- L Patients \n Stratified by IMRT',lwd=2,col=c(1:length(names(subset$strata))),lty=1,cex.axis = 2,cex.lab = 1.5,cex.main = 2)
legend(80, 1.0, c("IMRT","Other"),lwd=4,col=c(1:length(names(subset$strata))),lty=1)
logrank <- survdiff(Surv(followup_time, censor)~IMRT)
pval <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
text(30,0.07,paste("logrank ","p<0.001",sep = ' '))

# Cox Model
coxmodel <- coxph(Surv(followup_time, censor)~SEX + AGE + CDCC_TOTAL + as.factor(IMRT)+ as.factor(droplevels(TNM_CLIN_T))+ as.factor(droplevels(TNM_CLIN_N)))
coxmodel
detach(CCRT_)
#########################################
# survival analysis for SCRT w/ and w/o IMRT
SCRT_ <- larynx[larynx$RX == 'SCRT' | larynx$RX == 'SCRTL' ,]
attach(SCRT_)



#find the follow up and censor time
followup_time = DX_LASTCONTACT_DEATH_MONTHS
censor  = PUF_VITAL_STATUS

x11()

# create survival object
allpts <- survfit(Surv(followup_time, censor)~ 1)
plot(allpts, xlab="Time", ylab="Survival Probability (%)",lwd=2, main="Kaplan Meier Plot for Survival of All Patients with SCRT")

x11()

# stratify by IMRT
subset <- survfit(Surv(followup_time, censor)~IMRT)
plot(subset, xlab="Months", ylab="Survival Probability (%)", main='Overall Survival of All SCRT +/- L Patients \n Stratified by IMRT',lwd=2,col=c(1:length(names(subset$strata))),lty=1,cex.axis = 2,cex.lab = 1.5,cex.main = 2)
legend(80, 1.0, c("IMRT","Other"),lwd=4,col=c(1:length(names(subset$strata))),lty=1)
logrank <- survdiff(Surv(followup_time, censor)~IMRT)
pval <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
text(30,0.07,paste("IMRT vs. Other logrank ","p<0.001",sep = ' '))

# Cox Model
coxmodel <- coxph(Surv(followup_time, censor)~SEX + AGE + CDCC_TOTAL + as.factor(IMRT)+ as.factor(droplevels(TNM_CLIN_T))+ as.factor(droplevels(TNM_CLIN_N)))
coxmodel



#########################################
# propensity score matched cohort analysis

# KM survival analysis

# Cox Model

