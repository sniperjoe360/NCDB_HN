# NCDB_HN.R: Data analysis for the NCDB Larynx Project with Dr. Hu and Givi
read<-read.csv("H:/NCDB Files/Larynx.csv",strip.white=TRUE)

dim(read)

# filter for only T1-2 NX M0 cases
dat <- read[

(
(read$TNM_CLIN_T == 1 | 
 read$TNM_CLIN_T == '1A' |
 read$TNM_CLIN_T == '1B' |
  read$TNM_CLIN_T == 2)),]

dat2<-dat[!is.na(dat$TNM_CLIN_T),]

dat3 <- dat2[dat2$TNM_CLIN_N == 0 & !is.na(dat2$TNM_CLIN_N),] 

dat4 <- dat3[dat3$TNM_CLIN_M != 1,]

# only supraglottic larynx
dat5 <- dat4[dat4$PRIMARY_SITE == 'C321',]


# Exclude patients with previous cancer history
dats2 <- dat5[dat5$SEQUENCE == 1 | 
     dat5$SEQUENCE == 0 
,]
dim(dats2)


# Exclude patients with no radiation and those with chemo
dats3 <- dats2[dats2$RX_SUMM_RADIATION != 0 & dats2$RX_SUMM_CHEMO == 0 & dats2$RX_SUMM_SURG_PRIM_SITE == 0,]
dim(dats3)


# include patients with complete survival data
dat7 <- dats3[!is.na(dats3$DX_LASTCONTACT_DEATH_MONTHS)& !is.na(dats3$PUF_VITAL_STATUS),]

# Exclude patients who died within 3 months
dat8 <- dat7[(dat7$DX_LASTCONTACT_DEATH_MONTHS < 3 & (dat7$PUF_VITAL_STATUS != 0))|dat7$DX_LASTCONTACT_DEATH_MONTHS >= 3 ,]
dim(dat8)

larynx <- dat8
dim(larynx)

#install.packages("party")
 # install.packages("ReporteRs")
 # install.packages("magrittr")
#install.packages("forestmodel")
#install.packages("tableone")
#install.packages("survival")
#install.packages('MatchIt')
 # FOREST PLOT
library("forestmodel")
library(tableone)
library(party)
library(survival)
library(MatchIt)
#library(ReporteRs)

############################

#recode for IMRT
larynx$IMRT <- "1. Non-IMRT Radiation"
larynx$IMRT[larynx$RAD_REGIONAL_RX_MODALITY == 31] = "2. IMRT"

larynx$IMRT_Bin <- FALSE
larynx$IMRT_Bin[larynx$RAD_REGIONAL_RX_MODALITY == 31] = TRUE



larynx$NODE_CAT="Unknown SURVIVAL DATA"
larynx[(!is.na(larynx$DX_LASTCONTACT_DEATH_MONTHS)& !is.na(larynx$PUF_VITAL_STATUS)),'NODE_CAT'] = 'FULL SURVIVAL DATA'  

larynx$RX_SUMM_RADIATION_Recode = NULL
larynx[larynx$RX_SUMM_RADIATION ==1 ,'RX_SUMM_RADIATION_Recode'] = "EBRT"
larynx[larynx$RX_SUMM_RADIATION ==2 ,'RX_SUMM_RADIATION_Recode'] = "Unknown"
larynx[larynx$RX_SUMM_RADIATION == 3 ,'RX_SUMM_RADIATION_Recode'] = "Unknown"
larynx[larynx$RX_SUMM_RADIATION == 4 ,'RX_SUMM_RADIATION_Recode'] = "EBRT"
larynx[larynx$RX_SUMM_RADIATION == 5 ,'RX_SUMM_RADIATION_Recode'] = "Unknown"
larynx[larynx$RX_SUMM_RADIATION == 9 ,'RX_SUMM_RADIATION_Recode'] = "Unknown"
larynx[larynx$RX_SUMM_RADIATION == 0 ,'RX_SUMM_RADIATION_Recode'] = "No RT"


# T-stage
larynx$TNM_CLIN_T_Recode  = larynx$TNM_CLIN_T
larynx[
larynx$TNM_CLIN_T ==1    |
larynx$TNM_CLIN_T =='1A' |
larynx$TNM_CLIN_T =='1B' |
larynx$TNM_CLIN_T =='1C' |
larynx$TNM_CLIN_T =='1MI'
,'TNM_CLIN_T_Recode'] = 1

# N-stage
larynx$TNM_CLIN_N_Recode  = larynx$TNM_CLIN_N
larynx[
larynx$TNM_CLIN_N ==1    |
larynx$TNM_CLIN_N =='1A' |
larynx$TNM_CLIN_N =='1B' |
larynx$TNM_CLIN_N =='1C' |
larynx$TNM_CLIN_N =='1MI'
,'TNM_CLIN_N_Recode'] = 1

larynx[
larynx$TNM_CLIN_N ==0    |
larynx$TNM_CLIN_N =='0M+' |
larynx$TNM_CLIN_N =='0M-' |
larynx$TNM_CLIN_N =='0I+' |
larynx$TNM_CLIN_N =='0I-'
,'TNM_CLIN_N_Recode'] = 0

# M-stage
larynx$TNM_CLIN_M_Recode  = larynx$TNM_CLIN_M
larynx[
larynx$TNM_CLIN_M =='X'    |
larynx$TNM_CLIN_M ==' '    |
larynx$TNM_CLIN_M ==88 
,'TNM_CLIN_M_Recode'] = 'X'

# Age
larynx$AGE_Recode=NULL
larynx[larynx$AGE <50 ,'AGE_Recode']                 = "1. <50"
larynx[larynx$AGE <65 & larynx$AGE >=50,'AGE_Recode'] = "2. 50-65"
larynx[larynx$AGE >=65,'AGE_Recode']                 = "3. > or = 65"

# recode radiation sequencing
larynx$RX_SUMM_SURGRAD_SEQ_Recode = NULL
larynx[larynx$RX_SUMM_SURGRAD_SEQ ==0 ,'RX_SUMM_SURGRAD_SEQ_Recode'] = "No RT"
larynx[larynx$RX_SUMM_SURGRAD_SEQ ==2 ,'RX_SUMM_SURGRAD_SEQ_Recode'] = "Neoadjuvant RT"
larynx[larynx$RX_SUMM_SURGRAD_SEQ ==3 ,'RX_SUMM_SURGRAD_SEQ_Recode'] = "Adjuvant RT"
larynx[larynx$RX_SUMM_SURGRAD_SEQ ==4 ,'RX_SUMM_SURGRAD_SEQ_Recode'] = "Before and After Surgery"
larynx[larynx$RX_SUMM_SURGRAD_SEQ ==5 ,'RX_SUMM_SURGRAD_SEQ_Recode'] = "Intra-Operative"
larynx[larynx$RX_SUMM_SURGRAD_SEQ ==6 ,'RX_SUMM_SURGRAD_SEQ_Recode'] = "Intra-operative and additional RT"
larynx[larynx$RX_SUMM_SURGRAD_SEQ ==9 ,'RX_SUMM_SURGRAD_SEQ_Recode'] = "Unknown"

# Race # DO NOT CHANGE ORDER, It matters for White --> Hispanic
larynx$RACE_Recode = '5. Other'
larynx[larynx$RACE ==1 ,'RACE_Recode'] = "1. White"
larynx[(larynx$SPANISH_HISPANIC_ORIGIN >=1 & larynx$SPANISH_HISPANIC_ORIGIN <=6) | (larynx$SPANISH_HISPANIC_ORIGIN ==8),'RACE_Recode'] = "3. Hispanic"
larynx[larynx$RACE ==2 ,'RACE_Recode'] = "2. Black"

larynx[larynx$RACE >=4 & larynx$RACE <=97,'RACE_Recode'] = "4. East Asian, South Asian, Pacific Islander"
# Race
larynx$HISPANIC_Recode = NULL
larynx[larynx$SPANISH_HISPANIC_ORIGIN ==1 ,'HISPANIC_Recode'] = "Hispanic"
larynx[larynx$SPANISH_HISPANIC_ORIGIN ==2 ,'HISPANIC_Recode'] = "Hispanic"
larynx[larynx$SPANISH_HISPANIC_ORIGIN ==3 ,'HISPANIC_Recode'] = "Hispanic"
larynx[larynx$SPANISH_HISPANIC_ORIGIN ==4 ,'HISPANIC_Recode'] = "Hispanic"
larynx[larynx$SPANISH_HISPANIC_ORIGIN ==5 ,'HISPANIC_Recode'] = "Hispanic"
larynx[larynx$SPANISH_HISPANIC_ORIGIN ==6 ,'HISPANIC_Recode'] = "Hispanic"
larynx[larynx$SPANISH_HISPANIC_ORIGIN ==7 ,'HISPANIC_Recode'] = "Other"
larynx[larynx$SPANISH_HISPANIC_ORIGIN ==8 ,'HISPANIC_Recode'] = "Hispanic"
larynx[larynx$SPANISH_HISPANIC_ORIGIN ==9 ,'HISPANIC_Recode'] = "Other"
larynx[larynx$SPANISH_HISPANIC_ORIGIN ==0 ,'HISPANIC_Recode'] = "Non-Hispanic"


# ER Status
larynx$ERSTATUS = '3. Unknown'
larynx[larynx$CS_SITESPECIFIC_FACTOR_1 ==10 ,'ERSTATUS'] = "1. Positive"
larynx[larynx$CS_SITESPECIFIC_FACTOR_1 ==30 ,'ERSTATUS'] = "3. Unknown"
larynx[larynx$CS_SITESPECIFIC_FACTOR_1 ==20 ,'ERSTATUS'] = "2. Negative"

# YEAR OF DIAGNOSIS
larynx$YEAR_OF_DIAGNOSIS_Recode = "2008+"
larynx[larynx$YEAR_OF_DIAGNOSIS <2008 ,'YEAR_OF_DIAGNOSIS_Recode'] = "Before 2008"


# Grade
larynx$GRADE_Recode = NULL
larynx[larynx$GRADE ==1 ,'GRADE_Recode'] = "1. Well Differentiated"
larynx[larynx$GRADE ==2 ,'GRADE_Recode'] = "2. Moderately Differentiated"
larynx[larynx$GRADE ==3 ,'GRADE_Recode'] = "3. Poorly Differentiated"
larynx[larynx$GRADE ==4 ,'GRADE_Recode'] = "3. Poorly Differentiated"
larynx[larynx$GRADE ==9 ,'GRADE_Recode'] = "4. Unknown"

# LVI
larynx$LYMPH_VASCULAR_INVASION_Recode = NULL
larynx[larynx$LYMPH_VASCULAR_INVASION ==0 & !is.na(larynx$LYMPH_VASCULAR_INVASION) ,'LYMPH_VASCULAR_INVASION_Recode'] = "Absent"
larynx[larynx$LYMPH_VASCULAR_INVASION ==1 & !is.na(larynx$LYMPH_VASCULAR_INVASION) ,'LYMPH_VASCULAR_INVASION_Recode'] = "Present"
larynx[larynx$LYMPH_VASCULAR_INVASION ==8 & !is.na(larynx$LYMPH_VASCULAR_INVASION) ,'LYMPH_VASCULAR_INVASION_Recode'] = "Unknown/Indeterminate"
larynx[larynx$LYMPH_VASCULAR_INVASION ==9 & !is.na(larynx$LYMPH_VASCULAR_INVASION) ,'LYMPH_VASCULAR_INVASION_Recode'] = "Unknown/Indeterminate"
larynx[is.na(larynx$LYMPH_VASCULAR_INVASION) ,'LYMPH_VASCULAR_INVASION_Recode'] = "Unknown/Indeterminate"


# Charlson-Deyo
larynx$CDCC_TOTAL_Recode = NULL
larynx[larynx$CDCC_TOTAL ==0 ,'CDCC_TOTAL_Recode'] = "0"
larynx[larynx$CDCC_TOTAL ==1 ,'CDCC_TOTAL_Recode'] = "1"
larynx[larynx$CDCC_TOTAL ==2 ,'CDCC_TOTAL_Recode'] = "2+"


# RX SUMMARY CHEMO
larynx$RX_SUMM_CHEMO_Recode = NULL
larynx[larynx$RX_SUMM_CHEMO ==0 ,'RX_SUMM_CHEMO_Recode'] = "No Chemotherapy"
larynx[larynx$RX_SUMM_CHEMO ==1 ,'RX_SUMM_CHEMO_Recode'] = "Chemotherapy"
larynx[larynx$RX_SUMM_CHEMO ==2 ,'RX_SUMM_CHEMO_Recode'] = "Chemotherapy"
larynx[larynx$RX_SUMM_CHEMO ==3 ,'RX_SUMM_CHEMO_Recode'] = "Chemotherapy"
larynx[larynx$RX_SUMM_CHEMO ==82 ,'RX_SUMM_CHEMO_Recode'] = "No Chemotherapy"
larynx[larynx$RX_SUMM_CHEMO ==85 ,'RX_SUMM_CHEMO_Recode'] = "No Chemotherapy"
larynx[larynx$RX_SUMM_CHEMO ==86 ,'RX_SUMM_CHEMO_Recode'] = "No Chemotherapy"
larynx[larynx$RX_SUMM_CHEMO ==87 ,'RX_SUMM_CHEMO_Recode'] = "No Chemotherapy"
larynx[larynx$RX_SUMM_CHEMO ==88 ,'RX_SUMM_CHEMO_Recode'] = "Unknown"
larynx[larynx$RX_SUMM_CHEMO ==99 ,'RX_SUMM_CHEMO_Recode'] = "Unknown"


# RX SUMM SYSTEMIC SUR SEQ 
larynx$RX_SUMM_SYSTEMIC_SUR_SEQ_Recode = NULL
larynx[larynx$RX_SUMM_SYSTEMIC_SUR_SEQ ==0 & !is.na(larynx$RX_SUMM_SYSTEMIC_SUR_SEQ),'RX_SUMM_SYSTEMIC_SUR_SEQ_Recode'] = "No systemic therapy"
larynx[larynx$RX_SUMM_SYSTEMIC_SUR_SEQ ==2 & !is.na(larynx$RX_SUMM_SYSTEMIC_SUR_SEQ),'RX_SUMM_SYSTEMIC_SUR_SEQ_Recode'] = "Neoadjuvant"
larynx[larynx$RX_SUMM_SYSTEMIC_SUR_SEQ ==3 & !is.na(larynx$RX_SUMM_SYSTEMIC_SUR_SEQ),'RX_SUMM_SYSTEMIC_SUR_SEQ_Recode'] = "Adjuvant"
larynx[larynx$RX_SUMM_SYSTEMIC_SUR_SEQ ==4 & !is.na(larynx$RX_SUMM_SYSTEMIC_SUR_SEQ),'RX_SUMM_SYSTEMIC_SUR_SEQ_Recode'] = "Neoadjuvant"
larynx[larynx$RX_SUMM_SYSTEMIC_SUR_SEQ ==5 & !is.na(larynx$RX_SUMM_SYSTEMIC_SUR_SEQ),'RX_SUMM_SYSTEMIC_SUR_SEQ_Recode'] = "Unknown"
larynx[larynx$RX_SUMM_SYSTEMIC_SUR_SEQ ==6 & !is.na(larynx$RX_SUMM_SYSTEMIC_SUR_SEQ),'RX_SUMM_SYSTEMIC_SUR_SEQ_Recode'] = "Unknown"
larynx[larynx$RX_SUMM_SYSTEMIC_SUR_SEQ ==9 & !is.na(larynx$RX_SUMM_SYSTEMIC_SUR_SEQ),'RX_SUMM_SYSTEMIC_SUR_SEQ_Recode'] = "Unknown"
larynx[is.na(larynx$RX_SUMM_SYSTEMIC_SUR_SEQ),'RX_SUMM_SYSTEMIC_SUR_SEQ_Recode'] = "Unknown"

# RX_SUMM_HORMONE
larynx$RX_SUMM_HORMONE_Recode = NULL
larynx[larynx$RX_SUMM_HORMONE ==0 ,'RX_SUMM_HORMONE_Recode'] = "None"
larynx[larynx$RX_SUMM_HORMONE ==1 ,'RX_SUMM_HORMONE_Recode'] = "Hormonal Therapy"
larynx[larynx$RX_SUMM_HORMONE ==82 ,'RX_SUMM_HORMONE_Recode'] = "None"
larynx[larynx$RX_SUMM_HORMONE ==85 ,'RX_SUMM_HORMONE_Recode'] = "None"
larynx[larynx$RX_SUMM_HORMONE ==86 ,'RX_SUMM_HORMONE_Recode'] = "None"
larynx[larynx$RX_SUMM_HORMONE ==87 ,'RX_SUMM_HORMONE_Recode'] = "None"
larynx[larynx$RX_SUMM_HORMONE ==88 ,'RX_SUMM_HORMONE_Recode'] = "Unknown"
larynx[larynx$RX_SUMM_HORMONE ==99 ,'RX_SUMM_HORMONE_Recode'] = "Unknown"


# RX_SUMM_SURGICAL_MARGINS
larynx$RX_SUMM_SURGICAL_MARGINS_Recode = NULL
larynx[larynx$RX_SUMM_SURGICAL_MARGINS ==0 ,'RX_SUMM_SURGICAL_MARGINS_Recode'] = "Negative"
larynx[larynx$RX_SUMM_SURGICAL_MARGINS ==1 ,'RX_SUMM_SURGICAL_MARGINS_Recode'] = "Positive Macroscopic/Microscopic"
larynx[larynx$RX_SUMM_SURGICAL_MARGINS ==2 ,'RX_SUMM_SURGICAL_MARGINS_Recode'] = "Positive Macroscopic/Microscopic"
larynx[larynx$RX_SUMM_SURGICAL_MARGINS ==3 ,'RX_SUMM_SURGICAL_MARGINS_Recode'] = "Positive Macroscopic/Microscopic"
larynx[larynx$RX_SUMM_SURGICAL_MARGINS ==7 ,'RX_SUMM_SURGICAL_MARGINS_Recode'] = "Unknown"
larynx[larynx$RX_SUMM_SURGICAL_MARGINS ==9 ,'RX_SUMM_SURGICAL_MARGINS_Recode'] = "Unknown"
larynx[is.na(larynx$RX_SUMM_SURGICAL_MARGINS),'RX_SUMM_SURGICAL_MARGINS_Recode'] = "Unknown"


#recode facility
# FACILITY_TYPE_CD_Recode
larynx$FACILITY_TYPE_CD_Recode = NULL
larynx[larynx$FACILITY_TYPE_CD ==1 & !is.na(larynx$FACILITY_TYPE_CD),'FACILITY_TYPE_CD_Recode'] = "Community Cancer Program"
larynx[larynx$FACILITY_TYPE_CD ==2 & !is.na(larynx$FACILITY_TYPE_CD),'FACILITY_TYPE_CD_Recode'] = "Comprehensive Community Cancer Program"
larynx[larynx$FACILITY_TYPE_CD ==3 & !is.na(larynx$FACILITY_TYPE_CD),'FACILITY_TYPE_CD_Recode'] = "Academic/Research Program"
larynx[larynx$FACILITY_TYPE_CD ==4 & !is.na(larynx$FACILITY_TYPE_CD),'FACILITY_TYPE_CD_Recode'] = "Integrated Network Cancer Program"
larynx[larynx$FACILITY_TYPE_CD ==9 & !is.na(larynx$FACILITY_TYPE_CD),'FACILITY_TYPE_CD_Recode'] = "Community Cancer Program"
larynx[is.na(larynx$FACILITY_TYPE_CD),'FACILITY_TYPE_CD_Recode'] = "Community Cancer Program"


#recode facility type
# FACILITY_LOCATION_CD_Recode
larynx$FACILITY_LOCATION_CD_Recode = NULL
larynx[larynx$FACILITY_LOCATION_CD ==1 & !is.na(larynx$FACILITY_LOCATION_CD),'FACILITY_LOCATION_CD_Recode'] = "New England"
larynx[larynx$FACILITY_LOCATION_CD ==2 & !is.na(larynx$FACILITY_LOCATION_CD),'FACILITY_LOCATION_CD_Recode'] = "Middle Atlantic"
larynx[larynx$FACILITY_LOCATION_CD ==3 & !is.na(larynx$FACILITY_LOCATION_CD),'FACILITY_LOCATION_CD_Recode'] = "South Atlantic"
larynx[larynx$FACILITY_LOCATION_CD ==4 & !is.na(larynx$FACILITY_LOCATION_CD),'FACILITY_LOCATION_CD_Recode'] = "East North Central"
larynx[larynx$FACILITY_LOCATION_CD ==5 & !is.na(larynx$FACILITY_LOCATION_CD),'FACILITY_LOCATION_CD_Recode'] = "East South Central"
larynx[larynx$FACILITY_LOCATION_CD ==6 & !is.na(larynx$FACILITY_LOCATION_CD),'FACILITY_LOCATION_CD_Recode'] = "West North Central"
larynx[larynx$FACILITY_LOCATION_CD ==7 & !is.na(larynx$FACILITY_LOCATION_CD),'FACILITY_LOCATION_CD_Recode'] = "West South Central"
larynx[larynx$FACILITY_LOCATION_CD ==8 & !is.na(larynx$FACILITY_LOCATION_CD),'FACILITY_LOCATION_CD_Recode'] = "Mountain"
larynx[larynx$FACILITY_LOCATION_CD ==9 & !is.na(larynx$FACILITY_LOCATION_CD),'FACILITY_LOCATION_CD_Recode'] = "Pacific"
larynx[is.na(larynx$FACILITY_LOCATION_CD),'FACILITY_LOCATION_CD_Recode'] = "Unknown"

larynx$HISTOLOGY_Recode = 'Other'
larynx[larynx$HISTOLOGY >= 8500 &  larynx$HISTOLOGY <= 8500 & !is.na(larynx$HISTOLOGY),'HISTOLOGY_Recode'] = "Ductal"
larynx[larynx$HISTOLOGY >= 8520 &  larynx$HISTOLOGY <= 8520 & !is.na(larynx$HISTOLOGY),'HISTOLOGY_Recode'] = "Lobular"
larynx[larynx$HISTOLOGY >= 8522 &  larynx$HISTOLOGY <= 8523 & !is.na(larynx$HISTOLOGY),'HISTOLOGY_Recode'] = "Mixed Ductal and Lobular"


# CROWFLY
larynx$CROWFLY_Recode = '1. <10'
larynx[larynx$CROWFLY >= 10 &  larynx$CROWFLY <= 20 & !is.na(larynx$CROWFLY),'CROWFLY_Recode'] = "2. 10-20"
larynx[larynx$CROWFLY >= 20 &  larynx$CROWFLY <= 50 & !is.na(larynx$CROWFLY),'CROWFLY_Recode'] = "3. 20-50"
larynx[larynx$CROWFLY >= 50 &  larynx$CROWFLY <= 100 & !is.na(larynx$CROWFLY),'CROWFLY_Recode'] = "4. 50-100"
larynx[larynx$CROWFLY >= 100 & !is.na(larynx$CROWFLY),'CROWFLY_Recode'] = "5. >100"


# Education
larynx$EDUCATION_Recode = NULL
larynx[larynx$NO_HSD_QUAR_12 == 1 & !is.na(larynx$NO_HSD_QUAR_12),'EDUCATION_Recode'] = "1. >21%"
larynx[larynx$NO_HSD_QUAR_12 == 2 & !is.na(larynx$NO_HSD_QUAR_12),'EDUCATION_Recode'] = "2. 13-21%"
larynx[larynx$NO_HSD_QUAR_12 == 3 & !is.na(larynx$NO_HSD_QUAR_12),'EDUCATION_Recode'] = "3. 7-12.9%"
larynx[larynx$NO_HSD_QUAR_12 == 4 & !is.na(larynx$NO_HSD_QUAR_12),'EDUCATION_Recode'] = "4. <7%"
larynx[is.na(larynx$NO_HSD_QUAR_12),'EDUCATION_Recode'] = "Unknown"


# Education2000
larynx$EDUCATION2000_Recode = NULL
larynx[larynx$NO_HSD_QUAR_00 == 1 & !is.na(larynx$NO_HSD_QUAR_00),'EDUCATION2000_Recode'] = "1. >21%"
larynx[larynx$NO_HSD_QUAR_00 == 2 & !is.na(larynx$NO_HSD_QUAR_00),'EDUCATION2000_Recode'] = "2. 13-21%"
larynx[larynx$NO_HSD_QUAR_00 == 3 & !is.na(larynx$NO_HSD_QUAR_00),'EDUCATION2000_Recode'] = "3. 7-12.9%"
larynx[larynx$NO_HSD_QUAR_00 == 4 & !is.na(larynx$NO_HSD_QUAR_00),'EDUCATION2000_Recode'] = "4. <7%"
larynx[is.na(larynx$NO_HSD_QUAR_00),'EDUCATION2000_Recode'] = "Unknown"


# Income
larynx$INCOME_Recode = NULL
larynx[larynx$MED_INC_QUAR_12 ==1 & !is.na(larynx$MED_INC_QUAR_12),'INCOME_Recode'] = "1. <38K"
larynx[larynx$MED_INC_QUAR_12 == 2 & !is.na(larynx$MED_INC_QUAR_12),'INCOME_Recode'] = "2. 38K-48K"
larynx[larynx$MED_INC_QUAR_12 == 3 & !is.na(larynx$MED_INC_QUAR_12),'INCOME_Recode'] = "4. 48K-63K"
larynx[larynx$MED_INC_QUAR_12 == 4 & !is.na(larynx$MED_INC_QUAR_12),'INCOME_Recode'] = "5. >63K"
larynx[is.na(larynx$MED_INC_QUAR_12),'INCOME_Recode'] = "Unknown"

# Income2000
larynx$INCOME2000_Recode = NULL
larynx[larynx$MED_INC_QUAR_00 ==1 & !is.na(larynx$MED_INC_QUAR_00),'INCOME2000_Recode'] = "1. <30K"
larynx[larynx$MED_INC_QUAR_00 == 2 & !is.na(larynx$MED_INC_QUAR_00),'INCOME2000_Recode'] = "2. 30K-35K"
larynx[larynx$MED_INC_QUAR_00 == 3 & !is.na(larynx$MED_INC_QUAR_00),'INCOME2000_Recode'] = "4. 35K-46K"
larynx[larynx$MED_INC_QUAR_00 == 4 & !is.na(larynx$MED_INC_QUAR_00),'INCOME2000_Recode'] = "5. >46K"
larynx[is.na(larynx$MED_INC_QUAR_00),'INCOME2000_Recode'] = "Unknown"

# Insurance
larynx$INSURANCE_STATUS_Recode = NULL
larynx[larynx$INSURANCE_STATUS ==1,'INSURANCE_STATUS_Recode'] = "1. Private"
larynx[larynx$INSURANCE_STATUS == 2 ,'INSURANCE_STATUS_Recode'] = "2. Medicaid"
larynx[larynx$INSURANCE_STATUS == 3 ,'INSURANCE_STATUS_Recode'] = "3. Medicare"
larynx[larynx$INSURANCE_STATUS == 4 ,'INSURANCE_STATUS_Recode'] = "4. Other Gvt"
larynx[larynx$INSURANCE_STATUS == 0 ,'INSURANCE_STATUS_Recode'] = "5. Not Insured"
larynx[larynx$INSURANCE_STATUS == 9 ,'INSURANCE_STATUS_Recode'] = "6. Unknown"

#survival
larynx$followup_time = larynx["DX_LASTCONTACT_DEATH_MONTHS"][[1]]
censor  = larynx["PUF_VITAL_STATUS"][[1]]
larynx$censor = 1-censor 

#############################

# select variables
myvars <- c('SEX','IMRT_Bin','IMRT',"AGE_Recode","YEAR_OF_DIAGNOSIS_Recode","TNM_CLIN_T_Recode","TNM_CLIN_N_Recode","RACE_Recode","GRADE_Recode",'CDCC_TOTAL_Recode','RX_SUMM_CHEMO_Recode','FACILITY_TYPE_CD_Recode','CROWFLY_Recode','followup_time','censor','HISPANIC_Recode','EDUCATION_Recode','INCOME_Recode','INSURANCE_STATUS_Recode')
matchdata <- larynx[myvars]

# Create the match based on stage, age, chemo using the nearest neighbor method etc...)
m.out=matchit(IMRT_Bin ~ SEX 
+ as.factor(AGE_Recode) 
+ as.factor(RACE_Recode) 
+ as.factor(CDCC_TOTAL_Recode) 
+ as.factor(TNM_CLIN_T_Recode)
+ as.factor(EDUCATION_Recode) 
+ as.factor(INSURANCE_STATUS_Recode)
+ as.factor(FACILITY_TYPE_CD_Recode)
, data = matchdata, method = "nearest");
#summary(m.out);

larynx_Matched <- match.data(m.out)

# descriptive statistics
listVars<-(c('SEX','IMRT_Bin','IMRT',"AGE_Recode","YEAR_OF_DIAGNOSIS_Recode","TNM_CLIN_T_Recode","TNM_CLIN_N_Recode","RACE_Recode","GRADE_Recode",'CDCC_TOTAL_Recode','RX_SUMM_CHEMO_Recode','FACILITY_TYPE_CD_Recode','CROWFLY_Recode','EDUCATION_Recode','INCOME_Recode','INSURANCE_STATUS_Recode'))

catVars<-(c('SEX','IMRT_Bin','IMRT',"AGE_Recode","YEAR_OF_DIAGNOSIS_Recode","TNM_CLIN_T_Recode","TNM_CLIN_N_Recode","RACE_Recode","GRADE_Recode",'CDCC_TOTAL_Recode','RX_SUMM_CHEMO_Recode','FACILITY_TYPE_CD_Recode','CROWFLY_Recode','EDUCATION_Recode','INCOME_Recode','INSURANCE_STATUS_Recode'))

table1 <- CreateTableOne(vars = listVars, data = larynx_Matched, factorVars = catVars, strata = 'IMRT',includeNA=TRUE)

 #############################
write.csv(print(table1),'H:\\NCDB Files\\larynx_IMRT3D_matched.csv')

larynx_test <-larynx_Matched

# create survival object
allpts <- survfit(Surv(followup_time, censor)~ 1,data=larynx_test)
plot(allpts, xlab="Time", ylab="Survival Probability (%)",lwd=2, main="Kaplan Meier Plot for Survival of All T1-2N0 with RT Alone")

x11()

# stratify by IMRT

#survsum(Surv(followup_time, censor)~IMRT,sptmr=c(60),data=larynx_test)
subset <- survfit(Surv(followup_time, censor)~IMRT,data=larynx_test)
plot(subset, xlab="Months", ylab="Survival Probability (%)", main='Overall Survival of All T1-2N0 Supraglottic Larynx \n Treated with Definitive RT Stratified by IMRT',lwd=2,col=c(1:length(names(subset$strata))),lty=1,cex.axis = 2,cex.lab = 1.5,cex.main = 1.5)
legend(80, 1.0, c("IMRT","Non-IMRT"),lwd=4,col=c(1:length(names(subset$strata))),lty=1)
logrank <- survdiff(Surv(followup_time, censor)~IMRT,data=larynx_test)
pval <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
text(30,0.07,paste("logrank ","p =",round(pval,4),sep = ' '))



#SURVIVAL ANALYSIS #######
UNI_Results <- NULL


coxmodel <- coxph(Surv(followup_time, censor)~ as.factor(IMRT),data=larynx_test)
UNI_Results = rbind(UNI_Results,cbind(exp(coxmodel$coefficients),exp(confint(coxmodel)),round(coef(summary(coxmodel))[,5],2)))

coxmodel <- coxph(Surv(followup_time, censor)~ as.factor(SEX),data=larynx_test)
UNI_Results = rbind(UNI_Results,cbind(exp(coxmodel$coefficients),exp(confint(coxmodel)),round(coef(summary(coxmodel))[,5],2)))


coxmodel <- coxph(Surv(followup_time, censor)~ as.factor(AGE_Recode),data=larynx_test)
UNI_Results = rbind(UNI_Results,cbind(exp(coxmodel$coefficients),exp(confint(coxmodel)),round(coef(summary(coxmodel))[,5],2)))

coxmodel <- coxph(Surv(followup_time, censor)~ as.factor(YEAR_OF_DIAGNOSIS_Recode),data=larynx_test)
UNI_Results = rbind(UNI_Results,cbind(exp(coxmodel$coefficients),exp(confint(coxmodel)),round(coef(summary(coxmodel))[,5],2)))

coxmodel <- coxph(Surv(followup_time, censor)~ as.factor(RACE_Recode),data=larynx_test)
UNI_Results = rbind(UNI_Results,cbind(exp(coxmodel$coefficients),exp(confint(coxmodel)),round(coef(summary(coxmodel))[,5],2)))

coxmodel <- coxph(Surv(followup_time, censor)~ as.factor(droplevels(TNM_CLIN_N_Recode)),data=larynx_test)
UNI_Results = rbind(UNI_Results,cbind(exp(coxmodel$coefficients),exp(confint(coxmodel)),round(coef(summary(coxmodel))[,5],2)))

coxmodel <- coxph(Surv(followup_time, censor)~ as.factor(droplevels(TNM_CLIN_T_Recode)),data=larynx_test)
UNI_Results = rbind(UNI_Results,cbind(exp(coxmodel$coefficients),exp(confint(coxmodel)),round(coef(summary(coxmodel))[,5],2)))

coxmodel <- coxph(Surv(followup_time, censor)~ as.factor(GRADE_Recode),data=larynx_test)
UNI_Results = rbind(UNI_Results,cbind(exp(coxmodel$coefficients),exp(confint(coxmodel)),round(coef(summary(coxmodel))[,5],2)))

coxmodel <- coxph(Surv(followup_time, censor)~ as.factor(CDCC_TOTAL_Recode),data=larynx_test)
UNI_Results = rbind(UNI_Results,cbind(exp(coxmodel$coefficients),exp(confint(coxmodel)),round(coef(summary(coxmodel))[,5],2)))

coxmodel <- coxph(Surv(followup_time, censor)~ as.factor(FACILITY_TYPE_CD_Recode),data=larynx_test)
UNI_Results = rbind(UNI_Results,cbind(exp(coxmodel$coefficients),exp(confint(coxmodel)),round(coef(summary(coxmodel))[,5],2)))


coxmodel <- coxph(Surv(followup_time, censor)~ INCOME_Recode,data=larynx_test)
UNI_Results = rbind(UNI_Results,cbind(exp(coxmodel$coefficients),exp(confint(coxmodel)),round(coef(summary(coxmodel))[,5],2)))


coxmodel <- coxph(Surv(followup_time, censor)~ as.factor(INSURANCE_STATUS_Recode),data=larynx_test)
UNI_Results = rbind(UNI_Results,cbind(exp(coxmodel$coefficients),exp(confint(coxmodel)),round(coef(summary(coxmodel))[,5],2)))

coxmodel <- coxph(Surv(followup_time, censor)~ as.factor(CROWFLY_Recode),data=larynx_test)
UNI_Results = rbind(UNI_Results,cbind(exp(coxmodel$coefficients),exp(confint(coxmodel)),round(coef(summary(coxmodel))[,5],2)))

coxmodel <- coxph(Surv(followup_time, censor)~ EDUCATION_Recode,data=larynx_test)
UNI_Results = rbind(UNI_Results,cbind(exp(coxmodel$coefficients),exp(confint(coxmodel)),round(coef(summary(coxmodel))[,5],2)))

UNI_Results<-round(UNI_Results,3)
colnames(UNI_Results) = c('HR','2.5% CI','97.5% CI','p')
UNI_Results
write.csv(UNI_Results,"H:/NCDB Files/UNI_matched_T1-2N0Larynx.csv")
########

############# STEPWISE SELECTION #############

('SEX','IMRT_Bin','IMRT',"AGE_Recode","YEAR_OF_DIAGNOSIS","TNM_CLIN_T_Recode","TNM_CLIN_N_Recode","RACE_Recode","GRADE_Recode",'CDCC_TOTAL_Recode','FACILITY_TYPE_CD_Recode','CROWFLY_Recode','EDUCATION_Recode','INCOME_Recode','INSURANCE_STATUS_Recode')

full <- coxph(Surv(followup_time, censor)~ 
+ as.factor(IMRT) 
+ as.factor(SEX) 
+ as.factor(AGE_Recode) 
+ as.factor(CDCC_TOTAL_Recode)
+ as.factor(droplevels(TNM_CLIN_T_Recode))
+ as.factor(YEAR_OF_DIAGNOSIS_Recode)
+ as.factor(GRADE_Recode)
+ as.factor(RACE_Recode)
+ as.factor(INCOME_Recode)
+ as.factor(INSURANCE_STATUS_Recode)
+ as.factor(CROWFLY_Recode)
+ as.factor(FACILITY_TYPE_CD_Recode)
+ as.factor(EDUCATION_Recode)
,data=larynx_test)



forest_model(coxmodel)

null = coxph(Surv(followup_time, censor)~ as.factor(IMRT)) 



stepAIC(null,direction='forward',scope=list(lower=null,upper=full))
stepAIC(full,direction='backward',scope=list(lower=null,upper=full))

table(larynx$TNM_CLIN_T,larynx$TNM_CLIN_N)

############ Multivariate Survival ###############

coxmodel <- coxph(Surv(followup_time, censor)~ 
+ as.factor(IMRT) 
+ as.factor(SEX) 
+ as.factor(AGE_Recode) 
+ as.factor(CDCC_TOTAL_Recode)
+ as.factor(droplevels(TNM_CLIN_T_Recode))
+ as.factor(YEAR_OF_DIAGNOSIS_Recode)
+ as.factor(GRADE_Recode)
+ as.factor(INCOME_Recode)
+ as.factor(INSURANCE_STATUS_Recode)
+ as.factor(FACILITY_TYPE_CD_Recode)
,data=larynx_test)


forest_model(coxmodel)

#########################################
# propensity score matched cohort analysis

# KM survival analysis

# Cox Model

