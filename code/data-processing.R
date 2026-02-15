library('readxl')
dat <- read_excel("../data/Data extracted.xlsx", sheet = 'Data extracted')
dat <- data.frame(dat)
my_cols <- dat[2, ]

dat <- dat[-c(1:2), ]
colnames(dat) <- my_cols
dat <- dat[!is.na(dat$Author), ]

################################################################################
## Renaming columns
################################################################################
colnames(dat)[colnames(dat) == 'Author'] <- 'Authors_noID'
colnames(dat)[colnames(dat) == 'Study ID'] <- 'StudyID'
colnames(dat)[colnames(dat) == 'Testing site'] <- 'testingsite'
colnames(dat)[colnames(dat) == 'Target population'] <- 'population'
colnames(dat)[colnames(dat) == 'Location'] <- 'country'
colnames(dat)[colnames(dat) == 'Ag-RDT name'] <- 'AgName'
colnames(dat)[colnames(dat) == 'Frequency of Testing'] <- 'Testregime'
colnames(dat)[colnames(dat) == 'Times per week'] <- 'Testfrequency'
colnames(dat)[colnames(dat) == 'Overview symptom status'] <- 'Sympstat'

colnames(dat)[colnames(dat) == 'Number of people that collected / were sent\r\nself-tests'] <- 'UptakeN'
colnames(dat)[colnames(dat) == 'Number of people eligible / invited for self-testing'] <- 'UptakeD'

colnames(dat)[colnames(dat) == 'Number of people that reported self-testing results'] <- 'ReportN'
dat$ReportD <- dat$UptakeN

colnames(dat)[colnames(dat) == 'Number of results reported'] <- 'AdherentN'
colnames(dat)[colnames(dat) == 'Optimal # of tests performed as per testing scheme'] <- 'AdherentD'

colnames(dat)[colnames(dat) == 'Ag-RDT_n people_Positive'] <- 'PpositiveN'
dat$PpositiveD <- dat$ReportN

colnames(dat)[colnames(dat) == 'Ag-RDT_n tests_Positive'] <- 'TpositiveN'
dat$TpositiveD <- dat$AdherentN

colnames(dat)[colnames(dat) == 'Confirm_n tests_FP'] <- 'FPRN'
colnames(dat)[colnames(dat) == 'Confirm_n_tests_FP_ref'] <- 'FPRD'

colnames(dat)[colnames(dat) == 'Confirm_n people_FN'] <- 'MissedN'
colnames(dat)[colnames(dat) == 'Confirm_n people_FN_ref'] <- 'MissedD'

colnames(dat)[colnames(dat) == 'Ag-RDT_n tests_Invalid'] <- 'InvalidN'
dat$InvalidD <- dat$AdherentN

dat[is.na(dat$`SARS-CoV-2 status`) | dat$`SARS-CoV-2 status` != 'unknown', 
    c('PpositiveN', 'PpositiveD', 'TpositiveN', 'TpositiveD')] <- NA 

################################################################################
## Forming a small data set with only the necessary columns
################################################################################
dat_reduced <- dat[, c('Authors_noID', 'StudyID', 'testingsite', 'population', 'country', 
                       'AgName', 'Testregime', 'Testfrequency', 'Sympstat', 
                       'UptakeN', 'UptakeD', 'ReportN', 'ReportD', 
                       'AdherentN', 'AdherentD', 'PpositiveN', 'PpositiveD', 
                       'TpositiveN', 'TpositiveD', 'FPRN', 'FPRD', 
                       'MissedN', 'MissedD', 'InvalidN', 'InvalidD', 
                       'SA_Uptake', 'SA_Reporting', 
                       'FP_analysis', 'Missed_analysis')]

################################################################################
## Making columns numeric
################################################################################
my_numeric_cols <- c('UptakeN', 'UptakeD', 'ReportN', 'ReportD', 
                     'AdherentN', 'AdherentD', 'PpositiveN', 'PpositiveD', 
                     'TpositiveN', 'TpositiveD', 'FPRN', 'FPRD', 
                     'MissedN', 'MissedD', 'InvalidN', 'InvalidD')
for (col in my_numeric_cols){
  dat_reduced[, eval(col)] <- as.numeric(dat_reduced[, eval(col)])
}

################################################################################
## Fixing ReportN and ReportD issue with a80_1/a80_2
################################################################################
myrow <- dat_reduced[dat_reduced$StudyID == 'a79_1', ]
myrow[, my_numeric_cols] <- NA
myind <- dat_reduced$StudyID %in% c('a79_1', 'a79_2')
ReportN_val <- sum(dat_reduced[myind, 'ReportN'], na.rm = T)
ReportD_val <- sum(dat_reduced[myind, 'ReportD'], na.rm = T)
myrow[, c('StudyID', 'ReportN', 'ReportD')] <- c('a79_0', ReportN_val, ReportD_val)
dat_reduced[myind, c('ReportN', 'ReportD')] <- c(NA, NA)
dat_reduced <- rbind(dat_reduced, myrow)

################################################################################
## Making entries NA if we don't have numeric values for both _N and _D column types
## of if other restrictions apply (eg, FPR and Missed)
################################################################################
na_Uptake <- is.na(dat_reduced$UptakeN) | is.na(dat_reduced$UptakeD)
na_Report <- is.na(dat_reduced$ReportN) | is.na(dat_reduced$ReportD)
na_Adherent <- is.na(dat_reduced$AdherentN) | is.na(dat_reduced$AdherentD)
na_Ppositive <- is.na(dat_reduced$PpositiveN) | is.na(dat_reduced$PpositiveD)
na_Tpositive <- is.na(dat_reduced$TpositiveN) | is.na(dat_reduced$TpositiveD)
na_FPR <- is.na(dat_reduced$FPRN) | is.na(dat_reduced$FPRD) | dat_reduced$FP_analysis != 'include'
na_Missed <- is.na(dat_reduced$MissedN) | is.na(dat_reduced$MissedD) | dat_reduced$Missed_analysis != 'include'
na_Invalid <- is.na(dat_reduced$InvalidN) | is.na(dat_reduced$InvalidD)

dat_reduced[na_Uptake, c('UptakeN', 'UptakeD')] <- NA
dat_reduced[na_Report, c('ReportN', 'ReportD')] <- NA
dat_reduced[na_Adherent, c('AdherentN', 'AdherentD')] <- NA
dat_reduced[na_Ppositive, c('PpositiveN', 'PpositiveD')] <- NA
dat_reduced[na_Tpositive, c('TpositiveN', 'TpositiveD')] <- NA
dat_reduced[na_FPR, c('FPRN', 'FPRD')] <- NA
dat_reduced[na_Missed, c('MissedN', 'MissedD')] <- NA
dat_reduced[na_Invalid, c('InvalidN', 'InvalidD')] <- NA

all_na <- function(x) {return(all(is.na(x)))}
all_na_cols <- apply(dat_reduced[, my_numeric_cols], 1, all_na)
dat_reduced <- dat_reduced[!all_na_cols, ]

################################################################################
## Checking for the presence of 0/0 datasets and removing them if applicable
################################################################################
vars <- c("Uptake", "Report", "Adherent", "Ppositive", "Tpositive",
          "FPR", "Missed", "Invalid")
for (v in vars) {
  num_col <- paste0(v, "N"); den_col <- paste0(v, "D")
  zero_idx <- which(dat_reduced[[num_col]] == 0 & dat_reduced[[den_col]] == 0)
  if (length(zero_idx) > 0) {
    dat_reduced[zero_idx, num_col] <- NA
    dat_reduced[zero_idx, den_col] <- NA
    message(sprintf("Set %sN and %sD to NA in %d row(s)", v, v, length(zero_idx)))
  }
}

################################################################################
## Creating testing frequency variable
################################################################################

dat_reduced$Testfrequency[dat_reduced$Testfrequency == '#'] <- NA
dat_reduced$Testfrequency[dat_reduced$Testfrequency == 1] <- 'Once'
dat_reduced$Testfrequency[dat_reduced$Testfrequency == 2] <- 'Twice weekly'
dat_reduced$Testfrequency[!is.na(dat_reduced$Testfrequency)
                          & !dat_reduced$Testfrequency %in% c('Once', 'Twice weekly', 'Unclear')] <- 'Three or more times weekly'

dat_reduced$Testregime_new <- 
  ifelse(dat_reduced$Testregime == 'Single', 'Single', 
         ifelse(dat_reduced$Testregime == 'Unclear', 'Unclear', 
                ifelse(dat_reduced$Testregime == 'Routine' & dat_reduced$Testfrequency == 'Once', 'Routine (once weekly)', 
                       ifelse(dat_reduced$Testregime == 'Routine' & dat_reduced$Testfrequency == 'Twice weekly', 'Routine (twice weekly)', 
                              ifelse(dat_reduced$Testregime == 'Routine' & dat_reduced$Testfrequency == 'Three or more times weekly', 'Routine (three or more times weekly)', 
                                     ifelse(dat_reduced$Testregime == 'Routine' & dat_reduced$Testfrequency == 'Unclear', 'Routine (unclear frequency)', 'ERROR'))))))

################################################################################
## Creating Authors name variable
################################################################################

Authors_noperiod <- sub("\\.$", "", dat_reduced$Authors_noID)
StudyID_appendix <- sub(".*(_.*)$", "\\1", dat_reduced$StudyID)
dat_reduced$Authors <- paste0(Authors_noperiod, StudyID_appendix)

################################################################################
## Cleaning up symptom status variable
################################################################################

dat_reduced$Sympstat <- ifelse(dat_reduced$Sympstat == '1', 'Asymptomatic', 'Mixed/Symptomatic')

################################################################################
## Creating the main and sensitivity uptake and reporting variables
################################################################################

# Uptake
main_studies_uptake <- dat_reduced$SA_Uptake != 'include'
dat_reduced$UptakeMainN <- dat_reduced$UptakeMainD <- NA
dat_reduced[main_studies_uptake, ]$UptakeMainN <- dat_reduced[main_studies_uptake, ]$UptakeN
dat_reduced[main_studies_uptake, ]$UptakeMainD <- dat_reduced[main_studies_uptake, ]$UptakeD
dat_reduced$UptakeSensN <- dat_reduced$UptakeN
dat_reduced$UptakeSensD <- dat_reduced$UptakeD
dat_reduced$UptakeN <- dat_reduced$UptakeD <- NULL

# Reporting
main_studies_report <- dat_reduced$SA_Reporting != 'include'
dat_reduced$ReportMainN <- dat_reduced$ReportMainD <- NA
dat_reduced[main_studies_report, ]$ReportMainN <- dat_reduced[main_studies_report, ]$ReportN
dat_reduced[main_studies_report, ]$ReportMainD <- dat_reduced[main_studies_report, ]$ReportD
dat_reduced$ReportSensN <- dat_reduced$ReportN
dat_reduced$ReportSensD <- dat_reduced$ReportD
dat_reduced$ReportN <- dat_reduced$ReportD <- NULL

################################################################################
## Creating CSV file
################################################################################

row.names(dat_reduced) <- 1:nrow(dat_reduced)
write.csv(dat_reduced, file = '../results/SRInput1.csv')

################################################################################
## Sanity check on subgroup variables
################################################################################

table(dat_reduced$testingsite)
table(dat_reduced$population)
table(dat_reduced$country)
table(dat_reduced$AgName)
table(dat_reduced$Testregime)
table(dat_reduced$Testfrequency)
table(dat_reduced$Sympstat)
table(dat_reduced$Testregime_new)
