library(metafor)
library(dplyr)
library(openxlsx)
library(DescTools)

set.seed(1234)

# Load the dataset
data <- read.csv('../results/SRInput1.csv')

# Define the main and subgroup variables
main_variables <- c('UptakeMain', 'UptakeSens', 'ReportMain', 'ReportSens', 'Adherent', 'Ppositive', 'Tpositive', 'FPR', 'Missed', 'Invalid')
subgroup_variables <- c('testingsite', 'population', 'country', 'AgName', 'Testregime_new', 'Sympstat')
subgroup_variables_nice <- c('Testing Site', 'Target Population', 'Study Country', 'Rapid Ag Test', 'Testing Regime', 'Symptom Status')

myround <- function(x, digits){
  return(format(round(x, digits = digits), nsmall = digits))
}

# Function to perform GLMM meta-analysis
run_glmm_meta_analysis <- function(df, var) {
    n_col <- paste0(var, 'N')
    d_col <- paste0(var, 'D')

    if (n_col %in% names(df) && d_col %in% names(df)) {
        df_filtered <- df %>%
                       ungroup() %>%
                       filter(!is.na(df[[n_col]]), !is.na(df[[d_col]]), df[[n_col]] >= 0, df[[d_col]] > 0)

        res <- tryCatch({
            rma.glmm(measure="PLO", xi=df_filtered[[n_col]], ni=df_filtered[[d_col]],
                     method="ML", control=list(verbose=FALSE))
        }, error = function(e) {
            cat("GLMM meta-analysis failed for variable:", var, "\nError message:", e$message, "\n")
            return(NULL)
        }, warning = function(w) {
          rma.glmm(measure="PLO", xi=df_filtered[[n_col]], ni=df_filtered[[d_col]],
                   method="ML", control=list(verbose=FALSE), nAGQ = 10)
        })
        return(res)
    } else {
        return(NULL)
    }
}

# Function to apply inverse logit transformation with vector handling
inv_logit <- function(x) {
  exp(x) / (1 + exp(x))
}

# Initialize Excel workbook
wb <- createWorkbook()
addWorksheet(wb, "Overall Effects")
overall_effects <- data.frame(Variable=character(), Effect_Size=numeric(), CI_Lower=numeric(), CI_Upper=numeric(), I2=numeric(), PI_Lower=numeric(), PI_Upper=numeric())

# Perform meta-analysis and report results
for (var in main_variables) {
  print(paste("Running meta-analysis for:", var))
  overall_result <- run_glmm_meta_analysis(data, var)

  if (!is.null(overall_result)) {
    pred_result <- predict(overall_result, transf=transf.ilogit, addx=TRUE)
    overall_effects <- rbind(overall_effects, data.frame(
      Variable = var,
      Effect_Size = pred_result$pred * 100,
      CI_Lower = pred_result$ci.lb * 100,
      CI_Upper = pred_result$ci.ub * 100,
      I2_Nice = myround(overall_result$I2, 1),
      PI_Nice = paste0('(', 
                       myround(pred_result$pi.lb * 100, 1), 
                       ', ', 
                       myround(pred_result$pi.ub * 100, 1),  
                       ')'),
      I2 = overall_result$I2,
      PI_Lower = pred_result$pi.lb * 100,
      PI_Upper = pred_result$pi.ub * 100, 
      Number_Total=sum(overall_result$ni)
    ))
  } else {
    print(paste("No results for variable:", var))
  }
}

# Write overall effects to Excel
writeData(wb, "Overall Effects", overall_effects)

# Perform subgroup analysis and write results to Excel
for (var in main_variables) {
  addWorksheet(wb, paste(var, "Subgroups", sep="_"))
  subgroup_results <- data.frame(Subgroup_Variable=character(), Subgroup=character(), Effect_Size=numeric(), CI_Lower=numeric(), CI_Upper=numeric(), I2=numeric(), PI_Lower=numeric(), PI_Upper=numeric())
  
  for (subgroup in subgroup_variables) {
    subgroup_nice <- subgroup_variables_nice[which(subgroup == subgroup_variables)]
    print(paste("Running subgroup analysis for:", var, "with subgroup variable:", subgroup))
    unique_subgroups <- unique(na.omit(data[[subgroup]]))
    
    if (length(unique_subgroups) == 0) {
      print(paste("No subgroups found for variable:", subgroup))
    }
    if (subgroup == 'Testregime_new'){
      unique_subgroups <- c('Single', 'Unclear', 'Routine (once weekly)', 'Routine (twice weekly)', 'Routine (three or more times weekly)', 'Routine (unclear frequency)')
    }
    
    if ('General population' %in% unique_subgroups){
      unique_subgroups <- 
        c(unique_subgroups[unique_subgroups == "General population"],
          unique_subgroups[!unique_subgroups %in% c("General population")])
    }
    
    i <- 0
    for (subgroup_val in unique_subgroups) {
      print(paste("Processing subgroup:", subgroup_val, "for variable:", var))
      data_subgroup <- data %>% filter(.[[subgroup]] == subgroup_val, !is.na(.[[paste0(var, 'N')]]), !is.na(.[[paste0(var, 'D')]]))
      
      if (nrow(data_subgroup) >= 4) {
        subgroup_result <- run_glmm_meta_analysis(data_subgroup, var)
        if (!is.null(subgroup_result)) {
          pred_subgroup_result <- predict(subgroup_result, transf=transf.ilogit, addx=TRUE)
          i2_value <- subgroup_result$I2
          subgroup_results <- rbind(subgroup_results, data.frame(
            Subgroup_Variable=subgroup, 
            Subgroup_Variable_Nice=subgroup_nice, 
            Subgroup=subgroup_val, 
            Effect_Size=pred_subgroup_result$pred * 100, 
            CI_Lower=pred_subgroup_result$ci.lb * 100, 
            CI_Upper=pred_subgroup_result$ci.ub * 100, 
            I2_Nice = myround(i2_value, 1),
            PI_Nice = paste0('(', 
                             myround(pred_subgroup_result$pi.lb * 100, 1), 
                             ', ', 
                             myround(pred_subgroup_result$pi.ub * 100, 1),  
                             ')'),
            I2=i2_value, 
            PI_Lower=pred_subgroup_result$pi.lb * 100, 
            PI_Upper=pred_subgroup_result$pi.ub * 100, 
            Number_Total=sum(subgroup_result$ni)))
          
          # Creates separate worksheets for each outcome and subgroup combination, printing the study-specific summary data
          i <- i + 1
          worksheet_title <- paste(var, subgroup, paste0('sub', i), 'Studies', sep="_")
          if (nchar(worksheet_title) > 31) {
            worksheet_title <- substr(worksheet_title, 1, 31)
          }
          addWorksheet(wb, worksheet_title)
          
          # Extract study details for this subgroup
          study_details_subgroup <- data_subgroup %>%
            mutate(Estimate = inv_logit(subgroup_result$yi) * 100, 
                   Subgroup_Value = subgroup_val) %>%
            select(Authors, StudyID, Estimate)
          study_details_subgroup$Subgroup_Value <- subgroup_val
          
          # Write study details to the worksheet
          writeData(wb, worksheet_title, study_details_subgroup)
        }
      }
    }
  }
  writeData(wb, paste(var, "Subgroups", sep="_"), subgroup_results)
}

# Get study-specific score CIs
for (var in main_variables){
  n_col <- paste0(var, 'N')
  d_col <- paste0(var, 'D')
  
  if (n_col %in% names(data) && d_col %in% names(data)) {
    df_filtered <- data %>%
      ungroup() %>%
      filter(!is.na(data[[n_col]]), !is.na(data[[d_col]]), data[[n_col]] >= 0, data[[d_col]] > 0)
    
    ci_lb <- ci_ub <- rep(NA, times = nrow(df_filtered))
    for (i in 1:nrow(df_filtered)){
      temp <- BinomCI(x = df_filtered[[n_col]][i], n = df_filtered[[d_col]][i])
      ci_lb[i] <- 100 * temp[2]
      ci_ub[i] <- 100 * temp[3]
    }
    
    study_details <- data.frame(
      Authors=df_filtered$Authors,
      StudyID=df_filtered$StudyID,
      est=100 * df_filtered[[n_col]] / df_filtered[[d_col]],
      Number=df_filtered[[d_col]],
      ll=ci_lb,
      ul=ci_ub, 
      groupvar = rep('Individual Studies', times = nrow(df_filtered)))
    
    # Create a worksheet for this subgroup
    worksheet_title <- paste(var, "Studies", sep="_")
    if (nchar(worksheet_title) > 31) {
      worksheet_title <- substr(worksheet_title, 1, 31)
    }
    addWorksheet(wb, worksheet_title)
    writeData(wb, worksheet_title, study_details)
  }
}

# Helper variables for re-ordering sheets and adding subgroup studies sheet
all_sheets <- names(wb)
initial_sheets <- c('Overall Effects', paste(main_variables, 'Subgroups', sep = '_'), paste(main_variables, 'Studies', sep = '_'))
remaining_sheets <- setdiff(all_sheets, initial_sheets)

# Adding subgroup studies sheet
df_all <- data.frame(Authors = NULL, StudyID = NULL, Estimate = NULL, Variable = NULL, 
                     Subgroup = NULL, Subgroup_Value = NULL)
for (sheet in remaining_sheets){
  df_temp <- readWorkbook(wb, sheet = sheet)
  df_temp$Subgroup <- sub(".*?_(.*?)_.*", "\\1", sheet)
  df_temp$Variable <- result <- sub("_.*", "", sheet)
  df_all <- rbind(df_all, df_temp)
}
df_all <- df_all[, c('Authors', 'StudyID', 'Estimate', 'Variable', 'Subgroup', 'Subgroup_Value')]

for (var in main_variables){
  df_UptakeMain <- df_all[df_all$Variable == var, ]
  df_UptakeMain$Variable <- NULL
  worksheet_title <- paste0(var, '_Studies_Subgroup')
  addWorksheet(wb, worksheet_title)
  writeData(wb, worksheet_title, df_UptakeMain)
}

worksheetOrder(wb) <- match(c(initial_sheets, 
                              paste0(main_variables, '_Studies_Subgroup'), 
                              remaining_sheets), names(wb))

# Save workbook
saveWorkbook(wb, '../results/meta_analysis_results.xlsx', overwrite = TRUE)

# Sensitivity analysis: Adding value of 1/1 for Lau et al
data_sens <- data[!is.na(data$MissedN) & !is.na(data$MissedD), ]
fit_sens <- rma.glmm(measure="PLO", xi=c(data_sens$MissedN, 1), ni=c(data_sens$MissedD, 1),
                     method="ML", control=list(verbose=FALSE))
predict(fit_sens, transf=transf.ilogit)


# Checking results with conventional LMM meta-analysis
run_lmm_meta_analysis <- function(df, var) {
  n_col <- paste0(var, 'N')
  d_col <- paste0(var, 'D')
  df_filtered <- df %>%
    ungroup() %>%
    filter(!is.na(df[[n_col]]), !is.na(df[[d_col]]), df[[n_col]] >= 0, df[[d_col]] > 0)
  res <-  rma.uni(measure = "PLO", xi=df_filtered[[n_col]], ni=df_filtered[[d_col]], 
                  add = 1)
  return(res)
}

overall_effects_lmm <- data.frame(Variable=character(), Effect_Size=numeric(), CI_Lower=numeric(), CI_Upper=numeric(), I2=numeric(), PI_Lower=numeric(), PI_Upper=numeric())
for (var in main_variables) {
  print(paste("Running meta-analysis for:", var))
  overall_result <- run_lmm_meta_analysis(data, var)
  
  if (!is.null(overall_result)) {
    pred_result <- predict(overall_result, transf=transf.ilogit, addx=TRUE)
    overall_effects_lmm <- rbind(overall_effects_lmm, data.frame(
      Variable = var,
      Effect_Size = pred_result$pred * 100,
      CI_Lower = pred_result$ci.lb * 100,
      CI_Upper = pred_result$ci.ub * 100,
      I2_Nice = myround(overall_result$I2, 1),
      PI_Nice = paste0('(', 
                       myround(pred_result$pi.lb * 100, 1), 
                       ', ', 
                       myround(pred_result$pi.ub * 100, 1),  
                       ')'),
      I2 = overall_result$I2,
      PI_Lower = pred_result$pi.lb * 100,
      PI_Upper = pred_result$pi.ub * 100, 
      Number_Total=sum(overall_result$ni)
    ))
  } else {
    print(paste("No results for variable:", var))
  }
}
overall_effects_lmm

df1 <- overall_effects[!overall_effects$Variable %in% c('UptakeSens', 'ReportSens'), ]
df2 <- overall_effects_lmm[!overall_effects_lmm$Variable %in% c('UptakeSens', 'ReportSens'), ]

k <- nrow(df1)
labels <- c('Test Uptake (%)', 'Test Reporting (%)', 
            'Adherent (%)','New Cases (%)', 
            'Tests Positive (%)', 'False Positives (%)', 
            'Missed Cases (%)', 'Invalid Tests (%)')

col1 <- "#D95F02"; col2 <- "#1B9E77" 
pch1 <- 15; pch2 <- 19 

y_pos <- k:1
y_labels <- labels

offset <- 0.15
xlim_data <- range(c(df1$CI_Lower, df1$CI_Upper, 
                     df2$CI_Lower, df2$CI_Upper), na.rm = TRUE)
xlim_pad <- (xlim_data[2] - xlim_data[1]) * 0.05
xlim_final <- c(xlim_data[1] - xlim_pad, xlim_data[2] + xlim_pad)

pdf('../results/figures/sensitivity-v2.pdf', width = 7, height = 5.5)
layout(matrix(c(1, 2), nrow = 2), heights = c(0.85, 0.15))
par(mar = c(4, 9, 3, 2))

plot(
  NA,
  xlim = xlim_final,
  ylim = c(0.5, k + 0.5), # y-range from 0.5 to k+0.5
  yaxt = "n",             # Turn off default y-axis
  xaxt = "n",             # Turn off default x-axis
  ylab = "",
  xlab = "",
  bty = "n"               # No box around the plot
)
abline(h = y_pos, col = "grey90", lty = "dotted")
abline(v = 0, col = "grey60", lty = "dashed", lwd = 1.5)

axis(1, at = pretty(xlim_final), cex.axis = 0.9)
title(xlab = "Pooled Estimate (%)", line = 2.5, cex.lab = 1)

axis(2, at = y_pos, labels = y_labels, las = 1, tick = FALSE, line = -1, cex.axis = 0.9)

segments(df1$CI_Lower, y_pos + offset,
         df1$CI_Upper, y_pos + offset,
         col = col1, lwd = 2)
points(df1$Effect_Size, y_pos + offset,
       pch = pch1, col = col1)

segments(df2$CI_Lower, y_pos - offset,
         df2$CI_Upper, y_pos - offset,
         col = col2, lwd = 2)
points(df2$Effect_Size, y_pos - offset,
       pch = pch2, col = col2)

par(mar = c(1, 9, 0, 2))
plot(0, 0, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
legend(
  "center",           
  legend = c("GLMM", "Inverse-Variance Method"),
  col = c(col1, col2),
  pch = c(pch1, pch2),
  lty = 1,            
  lwd = 2,            
  bty = "n",          
  horiz = TRUE,      
)

dev.off()


wb_sensitivity <- createWorkbook()
addWorksheet(wb_sensitivity, "Inverse Variance Method")
writeData(wb_sensitivity, "Inverse Variance Method", overall_effects_lmm)
saveWorkbook(wb_sensitivity, '../results/meta_analysis_results_sensitivity.xlsx', overwrite = TRUE)


## Funnel plots
fit_glmm <- run_glmm_meta_analysis(data, 'UptakeMain')
pdf('../results/figures/UptakeMain-Funnel.pdf', width = 5.5, height = 5.5)
par(mar = c(5, 5, 3, 2))
funnel(fit_glmm, yaxis = 'sqrtni', 
       xlab = 'Log Odds of Test Uptake Proportion', main = 'Test Uptake',
       digits = c(0, 1), cex = 1.1, cex.main = 1.25)
dev.off()

fit_glmm <- run_glmm_meta_analysis(data, 'ReportMain')
pdf('../results/figures/ReportMain-Funnel.pdf', width = 5.5, height = 5.5)
funnel(fit_glmm, yaxis = 'sqrtni', 
       xlab = 'Log Odds of Test Reporting Proportion', main = 'Testing Reporting',
       digits = c(0, 1), cex = 1.1, cex.main = 1.25)
dev.off()

fit_glmm <- run_glmm_meta_analysis(data, 'Adherent')
pdf('../results/figures/Adherent-Funnel.pdf', width = 5.5, height = 5.5)
funnel(fit_glmm, yaxis = 'sqrtni', 
       xlab = 'Log Odds of Adherent Proportion', main = 'Adherent', 
       digits = c(0, 1), cex = 1.1, cex.main = 1.25)
dev.off()

fit_glmm <- run_glmm_meta_analysis(data, 'Ppositive')
pdf('../results/figures/Ppositive-Funnel.pdf', width = 5.5, height = 5.5)
funnel(fit_glmm, yaxis = 'sqrtni', 
       xlab = 'Log Odds of New Cases Proportion', main = 'New Cases',
       digits = c(0, 1), cex = 1.1, cex.main = 1.25)
dev.off()

fit_glmm <- run_glmm_meta_analysis(data, 'Tpositive')
pdf('../results/figures/Tpositive-Funnel.pdf', width = 5.5, height = 5.5)
funnel(fit_glmm, yaxis = 'sqrtni', 
       xlab = 'Log Odds of Tests Positive Proportion', main = 'Tests Positive', 
       digits = c(0, 1), cex = 1.1, cex.main = 1.25)
dev.off()

fit_glmm <- run_glmm_meta_analysis(data, 'FPR')
pdf('../results/figures/FPR-Funnel.pdf', width = 5.5, height = 5.5)
funnel(fit_glmm, yaxis = 'sqrtni', 
       xlab = 'Log Odds of False Positives Proportion', main = 'False Positives',
       digits = c(0, 1), cex = 1.1, cex.main = 1.25)
dev.off()

fit_glmm <- run_glmm_meta_analysis(data, 'Missed')
pdf('../results/figures/Missed-Funnel.pdf', width = 5.5, height = 5.5)
funnel(fit_glmm, yaxis = 'sqrtni', 
       xlab = 'Log Odds of Missed Persons Proportion', main = 'Missed Persons',
       digits = c(0, 1), cex = 1.1, cex.main = 1.25)
dev.off()

fit_glmm <- run_glmm_meta_analysis(data, 'Invalid')
pdf('../results/figures/Invalid-Funnel.pdf', width = 5.5, height = 5.5)
funnel(fit_glmm, yaxis = 'sqrtni', 
       xlab = 'Log Odds of Invalid Tests Proportion', main = 'Invalid Tests',
       digits = c(0, 1), cex = 1.1, cex.main = 1.25)
dev.off()
