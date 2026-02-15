library(openxlsx)

# Load the meta-analysis results workbook
wb_input <- loadWorkbook('../results/meta_analysis_results.xlsx')

# Read the Overall Effects sheet
overall_effects <- readWorkbook(wb_input, sheet = 'Overall Effects')

# Get variables
all_sheets <- names(wb_input)
studies_sheets <- all_sheets[grepl('_Studies$', all_sheets) & 
                             !grepl('_Subgroup$', all_sheets) &
                             !grepl('_sub[0-9]+_Studies', all_sheets)]
variables <- sub('_Studies$', '', studies_sheets)

# Process each variable
wb_output <- createWorkbook()
for (var in variables) {
  # Get the corresponding row from Overall Effects
  overall_row <- overall_effects[overall_effects$Variable == var, ]
  
  # Skip if no overall effect found
  if (nrow(overall_row) == 0) {
    warning(paste("No overall effect found for variable:", var))
    next
  }
  
  # Read the individual studies data
  studies_sheet_name <- paste0(var, '_Studies')
  if (!studies_sheet_name %in% all_sheets) {
    warning(paste("Studies sheet not found for variable:", var))
    next
  }
  
  studies_data <- readWorkbook(wb_input, sheet = studies_sheet_name)
  if ("StudyID" %in% names(studies_data)) {
    studies_data <- studies_data[, !names(studies_data) %in% "StudyID"]
  }
  
  # Create the summary row
  summary_authors <- paste0("Summary (PI=", overall_row$PI_Nice, ", I2=", overall_row$I2_Nice, "%)")
  summary_row <- data.frame(
    Authors = summary_authors,
    est = overall_row$Effect_Size,
    Number = "_",
    ll = overall_row$CI_Lower,
    ul = overall_row$CI_Upper,
    groupvar = "Overall"
  )
  
  # Combine individual studies with summary row
  combined_data <- rbind(studies_data, summary_row)
  
  # Add worksheet to output workbook
  addWorksheet(wb_output, var)
  writeData(wb_output, var, combined_data)
}

# Save the workbook
saveWorkbook(wb_output, '../results/forest_plot_data.xlsx', overwrite = TRUE)
