#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 20:27:56 2024

@author: Sean
"""
import forestplot as fp
import pandas as pd 
import matplotlib.pyplot as plt


file_path_subgroup ="../results/meta_analysis_results.xlsx"
dataUMS= pd.read_excel(file_path_subgroup, sheet_name="UptakeMain_Subgroups")
dataUSS= pd.read_excel(file_path_subgroup, sheet_name="UptakeSens_Subgroups")
dataRMS= pd.read_excel(file_path_subgroup, sheet_name="ReportMain_Subgroups")
dataRSS= pd.read_excel(file_path_subgroup, sheet_name="ReportSens_Subgroups")
dataAS= pd.read_excel(file_path_subgroup, sheet_name="Adherent_Subgroups")
dataPS= pd.read_excel(file_path_subgroup, sheet_name="Ppositive_Subgroups")
dataTS= pd.read_excel(file_path_subgroup, sheet_name="Tpositive_Subgroups")
dataFS= pd.read_excel(file_path_subgroup, sheet_name="FPR_Subgroups")
dataMS= pd.read_excel(file_path_subgroup, sheet_name="Missed_Subgroups")
dataIS= pd.read_excel(file_path_subgroup, sheet_name="Invalid_Subgroups")

file_path_overall ="../results/forest_plot_data.xlsx"
dataUM= pd.read_excel(file_path_overall, sheet_name="UptakeMain")
dataUS= pd.read_excel(file_path_overall, sheet_name="UptakeSens")
dataRM= pd.read_excel(file_path_overall, sheet_name="ReportMain")
dataRS= pd.read_excel(file_path_overall, sheet_name="ReportSens")
dataA= pd.read_excel(file_path_overall, sheet_name="Adherent")
dataP= pd.read_excel(file_path_overall, sheet_name="Ppositive")
dataT= pd.read_excel(file_path_overall, sheet_name="Tpositive")
dataF= pd.read_excel(file_path_overall, sheet_name="FPR")
dataM= pd.read_excel(file_path_overall, sheet_name="Missed")
dataI= pd.read_excel(file_path_overall, sheet_name="Invalid")





def save_plot_overall(data, estimate, ll, hl, varlabel, annote, annoteheaders, xlabel, groupvar, plot_name, height):
    fig, ax = plt.subplots(figsize=(6, height))
    fp.forestplot(data, fig = fig, ax = ax, decimal_precision = 1,
                  estimate=estimate, ll=ll, hl=hl, varlabel=varlabel, 
                       annote=annote, annoteheaders=annoteheaders, xlabel=xlabel, 
                       groupvar=groupvar, table=True, **{"marker": "D", "markersize": 35, 
                       "xline": [0, 100], "xlinestyle": (0, (10, 5)), "xlinecolor": ".1", "xtick_size": 12, "table_label_size": 16})
    
    xmin, xmax = ax.get_xlim()
    if round(xmax) == 100:
        ax.set_xlim(left=xmin, right=100.25)  # Extend upper limit slightly so that the line shows

    
    plt.savefig(f"{plot_name}.png", bbox_inches='tight', transparent=True)
    plt.close()

def save_plot(data, estimate, ll, hl, varlabel, annote, annoteheaders, xlabel, groupvar, rightannote, right_annoteheaders, plot_name, height):
    fig, ax = plt.subplots(figsize=(6, height))
    fp.forestplot(data, fig = fig, ax = ax, decimal_precision = 1,
                  estimate=estimate, ll=ll, hl=hl, varlabel=varlabel, 
                       annote=annote, annoteheaders=annoteheaders, xlabel=xlabel, 
                       groupvar=groupvar, rightannote=rightannote, right_annoteheaders=right_annoteheaders, table=True, **{"marker": "D", "markersize": 35, 
                       "xline": [0, 100], "xlinestyle": (0, (10, 5)), "xlinecolor": ".1", "xtick_size": 12, "table_label_size": 16})
    
    xmin, xmax = ax.get_xlim()
    if round(xmax) == 100:
        ax.set_xlim(left=xmin, right=100.25)  # Extend upper limit slightly so that the line shows
    
    plt.savefig(f"{plot_name}.png", bbox_inches='tight', transparent=True)
    plt.close()



save_plot(dataUMS, "Effect_Size", "CI_Lower", "CI_Upper", "Subgroup", ["est_ci"], 
          ["Point Estimate (95% CI)"], "Test Uptake (%)", "Subgroup_Variable_Nice", ["PI_Nice", "I2_Nice"], ["95% PI", "I2"], "../results/figures/UptakeMainSubgroupPlot", 8)

save_plot(dataUSS, "Effect_Size", "CI_Lower", "CI_Upper", "Subgroup", ["est_ci"], 
          ["Point Estimate (95% CI)"], "Test Uptake (%)", "Subgroup_Variable_Nice", ["PI_Nice", "I2_Nice"], ["95% PI", "I2"], "../results/figures/UptakeSensSubgroupPlot", 8)

save_plot(dataRMS, "Effect_Size", "CI_Lower", "CI_Upper", "Subgroup", ["est_ci"], 
          ["Point Estimate (95% CI)"], "Test Reporting (%)", "Subgroup_Variable_Nice", ["PI_Nice", "I2_Nice"], ["95% PI", "I2"], "../results/figures/ReportMainSubgroupPlot", 9)

save_plot(dataRSS, "Effect_Size", "CI_Lower", "CI_Upper", "Subgroup", ["est_ci"], 
          ["Point Estimate (95% CI)"], "Test Reporting (%)", "Subgroup_Variable_Nice", ["PI_Nice", "I2_Nice"], ["95% PI", "I2"], "../results/figures/ReportSensSubgroupPlot", 9)

save_plot(dataAS, "Effect_Size", "CI_Lower", "CI_Upper", "Subgroup", ["est_ci"], 
          ["Point Estimate (95% CI)"], "Adherent (%)", "Subgroup_Variable_Nice", ["PI_Nice", "I2_Nice"], ["95% PI", "I2"], "../results/figures/AdherentSubgroupPlot", 10)

save_plot(dataPS, "Effect_Size", "CI_Lower", "CI_Upper", "Subgroup", ["est_ci"], 
          ["Point Estimate (95% CI)"], "New Cases (%)", "Subgroup_Variable_Nice", ["PI_Nice", "I2_Nice"], ["95% PI", "I2"], "../results/figures/PpositiveSubgroupPlot", 10)

save_plot(dataTS, "Effect_Size", "CI_Lower", "CI_Upper", "Subgroup", ["est_ci"], 
          ["Point Estimate (95% CI)"], "Tests Positive (%)", "Subgroup_Variable_Nice", ["PI_Nice", "I2_Nice"], ["95% PI", "I2"], "../results/figures/TpositiveSubgroupPlot", 10)

save_plot(dataFS, "Effect_Size", "CI_Lower", "CI_Upper", "Subgroup", ["est_ci"], 
          ["Point Estimate (95% CI)"], "False Positives (%)", "Subgroup_Variable_Nice", ["PI_Nice", "I2_Nice"], ["95% PI", "I2"], "../results/figures/FPRSubgroupPlot", 6)

# See plot below
# save_plot(dataMS, "Effect_Size", "CI_Lower", "CI_Upper", "Subgroup", ["est_ci"], 
#           ["Point Estimate (95% CI)"], "Missed Persons (%)", "Subgroup_Variable_Nice", ["PI_Nice", "I2_Nice"], ["95% PI", "I2"], "../results/figures/MissedSubgroupPlot", 4)

save_plot(dataIS, "Effect_Size", "CI_Lower", "CI_Upper", "Subgroup", ["est_ci"], 
          ["Point Estimate (95% CI)"], "Invalid Tests (%)", "Subgroup_Variable_Nice", ["PI_Nice", "I2_Nice"], ["95% PI", "I2"], "../results/figures/InvalidSubgroupPlot", 8)





save_plot_overall(dataUM, "est", "ll", "ul", "Authors", ["Number", "est_ci"], 
          ["N", "Point Estimate (95% CI)"], "Test Uptake (%) ", "groupvar", "../results/figures/UptakeMainPlot", 10)

save_plot_overall(dataUS, "est", "ll", "ul", "Authors", ["Number", "est_ci"], 
          ["N", "Point Estimate (95% CI)"], "Test Uptake (%) ", "groupvar", "../results/figures/UptakeSensPlot", 11)

save_plot_overall(dataRM, "est", "ll", "ul", "Authors", ["Number", "est_ci"], 
          ["N", "Point Estimate (95% CI)"], "Test Reporting (%)", "groupvar", "../results/figures/ReportMainPlot", 12)

save_plot_overall(dataRS, "est", "ll", "ul", "Authors", ["Number", "est_ci"], 
          ["N", "Point Estimate (95% CI)"], "Test Reporting (%)", "groupvar", "../results/figures/ReportSensPlot", 12)

save_plot_overall(dataA, "est", "ll", "ul", "Authors", ["Number", "est_ci"], 
          ["N", "Point Estimate (95% CI)"], "Adherent (%) ", "groupvar", "../results/figures/AdherentPlot", 13)

save_plot_overall(dataP, "est", "ll", "ul", "Authors", ["Number", "est_ci"], 
          ["N", "Point Estimate (95% CI)"], "New Cases (%)", "groupvar", "../results/figures/PpositivePlot", 15)

save_plot_overall(dataT, "est", "ll", "ul", "Authors", ["Number", "est_ci"], 
          ["N", "Point Estimate (95% CI)"], "Tests Positive (%) ", "groupvar", "../results/figures/TpositivePlot", 15)

save_plot_overall(dataF, "est", "ll", "ul", "Authors", ["Number", "est_ci"], 
          ["N", "Point Estimate (95% CI)"], "False Positives (%)", "groupvar", "../results/figures/FPRPlot", 5.5)

save_plot_overall(dataM, "est", "ll", "ul", "Authors", ["Number", "est_ci"], 
          ["N", "Point Estimate (95% CI)"], "Missed Persons (%)", "groupvar", "../results/figures/MissedPlot", 5)

save_plot_overall(dataI, "est", "ll", "ul", "Authors", ["Number", "est_ci"], 
          ["N", "Point Estimate (95% CI)"], "Invalid Tests (%)", "groupvar", "../results/figures/InvalidPlot", 10)





# Separate code for plotting the Missed subgroup analyses due to the bug in generating a header 
fig, ax = plt.subplots(figsize=(6, 3))

fp.forestplot(dataMS, fig=fig, ax=ax, decimal_precision=1,
              estimate="Effect_Size", ll="CI_Lower", hl="CI_Upper", 
              varlabel="Subgroup", annote=["est_ci"], 
              annoteheaders=["Point Estimate (95% CI)"], 
              xlabel="Missed Persons (%)", groupvar="Subgroup_Variable_Nice", 
              rightannote=["PI_Nice", "I2_Nice"], 
              right_annoteheaders=["95% PI", "I2"], table=True, 
              **{"marker": "D", "markersize": 35, "xline": [0, 100], 
                 "xlinestyle": (0, (10, 5)), "xlinecolor": ".1", 
                 "xtick_size": 12, "table_label_size": 16})

# Adjusting header position manually
ax.text(-17.5, 3.8, "Point Estimate (95% CI)", 
        fontsize=12, ha="center", va="bottom", 
        fontfamily="monospace",
        fontweight="semibold", 
        color="black")

ax.text(-60.75, 3.8, "Variable", 
        fontsize=12, ha="center", va="bottom", 
        fontfamily="monospace",
        fontweight="semibold",  
        color="black")

xmin, xmax = ax.get_xlim()
if round(xmax) == 100:
    ax.set_xlim(left=xmin, right=100.25)  # Extend upper limit slightly so that the line shows

plt.savefig("../results/figures/MissedSubgroupPlot.png", bbox_inches='tight', transparent=True)
plt.close()



