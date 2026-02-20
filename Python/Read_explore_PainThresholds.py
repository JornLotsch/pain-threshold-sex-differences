#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 19:22:32 2022

@author: joern
"""

# %% imports
import os
os.chdir("/home/joern/Aktuell/PainGenesDrugs/08AnalyseProgramme/Python/")
from explore_tukey_lop import explore_tukey_lop
from box_and_heatplots import box_and_heatplot
import copy
import numpy as np
import pandas as pd


# %% Read data
pfad_o = "/home/joern/Aktuell/PainGenesDrugs/"
pfad_u1 = "09Originale/"
filename = pfad_o + pfad_u1 + "Phaeno_125.csv"

dfPainThresholds = pd.read_csv(filename)
dfPainThresholds .columns
dfPainThresholds .set_index("Probandennummer", inplace=True)

# %% Assembling lists of variables

PainThresholdVarnames_nonsensitized = ["vonFrey", "Hitze", "Kaelte", "Druck", "Strom"]
PainThresholdVarnames_sensitized = ["vonFrey_C", "Hitze_C", "Kaelte_M"]

# %% Visual exploration of all variables

dfPainThresholdsAnalyzed = pd.concat([dfPainThresholds[PainThresholdVarnames_nonsensitized], 
                                      dfPainThresholds[PainThresholdVarnames_sensitized]],axis = 1)


for variable in dfPainThresholdsAnalyzed.columns:
    data_subset = copy.copy(
        dfPainThresholdsAnalyzed[[variable]])
    box_and_heatplot(data=data_subset, title=variable,
                     scale=True, cmap="viridis")
    
for variable in dfPainThresholdsAnalyzed.columns:
    data_subset = copy.copy(
        dfPainThresholdsAnalyzed[[variable]])
    explore_tukey_lop(data=data_subset, powers=[-1, -0.5, -0.333333, 0, 0.333333, 0.5, 1, 2])

# %% Transformation according to distribution majority vote

dfPainThresholdsAnalyzed_log = dfPainThresholdsAnalyzed.copy()
for variable in dfPainThresholdsAnalyzed_log.columns:
    data_min = np.nanmin(dfPainThresholdsAnalyzed[[variable]])
    s = np.sign(dfPainThresholdsAnalyzed_log[[variable]])
    dfPainThresholdsAnalyzed_log[[variable]] = np.absolute(dfPainThresholdsAnalyzed_log[[variable]]+1) 
    dfPainThresholdsAnalyzed_log[[variable]] = s*np.log(dfPainThresholdsAnalyzed_log[[variable]].astype("float"))


# %% Computing sensitizing effects

def calculate(x1,x2):
    return(x1 - x2)
    
dfSensEffect = pd.DataFrame(calculate(dfPainThresholdsAnalyzed_log["Hitze"], dfPainThresholdsAnalyzed_log["Hitze_C"]) )
dfSensEffect.columns = ["CapsHeat"]
dfSensEffect["CapsvFrey"] = calculate(dfPainThresholdsAnalyzed_log["vonFrey"], dfPainThresholdsAnalyzed_log["vonFrey_C"]) 
dfSensEffect["MenthCold"] = calculate(dfPainThresholdsAnalyzed_log["Kaelte"], dfPainThresholdsAnalyzed_log["Kaelte_M"]) 

for variable in dfSensEffect.columns:
    data_subset = copy.copy(
        dfSensEffect[[variable]])
    box_and_heatplot(data=data_subset, title=variable,
                     scale=True, cmap="viridis")
    
for variable in dfSensEffect.columns:
    data_subset = copy.copy(
        dfSensEffect[[variable]])
    explore_tukey_lop(data=data_subset, powers=[-1, -0.5, -0.333333, 0, 0.333333, 0.5, 1, 2])
    
# %% Write transformed data
dfPainThresholdsAnalyzed_log_eff = pd.concat([dfPainThresholdsAnalyzed_log, dfSensEffect],axis = 1)
dfPainThresholdsAnalyzed_log_eff.to_csv("dfPainThresholdsAnalyzed_log_eff.csv")
