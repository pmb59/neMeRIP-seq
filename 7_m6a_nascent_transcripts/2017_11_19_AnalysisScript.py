##################################################################
# Script Author: Chris Sadee
# start Date: 	16/11/17
# end Date: 	19/11/17
# Analysis of FC peaks within genes
##################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re


##################################################################
# Import data
##################################################################

namData='Subset3A_ann.bed'

#Reading data into pandas dataframe
Data=pd.read_table(namData,header=0,index_col=None)

#sort data by ENSEMBL_ID for simplicity
Data.sort_values(["ENSEMBL_ID"])


##################################################################
# Append columns
##################################################################

nameCol = ["Regulation"]

#Append
Data[nameCol[0]]= "None"

##################################################################
# calculate Up or Down reagulation
##################################################################

#Cut offs:
fhDown = np.log2(2.0/3) # <------ Select fresh hold for Downregulation
fhUp   = 1.5 # <------ Select fresh hold for Upregulation

#Significant values:
fhFDR = np.log10(0.05) # <------ Select fresh hold for FDR


#Conditions:
condFDR = Data["diff.lg.fdr"] < fhFDR

condDown= (Data["diff.log2.fc"] < fhDown) & condFDR
condUp  = (Data["diff.log2.fc"] > fhUp  ) & condFDR

Data.loc[condDown,nameCol[0]]="Down"
Data.loc[condUp,nameCol[0]]="up"

##################################################################
# Importing data frame into dictionary for ENSEMBL_ID key
##################################################################

IDs   = Data["ENSEMBL_ID"].unique()
nIDs, = IDs.shape

Dict=dict.fromkeys(IDs)

for i in range(0,nIDs):

	condIDs = Data["ENSEMBL_ID"] == IDs[i]

	Dict[IDs[i]] = Data[condIDs]


##################################################################
# Features per ENSEMBL ID compare min exon peak with mean of exon peaks
##################################################################

#Initializing OutPut Dataframe:
outDF1 = pd.DataFrame(columns=Dict[IDs[0]].columns.values)
#Adding Mean column
outDF1["Mean"]=[]

#Initializing iterator
iter1  = 0

#Looping through keys in Dict
for i in range(0,nIDs):

	#Re-casting into Data Frame for convinience (note this is slow)
	ensData=Dict[IDs[i]]

	##############################################################
	# Conditions to Pass
	
	#Condition that there are more than 1 exonic Peak
	condEx = ensData["Functional_Location"] == "Exon"
	condMultEx = condEx.sum() > 1
	if condMultEx == False: continue

	#Condition that there is at least one exonic peak Downregulated with fc>1.5 and fdr>0.05
	condDownReg = ((ensData[nameCol[0]] == "Down") & condEx).sum() > 0
	if condDownReg == False: continue


	##############################################################
	# Calaculations

	#Location of most downregulated peak in an exon and satisfying FDR
	idxMin = ensData["diff.log2.fc"][(ensData[nameCol[0]] == "Down") & condEx].idxmin()

	#Filling values into outPut Dataframe1
	outDF1.loc[iter1,:] = ensData.loc[idxMin,:]

	#Condition to find Values to use for mean, excluding min peak and using only exons
	condMean = (ensData.index.values != idxMin) & condEx


	outDF1.loc[iter1,"Mean"] = ensData["diff.log2.fc"][condMean].mean()


	iter1+=1

#Saving outDF1
outDF1.to_csv("Subset3A_ann_ChrisPoint1.bed",sep='\t',index=False,header=True, encoding='utf-8')


##################################################################
# Features per ENSEMBL ID compare min exon peak with mean of intron peaks
##################################################################


#Initializing OutPut Dataframe:
outDF2 = pd.DataFrame(columns=Dict[IDs[0]].columns.values)
#Adding Mean column
outDF2["Mean"]=[]

#Initializing iterator
iter2  = 0


#Looping through keys in Dict
for i in range(0,nIDs):

	#Re-casting into Data Frame for convinience (note this is slow)
	ensData=Dict[IDs[i]]

	##############################################################
	# Conditions to Pass
	
	#Condition that there is at least 1 exonic Peak
	condEx = ensData["Functional_Location"] == "Exon"
	condMultEx = condEx.sum() > 0
	if condMultEx == False: continue

	#Condition that there is at least 1 intronic Peak
	condIn = ensData["Functional_Location"] == "Intron"
	condMultIn = condIn.sum() > 0
	if condMultIn == False: continue


	#Condition that there is at least one exonic peak Downregulated with fc<fhDown and fdr<fhFDR
	condDownReg = ((ensData[nameCol[0]] == "Down") & condEx).sum() > 0
	if condDownReg == False: continue


	##############################################################
	# Calaculations

	#Location of most downregulated peak in an exon and satisfying FDR
	idxMin = ensData["diff.log2.fc"][(ensData[nameCol[0]] == "Down") & condEx].idxmin()

	#Filling values into outPut Dataframe2
	outDF2.loc[iter2,:] = ensData.loc[idxMin,:]

	#Condition to find Values to use for mean, using only introns
	condMean = condIn


	outDF2.loc[iter2,"Mean"] = ensData["diff.log2.fc"][condMean].mean()

	iter2+=1


outDF2.to_csv("Subset3A_ann_ChrisPoint2.bed",sep='\t',index=False,header=True, encoding='utf-8')


##################################################################
# Features per ENSEMBL ID compare min exon peak with mean everything else
##################################################################


#Initializing OutPut Dataframe:
outDF3 = pd.DataFrame(columns=Dict[IDs[0]].columns.values)
#Adding Mean column
outDF3["Mean"]=[]

#Initializing iterator
iter3  = 0


#Looping through keys in Dict
for i in range(0,nIDs):

	#Re-casting into Data Frame for convinience (note this is slow)
	ensData=Dict[IDs[i]]

	##############################################################
	# Conditions to Pass
	
	#Condition that there is at least 1 exonic Peak
	condEx = ensData["Functional_Location"] == "Exon"
	condMultEx = condEx.sum() > 0
	if condMultEx == False: continue

	#Condition that there are more than one peak of any kind
	condMultPeak = ensData["Functional_Location"].shape[0] > 1
	if condMultPeak == False: continue


	#Condition that there is at least one exonic peak Downregulated with fc<fhDown and fdr<fhFDR
	condDownReg = ((ensData[nameCol[0]] == "Down") & condEx).sum() > 0
	if condDownReg == False: continue


	##############################################################
	# Calaculations

	#Location of most downregulated peak in an exon and satisfying FDR
	idxMin = ensData["diff.log2.fc"][(ensData[nameCol[0]] == "Down") & condEx].idxmin()

	#Filling values into outPut Dataframe3
	outDF3.loc[iter3,:] = ensData.loc[idxMin,:]

	#Condition to find Values to use for mean, excluding min peak
	condMean = (ensData.index.values != idxMin)


	outDF3.loc[iter3,"Mean"] = ensData["diff.log2.fc"][condMean].mean()


	iter3+=1



outDF3.to_csv("Subset3A_ann_ChrisPoint3.bed",sep='\t',index=False,header=True, encoding='utf-8')


#################################################################
#Features per ENSEMBL ID find position of min exon peak and distance to any other peak
#################################################################


#Initializing OutPut Dataframe:
outDF4 = pd.DataFrame(columns=Dict[IDs[0]].columns.values)
#Adding Min Peak column
outDF4["Min_Peak"]=[]
#Centre_Peak_Position
outDF4["Centre_Peak_Position"]=[]
#Adding Distance column
outDF4["Distance"]=[]


#Looping through keys in Dict
for i in range(0,nIDs):

	#Re-casting into Data Frame for convinience (note this is slow)
	ensData=Dict[IDs[i]]

	##############################################################
	# Conditions to Pass
	
	#Condition that there is at least 1 exonic Peak
	condEx = ensData["Functional_Location"] == "Exon"
	condMultEx = condEx.sum() > 0
	if condMultEx == False: continue

	#Condition that there are more than one peak of any kind
	condMultPeak = ensData["Functional_Location"].shape[0] > 1
	if condMultPeak == False: continue


	#Condition that there is at least one exonic peak Downregulated with fc<fhDown and fdr<fhFDR
	condDownReg = ((ensData[nameCol[0]] == "Down") & condEx).sum() > 0
	if condDownReg == False: continue

	##############################################################
	# Min Peak Calculations

	#Making Temporary Dataframe with same columns as outDF4
	tempDF=pd.DataFrame(columns=outDF4.columns.values)

	#Filling values into temp Dataframe
	tempDF=tempDF.append(ensData, ignore_index=True)[outDF4.columns.tolist()]

	#Condition of most downregulated exon and satisfying FDR
	condEx = (tempDF[nameCol[0]] == "Down") & (tempDF["Functional_Location"] == "Exon")

	#Location of most downregulated peak in an exon 
	idxMin = tempDF["diff.log2.fc"][condEx].idxmin()

	#Condition to find other peaks, excluding min peak
	condPeaks = (tempDF.index.values != idxMin)

	#Indicating min peak
	tempDF.loc[idxMin,"Min_Peak"]="Y"

	#Indication all other peaks
	tempDF.loc[:,"Min_Peak"][condPeaks]="N"

	#Determining centre position of all peak
	meanPos=(tempDF.loc[:,"chromEnd"]-tempDF.loc[:,"chromStart"])/2
	tempDF.loc[:,"Centre_Peak_Position"] = np.floor(meanPos)+ tempDF.loc[:,"chromStart"]	

	#Determining distance from min peak to all other peaks
	tempDF["Distance"]= tempDF.loc[:,"Centre_Peak_Position"]-tempDF.loc[idxMin,"Centre_Peak_Position"]

	outDF4=outDF4.append(tempDF, ignore_index=True)[outDF4.columns.tolist()]


outDF4.to_csv("Subset3A_ann_ChrisPoint4.bed",sep='\t',index=False,header=True, encoding='utf-8')


