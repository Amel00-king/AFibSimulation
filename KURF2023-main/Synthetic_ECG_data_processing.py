# -*- coding: utf-8 -*-
"""
Created on Mon May 29 12:24:08 2023

@author: olegadmin
"""
import numpy as np
import pylab as plt
import os

#---------------------
# Get all file names in the folder with data:
# folder path
dir_path = 'C:/Users/olegadmin/Documents/KURF2023/Data_27May_1'
# list to store files
fileNames = []

# Iterate directory
for path in os.listdir(dir_path):
    # check if current path is a file
    if os.path.isfile(os.path.join(dir_path, path)):
        fileNames.append(path)
print(fileNames)

#--------------------------

# Loop through all filenames to check their ECGs are ok:
Spiral_or_FocalPoint = []
for f in range(len(fileNames)):
    rel_path = fileNames[f]
    # Get coordinates from title of .txt file:
    c1Start = rel_path.find('on_') + 3
    c1End = rel_path.find('com')
    c2Start = rel_path.find('mma') + 3
    c2End = rel_path.find('_Lii')
    C1 = int(rel_path[c1Start:c1End])
    C2 = int(rel_path[c2Start:c2End])
    Coords = [C1,C2]
    
    # Check if spiral or focal point:
    if (fileNames[f].count('FP') >  0):
        Spiral_or_FocalPoint.append('FP')
        title = 'Focal point at [' + str(Coords) + ']'
    elif (fileNames[f].count('SP') >  0):
        Spiral_or_FocalPoint.append('SP')
        title = 'Spiral initiated at [' + str(Coords) + ']'
    abs_file_path = os.path.join(dir_path, rel_path)
    # Opening file:
    ECGfile = open(abs_file_path, 'r');
    count = 0
    
    # Get total size of the text file:
    for lineCount in ECGfile:
        if count<1:
            listFloats = lineCount.split("\t")
            x = list(np.float_(listFloats))
            num_traces = len(x)
        count += 1
    
    # Pre-allocate size of wave array for faster assignment **MAY NOT BE NECESSARY**
    Pwaves_Lii = np.zeros(count)
    Pwaves_Liii = np.zeros(count)
    Pwaves_WCT = np.zeros(count)
    Pwaves_V1 = np.zeros(count)
    Pwaves_V2 = np.zeros(count)
    Pwaves_V3 = np.zeros(count)
    Pwaves_V4 = np.zeros(count)
    Pwaves_V5 = np.zeros(count)
    Pwaves_V6 = np.zeros(count)
    Pwaves_aVF = np.zeros(count)
     
    ECGfile = open(abs_file_path, 'r');
    # Using for loop to read in the data:
    print("Using for loop")
    iteration = 0
    for lineRead in ECGfile:
        listFloats = lineRead.split("\t")
        x = list(np.float_(listFloats))
        Pwaves_Lii[iteration] = x[0]
        Pwaves_Liii[iteration] = x[1]
        Pwaves_WCT[iteration] = x[2]
        Pwaves_V1[iteration] = x[3]
        Pwaves_V2[iteration] = x[4]
        Pwaves_V3[iteration] = x[5]
        Pwaves_V4[iteration] = x[6]
        Pwaves_V5[iteration] = x[7]
        Pwaves_V6[iteration] = x[8]
        Pwaves_aVF[iteration] = x[9]
        iteration += 1

    # Closing files
    ECGfile.close()
    
    plt.plot(np.arange(count),Pwaves_Lii)
    plt.plot(np.arange(count),Pwaves_Liii)
    plt.plot(np.arange(count),Pwaves_WCT)
    plt.plot(np.arange(count),Pwaves_V1)
    plt.plot(np.arange(count),Pwaves_V2)
    plt.plot(np.arange(count),Pwaves_V3)
    plt.plot(np.arange(count),Pwaves_V4)
    plt.plot(np.arange(count),Pwaves_V5)
    plt.plot(np.arange(count),Pwaves_V6)
    plt.plot(np.arange(count),Pwaves_aVF)
    plt.title(title)
    plt.show()