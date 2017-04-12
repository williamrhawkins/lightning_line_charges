# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 21:06:33 2015

Trims un-needed rows and columns from LDAR/LMA data to a format that works
with the MATLAB LCD_MAIN. Save the output from this program to a pos/negX.txt
file for analysis.

Uses Numpy functions.

"""
import numpy as np

fn = input("Enter the filename to trim to fit the MATLAB program. ")
save = input("Enter the name of the .txt file to save as (pos1.txt, neg4.txt, etc.) ")

matrix = np.loadtxt(fn, skiprows=1, usecols=(2,3,4))

np.savetxt(save, matrix)



