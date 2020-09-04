#!/usr/bin/env python3

"""
-----------------------------------------------------------------------------------------------
File:		  5_parsecsv.py
Program:	parsecsv
Version:	3.0
Created:	22 August 2020
Function:	This code takes the CSV file from earlier function to capture the CDR ranges
Author: 	Ruchira Sachdeva 
------------------------------------------------------------------------------------------------
"""
#***********************************************************************************************
#Import Libraries
import os, sys
import pandas as pd
import csv
import re

#***********************************************************************************************#

path = os.path.dirname(sys.argv[0])

CDR1 = range(0,30)
CDR2 = range(30-50)
CDR3 = range(50-100)

def parsecsv(infile):
    """Function below uses the csv files produced from earlier function to create a new
    column to capture the residue symbols to create comparisons with logo plots"""
    df = pd.read_csv(infile)
    df['Res Symbol'] = df['Res name'].astype(str).str[0]
    def  expert_level_check(num):
        if 0<= num < 30:
            return 'CDR1'
        elif 31<= num < 50:
            return 'CDR2'
        else:
            return 'CDR3'
    df['CDR'] = df['Res no'].apply(expert_level_check)
    new_rows = (df[['Res name','Res no','Res Symbol', 'CDR']])
    print(new_rows)
    with open('PDBoutput.csv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        for row in new_rows.items():
            csvwriter.writerow(row)


"""IGHV1-69*01"""
parsecsv('PDB-2cmr.csv')
parsecsv('PDB-2fx7.csv')
parsecsv('PDB-2fx8.csv')

"""IGHV1-69*02"""
parsecsv('PDB-1yym.csv')
parsecsv('PDB-2dd8.csv')

"""IGHV1-8*01"""
parsecsv('PDB-3x3f.csv')
parsecsv('PDB-4xmp.csv')

"""IGHV1-46*01"""
parsecsv('PDB-4lsp.csv')
parsecsv('PDB-5f9o.csv')

"""IGHV3-21*01"""
parsecsv('PDB-6cxg.csv')
parsecsv('PDB-6mu3.csv')
