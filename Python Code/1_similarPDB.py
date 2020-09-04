#!/usr/bin/env python3

"""
-----------------------------------------------------------------------------------------------
File:	1_similarPDB.py
Program:	similarPDB
Version:	3.0
Created:	28 May 2020
Function:	This code takes the CDR Sequences and searches for similarity within the
          IGHV1-69 germline sequence
Author: 	Ruchira Sachdeva 
------------------------------------------------------------------------------------------------
"""
#***********************************************************************************************
#Import Libraries
import os, sys
import pandas as pd
import csv
import numpy as np
import matplotlib.pyplot as plt
import re

#***********************************************************************************************#

path = os.path.dirname(sys.argv[0])

#Function below allows the user to input a filename and search the relevant CDR1/CDR2/CDR3 sequences
#within the IGHV1-69 germline sequence and create a plot to showcase similarities and differences'

def func(filename):
    """Function below allows the user to input a filename and search the relevant
    CDR1/CDR2/CDR3 sequences within the IGHV1-69 germline sequence"""
    seq = 'MDWTWRFLFVVAAATGVQSQVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCAR'
    pattern = pd.read_csv(filename)
    df = pd.read_csv(filename)
    CDR_Seq = (df[['PDB','CDR1 Sequence','CDR2 Sequence','CDR3 Sequence']])
    pattern_CDR1 = (df[['PDB','CDR1 Sequence','CDR2 Sequence','CDR3 Sequence']])
    PDB_CDR1 = pattern_CDR1['PDB']
    #print(PDB_CDR1)
    wordcount = 0
    for values in CDR_Seq['CDR1 Sequence']:
      if re.findall(values, seq):
          wordcount += 1
          CDR1 = values
          #print(CDR1)
      else:
          continue
    print(wordcount)
    #print(CDR1)
          
func("IGHV1-69_01 Data.csv")
func("IGHV1-69_02 Data.csv")
func("IGHV1-8_01 Data.csv")
func("IGHV1-46_01 Data.csv")
func("IGHV3-21_01 Data.csv")

def dist_graph(filename):
    """Function below plots a histogram to showcase differences between CDR1/CDR2/CDR3
    sequences within the IGHV1-69 germline sequence"""
    df = pd.read_csv(filename)
    CDR = ['CDR1 Len'], ['CDR2 Len', 'CDR3 Len']
    CDR1 = df['CDR1 Len'].values
    CDR2 = df['CDR2 Len'].values
    CDR3 = df['CDR3 Len'].values
    plt.hist((CDR1, CDR2, CDR3), bins = 5, edgecolor='black')
    plt.xlim(8,25)
    plt.legend(['CDR1', 'CDR2', 'CDR3'])
    plt.title('CDR Length Distribution within ' +filename[0:8]+'*'+filename[9:11])
    plt.xlabel('CDR Sequence Length')
    plt.ylabel('CDR Distribution Count')
    plt.tight_layout()
    plt.show()

dist_graph("IGHV1-69_01 Data.csv")
dist_graph("IGHV1-69_02 Data.csv")
dist_graph("IGHV1-8_01 Data.csv")
dist_graph("IGHV1-46_01 Data.csv")
dist_graph("IGHV3-21_01 Data.csv")


