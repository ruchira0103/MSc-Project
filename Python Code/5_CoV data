#!/usr/bin/env python3
"""
-----------------------------------------------------------------------------------------------
File:		  5_CoV data.py
Program:	CoV data
Version:	1.0
Created:	20 August 2020
Function:	Code below allows displays the distribution of IGHV germlines in SARS-CoV-2 
Author: 	Ruchira Sachdeva 
------------------------------------------------------------------------------------------------
"""
#***********************************************************************************************
#Import Libraries

import matplotlib.pyplot as plt
import pandas as pd

#***********************************************************************************************#


def pie_chart(infile):
    with open(infile, 'r') as f:
        """Function below allows to display the distribution of IGHV germlines in SARS-CoV-2"""
        df = pd.read_csv(f)
        table = pd.pivot_table(df, values = 'Ab or Nb', index = 'Heavy V Gene', aggfunc = 'count')
        df2 = table.reset_index()
        values = df2['Ab or Nb']
        labels = df2['Heavy V Gene']
        #fig, ax = plt.subplots(nrows=1, ncols=1)
        def autopct_more_than_1(pct):
            return ('%1.f%%' % pct) if pct > 3 else ''
        plt.pie(values,autopct=autopct_more_than_1,
                wedgeprops={'edgecolor': 'white'})
        plt.tight_layout()
        #plt.title('Title')
        plt.legend(labels, loc='best')
        plt.axis('equal')
        plt.show()

pie_chart('CoV-19 Data.csv')

