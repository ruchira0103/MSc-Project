#!/usr/bin/env python3

"""
-----------------------------------------------------------------------------------------------
File:	4_filtercsv.py
Program:	filtercsv
Version:	3.0
Created:	28 July 2020
Function:	This code takes the CDR Sequences and filters, downloads and parses files within
          the various IGHV germlines
Author: 	Ruchira Sachdeva 
------------------------------------------------------------------------------------------------
"""
#***********************************************************************************************
#Import Libraries
import pandas as pd
import csv
import re
import zipfile
import requests

#***********************************************************************************************#
        
def filter_csv(infile):
    """Function below allows filtering of IMGT data by capturing PDBs with Resolution < 2.5 A,
    capitalizing the IMGT ID column and saving to new CSV file"""
    with open(infile, 'r') as csvfile:
        df = pd.read_csv(csvfile, sep='\t', engine='python', delimiter = ',')
        df = df[df['Resolution'] <= 2.5]
        df['IMGT entry ID'] = df['IMGT entry ID'].str.upper() 
        compression_opts = dict(method='zip',archive_name='filtered_IMGT.csv') #save to zip file
        df.to_csv('out.zip', index=False, compression=compression_opts) 
        zip = zipfile.ZipFile('out.zip') #extract file here
        zip.extractall()

filter_csv('Ig-Ag Complex.csv')

def matchingPDB(infile1, infile2, outfile):
    """Function allows to search PDB codes from one file to another and save in a new CSV file"""
    PDBs_sliced = [] #create an empty array 
    #entering PDB codes into the array below
    with open(infile1, 'r') as f:
        df = pd.read_csv(f)
        PDBs = list(df['PDB'])
        PDBs_sliced = [x[:-1] for x in PDBs] 
        
    with open(infile2, 'r') as f, open(outfile, 'w') as outfile:
        writer = csv.DictWriter(outfile, fieldnames = ["Unnamed", "IMGT entry ID", "IMGT molecule name", "Species",
                                                       "IMGT entry type", "IMGT receptor description", "Ligand(s)",
                                                       "Experimental technique", "Resolution", "PDB release date"])
        writer.writeheader()
        for line in f:
            if str(line.split(',')[1]) in PDBs_sliced:
                outfile.write(line)

matchingPDB('IGHV1-69_01 Data.csv', 'filtered_IMGT.csv', 'PDBsmatched_V1-69_01.csv')
matchingPDB('IGHV1-69_02 Data.csv', 'filtered_IMGT.csv', 'PDBsmatched_V1-69_02.csv')
matchingPDB('IGHV1-8_01 Data.csv', 'filtered_IMGT.csv', 'PDBsmatched_V1-8_01.csv')
matchingPDB('IGHV1-46_01 Data.csv', 'filtered_IMGT.csv', 'PDBsmatched_V1-46_01.csv')
matchingPDB('IGHV3-21_01 Data.csv', 'filtered_IMGT.csv', 'PDBsmatched_V3-21_01.csv')


def download_file(infile):
    """Function below takes the PDB list from the earlier function of comparing PDBs with IMGT &
    PyIgClassify data and downloads the relevant 'List of atom-atom interactions across
    protein-protein interface' document for each PDB and saves as CSV file"""
    PDB_searchcodes = []
    df = pd.read_csv(infile, sep='\t', engine='python', delimiter = ',')
    for entry in (df['IMGT entry ID'].str.lower()):
        PDB_searchcodes = entry
        url = "http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb="+PDB_searchcodes + '&chain1=H&chain2=L'
        req = requests.get(url)
        filename = req.url[url.rfind('?')+5:-18]+'.csv'
        with open(filename, 'wb') as f:
            for chunk in req.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
            outfile_name = 'out-' + filename
            with open(filename,'r') as infile, open(outfile_name, 'w') as outfile:
                copy = False
                for line in infile:
                    #if line in "<":
                        #infile = line.replace(char,'')
                    if line.strip() == "Hydrogen bonds":
                        copy = True
                        continue
                    elif line.strip() == "Non-bonded contacts":
                        copy = False
                        continue
                    elif copy:
                        outfile.write(line)
                

download_file('PDBsmatched_V1-69_01.csv')
download_file('PDBsmatched_V1-69_02.csv')
download_file('PDBsmatched_V1-8_01.csv')
download_file('PDBsmatched_V1-46_01.csv')
download_file('PDBsmatched_V3-21_01.csv')

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
    with open('output.csv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        for row in new_rows.items():
            csvwriter.writerow(row)


"""IGHV1-69*01"""
parsecsv('out-2cmr.csv')
parsecsv('out-2fx7.csv')
parsecsv('out-2fx8.csv')

"""IGHV1-69*02"""
parsecsv('out-1yym.csv')
parsecsv('out-2dd8.csv')

"""IGHV1-8*01"""
parsecsv('out-3x3f.csv')
parsecsv('out-4xmp.csv')

"""IGHV1-46*01"""
parsecsv('out-4lsp.csv')
parsecsv('out-5f9o.csv')

"""IGHV3-21*01"""
parsecsv('out-6cxg.csv')
parsecsv('out-6mu3.csv')

