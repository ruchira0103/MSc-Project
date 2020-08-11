#!/usr/bin/env python3

"""
-----------------------------------------------------------------------------------------------
File:		biodata.py
Program:	Biodata
Version:	3.0
Created:	28 May 2020
Function:	This code takes the DNA sequences and converts to protein sequences
Author: 	Ruchira Sachdeva 
------------------------------------------------------------------------------------------------
"""
#***********************************************************************************************
#Import Libraries

import Bio
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_rna, generic_protein
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import SeqIO

#***********************************************************************************************#
"""Trial Code"""
##for record in SeqIO.parse('IGHV.fasta','fasta'):
##    IGHV1_id = str(record.id) + '\n'
##    IGHV1_seq = str(record.seq.ungap("."))+ '\n'
##    IGHV1_mrna = str(record.seq.ungap(".").transcribe().upper())
##print(IGHV1_id)
##print(IGHV1_seq)
##print(IGHV1_mrna)

"""Trial Code"""
##with open ("outfile.fasta", "r") as f:
##    data = f.readlines()
##    for record in SeqIO.parse('outfile.fasta','fasta'):
##        IGHV1_seq = str(record.seq.translate()) + '\n'
##        #print(IGHV1_seq)

def DNA_to_mRNA(infile):
    """Function below allows the transcription of DNA sequence to mRNA sequence"""
    with open("mRNA.fasta","w") as f:
        for record in SeqIO.parse(infile, "fasta"):
            f.write(">" + str(record.id) + "\n")
            f.write(str(record.seq.ungap(".").transcribe().upper()) + "\n")
      
def mRNA_to_protein(filename):
    """Function below allows the translation of mRNA sequence to protein sequence""" 
    """to use to create logoplots later"""
    with open("protein.fasta","w") as f:
        for record in SeqIO.parse(filename, "fasta"):
            f.write(">" + str(record.id) + "\n")
            f.write(str(record.seq.translate()) + "\n")

DNA_to_mRNA('IGHV.fasta')
mRNA_to_protein("mRNA.fasta")

"""Codon Table Trial Code"""
#standard_table = CodonTable.unambiguous_rna_by_name['Standard']
#print(standard_table)

