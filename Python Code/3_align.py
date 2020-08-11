#!/usr/bin/env python3
"""
-----------------------------------------------------------------------------------------------
File:	    align.py
Program:	align
Version:	1.0
Created:	11 Aug 2020
Function:	Code below allows alignment of fasta files 
Author: 	Ruchira Sachdeva 
------------------------------------------------------------------------------------------------
"""
#***********************************************************************************************
#Import Libraries

from Bio import AlignIO
from Bio import SeqIO
from Bio import Seq

#***********************************************************************************************#

def alignfasta(input_file):
    """Function below uses a fasta file to align sequences"""
    records = SeqIO.parse(input_file, 'fasta')
    records = list(records) # make a copy, otherwise our generator
                            # is exhausted after calculating maxlen
    maxlen = max(len(record.seq) for record in records)

    # pad sequences so that they all have the same length
    for record in records:
        if len(record.seq) != maxlen:
            sequence = str(record.seq).ljust(maxlen, '.')
            record.seq = Seq.Seq(sequence)
    assert all(len(record.seq) == maxlen for record in records)

    # write to temporary file and do alignment
    output_file = '{}_aligned.fasta'.format(os.path.splitext(input_file)[0])
    with open(output_file, 'w') as f:
        SeqIO.write(records, f, 'fasta')
    alignment = AlignIO.read(output_file, "fasta")
    #print (alignment) 

alignfasta('mRNA_ungapped.fasta')

