
MSc Project Code

Data availability
All data used in this study was download from the IMGT, PyIgClassify and PDBsum websites and can be found at: http://www.imgt.org/ and http://dunbrack2.fccc.edu/PyIgClassify/default.aspx and http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode=index.html respectively.

Code Availability
All custom scripts are written in Python 3 and are available open source in the MSc Project repository at https://github.com/ruchira0103/MSc-Project
The GitHub page is a collection of Python scripts intended to be used for this project to gain an understanding of the IGHV1-69 germline. The end goal is to facilitate the interpretation of potential germline variability.
The scripts do have some amount of overlap; however, it should run smoothly once everything is saved in the same folder when downloading the scripts.

Dependencies
It is all written in Python 3 – the scripts rely on Pandas, NumPy and matplotlib for visualization. 

Python Scripts
1. similarPDB.py Python script
    a.  Written in Python 3 was developed to understand the IGHV1-69 germlines  and its different alleles available, its CDR sequence frequencies and plot accordingly. For this         script, the modules pandas and matplotlib were heavily used. The data for this code was available from the PyIgClassify website.
    b.	The code takes in a CSV file and searches for similarity within the IGHV1-69 germline sequence and plots a histogram to showcase differences between the CDR1, CDR2 and           CDR3 sequence lengths.

2.	biodata.py Python script 
    a.	This code was developed to convert DNA sequences to protein sequences to be used to create logo plots later. This was done using the Biopython tool in Python, this               allowed the transcription (to mRNA sequence) and translation (to protein sequence) of the IGHV sequences for its 19 alleles. 
    b.	The code takes in one, or more, fasta file(s) as input and firstly converts the DNA sequence into an mRNA sequence and then into a protein sequence; saves the output             into the same folder as input file.

3.	align.py Python script 
    a.	This code was created to allow alignment of the mRNA (converted from part 3) sequences and display the results. 
    b.	Takes in one or more fasta file(s) as input and aligns the sequences and saves the output into the same folder as input file.

4.	filtercsv.py Python script
    a.	This code was created as an automated Python script to download the relevant information (using the PDB codes from the filtered list in the step below) from the PDBSum           website. Data is saved into text file(s) in the same folder as input file.
    b.	This code takes in one, or more, comma-separated file(s) as input and filters the PDBs with resolution of < 2.5A and saves to a new CSV file in the same folder as input         file. 
    c.	Then that list of filtered PDBs are used to compare from the IMGT file to PyIgClassify individual germline data files and creates a ‘new’ filtered list of PDB codes. 
    d.	The PDB codes are separately downloaded as their own individual files in CSV format to be used to parsing. 
    e.	The last function in this script, uses the individual CSV file and picks out specific columns to be displayed and at the same time creates a new column to capture the            amino acid residue symbols and identify the CDR1, 2, and 3 in the list. 

5.	CoV data.py Python script
    a.	This code was created to understand and display the distribution of the IGHV genes in SARS-CoV-2 outbreak.
    b.	This script takes in a comma-separated file(s) as input and creates a pivot table to count all the IGHV germlines separately and uses that to create the pie chart to             display the distribution. 
