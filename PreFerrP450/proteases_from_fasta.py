#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 15:34:28 2022

@author: friederike
#a tool to find all neighbouring P450/ferredoxins in fasta from ProGenome2 database
"""


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import os
import glob

#fill in correct filenames for input/output files
path = '/home/friederike/Documents/Databases/TryptorubinProducersGenomes/RawData/'

output_proteases="proteases_tryptorubin_producers.fasta"
proteases=[]
#process antismash file
for progenome_file in glob.glob(os.path.join(path, '*.fasta')):
        #for each genome
        for progenome_record in SeqIO.parse(progenome_file, "fasta"):
            if "protease" in progenome_record.description or"Protease" in progenome_record.description :
                proteases.append(progenome_record)
#safe files
SeqIO.write(proteases, output_proteases, 'fasta')
