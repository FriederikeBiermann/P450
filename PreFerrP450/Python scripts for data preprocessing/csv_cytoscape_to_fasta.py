#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 10:18:15 2022
create fasta from csv for PreFerrP450
@author: friederike
"""
import Bio
import pandas as pd
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
cytoscape_file=pd.read_csv("openframes_tryptorubinlike_peptides_withoutduplicates_protein_fas Full Network colorized default node.csv")
filename_output_frame="openframes_tryptorubinlike_peptides_withoutduplicates_protein_fas Full Network colorized default node"

for cluster in range(1,150):
    filename_output=filename_output_frame+str(cluster)+".fasta"
    fasta=[]
    for index, row in cytoscape_file.iterrows():
     
        if row["Sequence Count Cluster Number"]==cluster:
            fasta.append( SeqRecord(Seq.Seq(row["Sequence"]),id=row["Description"]))
    SeqIO.write(fasta, filename_output, 'fasta')
    