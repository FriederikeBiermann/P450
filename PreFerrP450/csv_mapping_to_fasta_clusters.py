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
cytoscape_file=pd.read_csv("cytoscape_ferredoxin_p450_mapping.csv")
filename_output_scaffold="respresentatives_threshold_40/respresentatives_threshold_40_cluster_"

for cluster in range(1,50):
    fasta=[]
    filename_output=filename_output_scaffold+str(cluster)+".fasta"
    for index, row in cytoscape_file.iterrows():
        if row["Ferredoxin_cluster"]==cluster:
            
            fasta.append( SeqRecord(Seq.Seq(row["p450_sequence"]),id=row["p450_id"]))
    
    SeqIO.write(fasta, filename_output, 'fasta')
    