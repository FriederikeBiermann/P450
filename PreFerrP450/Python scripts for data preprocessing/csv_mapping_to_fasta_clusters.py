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
cytoscape_file=pd.read_csv("whole_dataset_neighbouring_threshold_35/cytoscape_ferredoxin_p450_mapping_neighbouring_ferredoxins_threshold=35.csv")
filename_output_scaffold="whole_dataset_neighbouring_threshold_35/whole_dataset_threshold_35_"

for cluster in range(1,50):
    fasta=[]
    filename_output=filename_output_scaffold+str(cluster)+".fasta"
    print (filename_output)
    for index, row in cytoscape_file.iterrows():
        if row["Ferredoxin_cluster"]==cluster and "ferredoxin" not in row["p450_description"] and "Ferredoxin" not in row["p450_description"]:
            
            fasta.append( SeqRecord(Seq.Seq(row["p450_sequence"]),id=row["p450_id"]))
    
    SeqIO.write(fasta, filename_output, 'fasta')
    
