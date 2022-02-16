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
cytoscape_file=pd.read_csv("proteases_tryptorubin_producers_lysine_proteases_2_threshold_5_fasta_cluster1(serine_peptidases).csv")
filename_output="proteases_tryptorubin_producers_lysine_proteases_2_threshold_5_fasta_cluster1(serine_peptidases)"
fasta=[]
for cluster in range(1,150):
    for index, row in cytoscape_file.iterrows():
        fasta.append( SeqRecord(Seq.Seq(row["Sequence"]),id=row["Description"]))
    SeqIO.write(fasta, filename_output, 'fasta')
    