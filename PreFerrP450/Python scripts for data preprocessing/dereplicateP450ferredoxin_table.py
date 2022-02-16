# -*- coding: utf-8 -*-
# *** Spyder Python Console History Log ***
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
original_mapping_table=pd.read_csv("whole_dataset_fastatablep450ferredoxin_neighbouring_progenome2_neighbouring=4.csv")
deduplicated_mapping_table=pd.DataFrame()
p450s=[]
ferredoxins=[]
list_p450s_ferredoxins=[]
for index,row in enumerate(original_mapping_table.iterrows()):
    identifier=[original_mapping_table["ferredoxin_sequence"][index],original_mapping_table["p450_sequence"][index]]
    if identifier in list_p450s_ferredoxins:
        print (identifier)
    else:
        list_p450s_ferredoxins+=[identifier]
        deduplicated_mapping_table=deduplicated_mapping_table.append(row[1], ignore_index=True)
        p450s+=[SeqRecord(Seq(original_mapping_table["p450_sequence"][index]),id= original_mapping_table["p450_id"][index])]
        ferredoxins+=[SeqRecord(Seq(original_mapping_table["ferredoxin_sequence"][index]),id= original_mapping_table["ferredoxin_id"][index])]
SeqIO.write(p450s, "whole_dataset_fastatablep450ferredoxin_neighbouring_progenome2_neighbouring=4_dereplicated_p450s.fasta", "fasta")
SeqIO.write(ferredoxins, "whole_dataset_fastatablep450ferredoxin_neighbouring_progenome2_neighbouring=4_dereplicated_ferredoxins.fasta", "fasta")
deduplicated_mapping_table.to_csv("whole_dataset_fastatablep450ferredoxin_neighbouring_progenome2_neighbouring=4_deduplicated.csv")