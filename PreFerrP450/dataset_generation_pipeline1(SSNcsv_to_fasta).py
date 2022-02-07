# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np

cytoscape_ferredoxin_cluster_mapping=pd.read_csv("neighbouring=4_threshold=20/whole_dataset_fastatablep450ferredoxin_neighbouring_progenome2_neighbouring_4_dereplicated_ferredoxins_threshold_20_fasta Full Network colorized default node.csv")
ferredoxin_p450_mapping=pd.read_csv("whole_dataset_fastatablep450ferredoxin_neighbouring_progenome2_neighbouring=4_deduplicated.csv")
filename_output_scaffold="neighbouring=4_threshold=20/whole_dataset_neighbouring_4_threshold_20_deduplicated_"
reference_sequence=SeqRecord(Seq("MSAVALPRVSGGHDEHGHLEEFRTDPIGLMQRVRDECGDVGTFQLAGKQVVLLSGSHANEFFFRAGDDDLDQAKAYPFMTPIFGEGVVFDASPERRKEMLHNAALRGEQMKGHAATIEDQVRRMIADWGEAGEIDLLDFFAELTIYTSSACLIGKKFRDQLDGRFAKLYHELERGTDPLAYVDPYLPIESLRRRDEARNGLVALVADIMNGRIANPPTDKSDRDMLDVLIAVKAETGTPRFSADEITGMFISMMFAGHHTSSGTASWTLIELMRHRDAYAAVIDELDELYGDGRSVSFHALRQIPQLENVLKETLRLHPPLIILMRVAKGEFEVQGHRIHEGDLVAASPAISNRIPEDFPDPHDFVPARYEQPRQEDLLNRWTWIPFGAGRHRCVGAAFAIMQIKAIFSVLLREYEFEMAQPPESYRNDHSKMVVQLAQPACVRYRRRTGV"), id="Reference_P450")
ferredoxin_p450_mapping["Ferredoxin_cluster"]=np.nan
for index_f, row_ferre in ferredoxin_p450_mapping.iterrows():
    for index_c, row_cyto in cytoscape_ferredoxin_cluster_mapping.iterrows():
        if row_ferre["ferredoxin_id"] in row_cyto["Description"]:
            ferredoxin_p450_mapping["Ferredoxin_cluster"][index_f]=row_cyto["Node Count Cluster Number"]
            break
    print (index_f)
ferredoxin_p450_mapping.to_csv("neighbouring=4_threshold=20/whole_dataset_cytoscape_mapping_neighbouring_progenome2_neighbouring=4_threshold=20.csv")


for cluster in range(1,50):
    fasta=[]
    filename_output=filename_output_scaffold+str(cluster)+".fasta"
    print (filename_output)
    for index, row in ferredoxin_p450_mapping.iterrows():
        if row["Ferredoxin_cluster"]==cluster and "ferredoxin" not in row["p450_description"] and "Ferredoxin" not in row["p450_description"]:
            
            fasta.append( SeqRecord(Seq(row["p450_sequence"]),id=row["p450_id"]))
           
    fasta.append(reference_sequence)
    SeqIO.write(fasta, filename_output, 'fasta')
    


