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
filename_output_scaffold="neighbouring=4_threshold=20/whole_dataset_neighbouring_4_threshold_20_deduplicated_ferredoxin_"
reference_sequence=SeqRecord(Seq("MSAVALPRVSGGHDEHGHLEEFRTDPIGLMQRVRDECGDVGTFQLAGKQVVLLSGSHANEFFFRAGDDDLDQAKAYPFMTPIFGEGVVFDASPERRKEMLHNAALRGEQMKGHAATIEDQVRRMIADWGEAGEIDLLDFFAELTIYTSSACLIGKKFRDQLDGRFAKLYHELERGTDPLAYVDPYLPIESLRRRDEARNGLVALVADIMNGRIANPPTDKSDRDMLDVLIAVKAETGTPRFSADEITGMFISMMFAGHHTSSGTASWTLIELMRHRDAYAAVIDELDELYGDGRSVSFHALRQIPQLENVLKETLRLHPPLIILMRVAKGEFEVQGHRIHEGDLVAASPAISNRIPEDFPDPHDFVPARYEQPRQEDLLNRWTWIPFGAGRHRCVGAAFAIMQIKAIFSVLLREYEFEMAQPPESYRNDHSKMVVQLAQPACVRYRRRTGV"), id="Reference_P450")


ferredoxin_p450_mapping=pd.read_csv("neighbouring=4_threshold=20/whole_dataset_cytoscape_mapping_neighbouring_progenome2_neighbouring=4_threshold=20.csv")


for cluster in range(1,11):
    fasta=[]
    filename_output=filename_output_scaffold+str(cluster)+".fasta"
    print (filename_output)
    for index, row in ferredoxin_p450_mapping.iterrows():
        if row["Ferredoxin_cluster"]==cluster and "ferredoxin" not in row["p450_description"] and "Ferredoxin" not in row["p450_description"]:
            
            fasta.append( SeqRecord(Seq(row["ferredoxin_sequence"]),id=row["ferredoxin_description"]))
           
    SeqIO.write(fasta, filename_output, 'fasta')
    


