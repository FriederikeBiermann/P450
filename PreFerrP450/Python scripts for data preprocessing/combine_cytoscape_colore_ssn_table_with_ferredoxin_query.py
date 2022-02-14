# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import numpy as np

cytoscape_ferredoxin_cluster_mapping=pd.read_csv("whole_dataset_neighbouring_ferredoxins_progenome2_treshold_30_50-125_fasta Full Network colorized default node.csv")
ferredoxin_p450_mapping=pd.read_csv("whole_dataset_neighbouring_ferredoxins_progenome2.csv")
ferredoxin_p450_mapping["Ferredoxin_cluster"]=np.nan
for index_f, row_ferre in ferredoxin_p450_mapping.iterrows():
    for index_c, row_cyto in cytoscape_ferredoxin_cluster_mapping.iterrows():
        if row_ferre["ferredoxin_id"] in row_cyto["Description"]:
            ferredoxin_p450_mapping["Ferredoxin_cluster"][index_f]=row_cyto["Node Count Cluster Number"]
ferredoxin_p450_mapping.to_csv("cytoscape_ferredoxin_p450_mapping_neighbouring_ferredoxins_threshold=30.csv")
