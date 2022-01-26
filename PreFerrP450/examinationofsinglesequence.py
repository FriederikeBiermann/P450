#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 18:24:41 2022
GUI
@author: friederike
"""


import pandas as pd
from Bio import pairwise2
import re
import pickle
import math
import itertools
import numpy as np
import pandas
from tkinter import *


def charge(sequence):
    AACharge = {"C":-.045,"D":-.999,  "E":-.998,"H":.091,"K":1,"R":1,"Y":-.001}
    charge = -0.002
    seqstr=str(sequence)
    seqlist=list(seqstr)
    for aa in seqlist:
        if aa in AACharge:
            charge += AACharge[aa]
    return charge
def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""
def easysequence (sequence):
    #creates a string out of the sequence file, that only states if AA is acidic (a), basic (b), polar (p), neutral/unpolar (n),aromatic (r),Cystein (s) or a Prolin (t)
    seqstr=str(sequence)
    seqlist=list(seqstr)
    easylist=[]
    for i in seqlist:
        if i == 'E' or i== 'D':           
            easylist=easylist+['a']
        if i == 'K' or i=='R' or i=='H':      
            easylist=easylist+['b']
        if i == 'S' or i=='T' or i=='N' or i=='Q':
            easylist=easylist+['p']
        if i == 'F' or i=='Y' or i=='W':
            easylist=easylist+['r']
        if i == 'C':
            easylist=easylist+['s']
        if i == 'P':
            easylist=easylist+['t']
        if i == 'G' or i=='A' or i=='V' or i=='L' or i=='I' or i=='M':
            easylist=easylist+['n']
            
    seperator=''
    easysequence=seperator.join(easylist)
    return easysequence
def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z        


#getting Data for p450 associated with terpenes
include_parts= False
def ferredoxin_calculation (seq_record):
    # fill in names of files here!
    foldernameoutput=""
    filename_permutations=foldernameoutput+"permutations.txt"
    classifier=foldernameoutput+'forestclassifieronlynext.sav'
    forest = pickle.load(open(classifier, 'rb')) 
    # creates tables
    testsequence=pd.DataFrame()
    
    # S#14703] cytochrome P450 51 Cyp51 ([P#10800] [CYP51] cytochrome P450, family 51 )  from https://cyped.biocatnet.de/sequence/14703
    alignmentfa=("MSAVALPRVSGGHDEHGHLEEFRTDPIGLMQRVRDECGDVGTFQLAGKQVVLLSGSHANEFFFRAGDDDLDQAKAYPFMTPIFGEGVVFDASPERRKEMLHNAALRGEQMKGHAATIEDQVRRMIADWGEAGEIDLLDFFAELTIYTSSACLIGKKFRDQLDGRFAKLYHELERGTDPLAYVDPYLPIESLRRRDEARNGLVALVADIMNGRIANPPTDKSDRDMLDVLIAVKAETGTPRFSADEITGMFISMMFAGHHTSSGTASWTLIELMRHRDAYAAVIDELDELYGDGRSVSFHLRQIPQLENVLKETLRLHPPLIILMRVAKGEFEVQGHRIHEGDLVAASPAISNRIPEDFPDPHDFVPARYEQPRQEDLLNRWTWIPFGAGRHRCVGAAFAIMQIKAIFSVLLREYEFEMAQPPESYRNDHSKMVVQLAQPACVRYRRRTGV")
    fragments=("begin","sbr1","sbr2","core","end","fes1","fes2")
    with open(filename_permutations, 'r') as f:
        permutations = [line.rstrip('\n') for line in f]
    for f in fragments:
        for i in permutations:
            n=i+f
            testsequence[n]=[]
    seq_record =input_box.get("1.0", "end")

    fewgaps = lambda x, y: -20 - y
    specificgaps = lambda x, y: (-2 - y)
    
    alignments = pairwise2.align.globalmc(alignmentfa,  seq_record, 1, -1, fewgaps, specificgaps)

    seqa=alignments[0][0]

    
    seqb=alignments[0][1]

    #begin
    for match in re.finditer("DAS",seqa):
                            end_begin=match.span()[1]
                          
    begbne=seqb[:end_begin+1].replace("-","")
    begb=easysequence(seqb[:end_begin+1].replace("-",""))

    #sbr1
    for match in re.finditer("ESL",seqa):
                         
            end_sbr1=match.span()[1]
    sbr1bne=seqb[end_begin:end_sbr1+1].replace("-","")
    sbr1b=easysequence(seqb[end_begin:end_sbr1+1].replace("-",""))

    #sbr2
    for match in re.finditer("MRH",seqa):
                     
            end_sbr2=match.span()[1]
    sbr2bne=seqb[end_sbr1:end_sbr2].replace("-","")
    sbr2b=easysequence(seqb[end_sbr1:end_sbr2].replace("-",""))

    #core
    for match in re.finditer("RCV",seqa):
        
            end_core=match.span()[1]
    corebne=seqb[end_sbr2:end_core+1].replace("-","")
    coreb=easysequence(seqb[end_sbr2:end_core+1].replace("-",""))
    
    endbne=seqb[end_core+1:].replace("-","")
    endb= easysequence(seqb[end_core+1:].replace("-",""))
    
    #fes1
    for match in re.finditer("GRH",seqa):

            end_fes1=match.span()[1]
    for match in re.finditer("FRA",seqa):
            begin_fes1= match.span()[0] 
    
    fes1bne=seqb[begin_fes1-1:end_fes1+1].replace("-","")
    fes1b=easysequence(seqb[begin_fes1-1:end_fes1+1].replace("-",""))
    
    
    #fes2
    for match in re.finditer("GRH",seqa):

            end_fes2=match.span()[1]
    for match in re.finditer("QED",seqa):
            begin_fes2= match.span()[0] 

    fes2bne=seqb[begin_fes2-1:end_fes2+1].replace("-","")
    fes2b=easysequence(seqb[begin_fes2-1:end_fes2+1].replace("-",""))
    
    listfragmentsalt=(begb,sbr1b,sbr2b,coreb,endb,fes1b,fes2b)
    nealt=(begbne,sbr1bne,sbr2bne,corebne,endbne,fes1bne,fes2bne)
    new_row={}
    c=0
    listfragments=[]
    ne=[]
    for frag in listfragmentsalt:
        if frag=='':
            listfragments.append("n")
        else:
            listfragments.append(frag)
    for frag in nealt:       
        if frag=='':
            ne.append("A")
        else:
            ne.append(frag)
    for frag in listfragments:
        f=fragments[c]
        for i in permutations:
            n=i+f
            new_row =merge_two_dicts(new_row,{n:frag.count(i)})
        a=f+"acidic"
        new_row =merge_two_dicts(new_row,{a:(frag.count("a")/len(frag))})
        aabs=f+"acidic absolute"
        new_row =merge_two_dicts(new_row,{aabs:(frag.count("a"))})
        char=f+"charge"
        new_row =merge_two_dicts(new_row,{char:charge(ne[c])})
        b=f+"basic"
        babs=f+"basic absolute"
        new_row =merge_two_dicts(new_row,{b:(frag.count("b")/len(frag))})
        new_row =merge_two_dicts(new_row,{babs:(frag.count("b"))})
        if include_parts==True:
            if f=="fes1" or f=="begin" or f== "end":
                how_many_parts=4
            else:
                how_many_parts=3
            each_part_length = math.ceil(len(frag)/how_many_parts)
            iterator = [iter(frag)] * each_part_length
            parts = list(itertools.zip_longest(*iterator, fillvalue='X'))
            for i in range (0,len(parts)):
                    part=f+"acidic"+ str(i)
                    par=parts[i]
                    new_row =merge_two_dicts(new_row,{part:(par.count("a")/(len(par)+1))})
                    part=f+"basic"+ str(i)
                    par=parts[i]
                    new_row =merge_two_dicts(new_row,{part:(par.count("b")/(len(par)+1))})
            if len(parts)<3:
                for dif in range ((len(parts)),4):
                    part=f+"acidic"+ str(dif)
                    new_row =merge_two_dicts(new_row,{part:0})
                    part=f+"basic"+ str(dif)
                    new_row =merge_two_dicts(new_row,{part:0})
        c=c+1  
    testsequence = testsequence.append(new_row, ignore_index=True)
    chargerows=[]
    acidicrows=[]
    basicrows=[]
    absacidicrows=[]
    absbasicrows=[]
    for f in fragments:
        chargerows.append(str(f)+"charge")
        acidicrows.append(str(f)+"acidic")
        basicrows.append(str(f)+"basic")
        absacidicrows.append(str(f)+"acidic absolute")
        absbasicrows.append(str(f)+"basic absolute")
    testsequence['complete charge']=testsequence[chargerows].sum(axis=1)
    testsequence['mean acidic']=testsequence[acidicrows].mean(axis=1)  
    testsequence['mean basic']=testsequence[basicrows].mean(axis=1)  
    testsequence['absolute acidic']=testsequence[absacidicrows].sum(axis=1)  
    testsequence['absolute basic']=testsequence[absbasicrows].sum(axis=1)
    testsequence = testsequence.replace('NaN', np.nan).fillna(0)
    
    foldernameoutput=""
    pathcompletetable=foldernameoutput+"tablep450onlyindex.csv"
    # fill in names of files here!
    completetable=pd.read_csv(pathcompletetable)
    completetable.drop(columns=['target'],inplace=True)

    completetable=testsequence.reindex(columns=completetable.columns.values).fillna(0)

    test_predictf=forest.predict(completetable)
    print (test_predictf)
    # prediciton probability
    y_score = forest.predict_proba(completetable)
    print (y_score)
    
    return messagebox.showinfo('Result', f'The matching family of ferredoxins is probably {test_predictf}' )

gui = Tk()
gui.title('PreFerrP450')
gui.geometry('500x500')

Label(gui, text='Enter amino acid sequence of CYP P450 & hit Enter Key').pack(pady=20)
input_box = Text(gui,height=2000, width= 450, yscrollcommand=TRUE )
input_box.bind('<Return>',ferredoxin_calculation)
input_box.pack()

gui.mainloop()
