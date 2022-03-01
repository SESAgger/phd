#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import argparse
import os
import time
import sys
import csv
import pandas as pd 
import numpy as np

parser=argparse.ArgumentParser(description="Finds AA and TSTV for all positions.")
parser.add_argument("-s","--sample_file",help="sample file, can be in any format, path can be included (/proj/../fil.txt) ",default="samples.txt")
parser.add_argument("-n","--name",help="What you want the output to be called", default="AA")

args = parser.parse_args()
df=pd.read_csv(args.sample_file,sep=',')

canids=df.columns[4:29]

df["homref"]=df[df[canids]=="0/0"].count(axis=1)
df["het"]=df[df[canids]=="0/1"].count(axis=1)
df["homalt"]=df[df[canids]=="1/1"].count(axis=1)


conditions = [
    (df.BerlinZoo=="0/0")&(df.C_adustus=="0/0")|(df.C_mesomelas=="0/0")&(df.homref >=10),
    (df.BerlinZoo=="1/1")&(df.C_adustus=="1/1")|(df.C_mesomelas=="1/1")&(df.homalt >=10),
]
values = [
    df.REF,df.ALT
]


df['AA']=np.select(conditions,values,default="NA")

condlist=[
    (((df['REF']=='A')|(df['REF']=='G'))&
        ((df['ALT']=='C')|(df['ALT']=='T'))
    )|
    (
        ((df['ALT']=='A')|(df['ALT']=='G'))&
        ((df['REF']=='C')|(df['REF']=='T'))
    ),
    (
        ((df['REF']=='A')|(df['REF']=='G'))&
        ((df['ALT']=='A')|(df['ALT']=='G'))
    )|
    (
        ((df['ALT']=='C')|(df['ALT']=='T'))&
        ((df['REF']=='C')|(df['REF']=='T'))
    )
]

df['Type'] = np.select(condlist,choicelist=['V','S'],default="NA")

name=args.name+"_annotated_AA_TSTV.tsv"

df.to_csv(name,sep="\t",columns=["#CHROM","POS","AA","Type"],index=False)
