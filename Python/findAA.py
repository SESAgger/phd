#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import argparse
import os
import time
import sys
import csv
import pandas as pd 
import numpy as np
from sigfig import round

parser=argparse.ArgumentParser(description="Finds AA and TSTV for all positions.")
parser.add_argument("-i","--input_file",help="sample file, must be tsv, path can be included (/proj/../fil.txt) ",default="samples.tsv")
parser.add_argument("-n","--name",help="What you want the output to be called", default="AA")
parser.add_argument("-p","--phylop",help="phylopfile")
parser.add_argument("-s","--sift",help="siftfile")


args = parser.parse_args()
aa = pd.read_csv(args.input_file,sep='\t')
number_of_canids=len(aa.columns[aa.columns.str.find("GT")!=-1])

aa.columns=aa.columns.str.lstrip(" # [1234567890]").str.replace(":GT","")

#Define columns containing canids 
canids=aa.columns[4:4+number_of_canids]

#Count number of each category. Could instead count the ones that aren't the 3, but I think this is fine
aa["homref"] = aa[aa[canids] == "0/0"].count(axis = 1)
aa["het"]    = aa[aa[canids] == "0/1"].count(axis = 1)
aa["homalt"] = aa[aa[canids] == "1/1"].count(axis = 1)

#Infer the ancestral allele
#If the dhole is {REF|ALT} and either Adu or Mes is {REF|ALT} and at least half of the rest of the canids are {REF|ALT} the {REF|ALT} is the ancestral allele
#half is 13: 24 canids in total, we know at least 2 are {REF|ALT} (first two requirements), so half of the rest is 11 (22/2) and then we add the 2 back=13
conditions = [
    (aa.BerlinZoo == "0/0") & 
        ((aa.C_adustus == "0/0") | (aa.C_mesomelas == "0/0")) & 
        (aa.homref >= (number_of_canids-2)/2+2),
    (aa.BerlinZoo == "1/1") & 
        ((aa.C_adustus == "1/1") | (aa.C_mesomelas == "1/1")) & 
        (aa.homalt >= (number_of_canids-2)/2+2),
]
#if above is REF, then anc=ref and if above is alt, then anc=alt 
anc_all = [
    aa.REF, aa.ALT
]
#if above is REF, then der=alt and if above is alt, then der=ref
der_all = [
    aa.ALT, aa.REF
]

# make the columns
aa['AA']  = np.select(conditions, anc_all, default = np.nan)
aa['DER'] = np.select(conditions, der_all, default = np.nan)

#Determine if it is a transversion or not
tstv = [
    (((aa['REF']=='A')|(aa['REF']=='G'))&
        ((aa['ALT']=='C')|(aa['ALT']=='T'))
    )|
    (
        ((aa['ALT']=='A')|(aa['ALT']=='G'))&
        ((aa['REF']=='C')|(aa['REF']=='T'))
    ),
    (
        ((aa['REF']=='A')|(aa['REF']=='G'))&
        ((aa['ALT']=='A')|(aa['ALT']=='G'))
    )|
    (
        ((aa['ALT']=='C')|(aa['ALT']=='T'))&
        ((aa['REF']=='C')|(aa['REF']=='T'))
    )
]

#Make the column
aa['Type'] = np.select(tstv,choicelist = ['V','S'],default=np.nan)

#Import the PhyloP scores
phylop = pd.read_csv(args.phylop, sep = '\t', usecols = [0,2,4], names = ['CHROM','POS','PhyloP'], header = 0)

#Make reference file (should also contain phastcon and sift at some point)
ref = aa.merge(phylop, how = "left")

sift = pd.read_csv(args.sift,sep="\t",comment="#",usecols=[1,6,13],names=['Location','Consequence','Extra'])

sift[['CHROM','POS']]=sift.Location.str.split(":",expand=True)
sift['POS']=sift.POS.astype(int)
sift['SIFT_txt']=sift.Extra.str.extract(r'SIFT=(\w+)')
sift['SIFT_score']=sift.Extra.str.extract(r'SIFT=.*(\d+\.\d*)').astype(float)
sift.drop(columns=["Extra","Location"])

ref1=ref.merge(sift,how="outer")

#Make reference file
ref1.to_csv(args.name+"_annotated_AA_TSTV_PhyloP_SIFT.tsv",sep="\t",columns=["CHROM","POS","AA","DER","Type","PhyloP","SIFT_score","SIFT_txt","Consequence"],index=False)

ref1[(ref1["AA"]!="")&(ref1["Type"]=="V")&(ref1.SIFT_txt.notna())].to_csv(args.name+"_der_transversions_sift_pp.tsv",sep="\t",columns=["CHROM","POS","DER","Type","PhyloP","SIFT_score","SIFT_txt","Consequence"],index=False)