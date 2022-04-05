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
parser.add_argument("-i","--input_file",help="sample file, must be tsv, path can be included (/proj/../fil.txt) ",default="samples.tsv")
parser.add_argument("-n","--name",help="What you want the output to be called", default="AA")
parser.add_argument("-p","--phylop",help="phylopfile")
parser.add_argument("-s","--sift",help="siftfile")
parser.add_argument("-t","--phast",help="phastfile")

args = parser.parse_args()

# Import files 
## Import the gts
aa = pd.read_csv(args.input_file,sep='\t')
## Import the PhyloP scores
phylop = pd.read_csv(args.phylop, sep = '\t', usecols = [0,2,4], names = ['CHROM','POS','PhyloP'], header = 0)
## Import the Sift scores
sift = pd.read_csv(args.sift,sep="\t",comment="#",usecols=[1,6,13],names=['Location','Consequence','Extra'])
sift[['CHROM','POS']]=sift.Location.str.split(":",expand=True)
sift['POS']=sift.POS.astype(int)
sift['SIFT_txt']=sift.Extra.str.extract(r'SIFT=(\w+)')
sift['SIFT_score']=sift.Extra.str.extract(r'SIFT=.*(\d+\.\d*)').astype(float)
sift.drop(columns=["Extra","Location"])

## Import the PhastCon scores
phast=pd.read_csv(args.phast,sep="\t",names=["CHROM","start","POS","PhastCon"])

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

#Make reference file
ref = aa.merge(phylop, how = "left")

ref1=ref.merge(sift,how="outer")

ref2=ref1.merge(phast,how="left")

#Make reference file
ref2.to_csv(args.name+"ref_file.tsv", sep = "\t", columns = ["CHROM","POS","AA","DER","Type","PhyloP","SIFT_txt","SIFT_score","Consequence","PhastCon"], index = False)

print("File is called: "+args.name+"_annotated_AA_TSTV_PhyloP_SIFT.tsv")

# Summary
## Basic
with_AA=ref1.AA.notna().sum()
tv=len(ref1[ref1["Type"]=="V"])
with_AA_tv=len(ref1[(ref1["Type"]=="V")&(ref1.AA.notna())])

## Values
with_PP=ref1.PhyloP.notna().sum()
with_sift=ref1.SIFT_score.notna().sum()
with_phastcon=ref1.PhastCon.notna.sum()

with_all_c=len(ref1[(ref1.PhyloP.notna())&(ref1.SIFT_score.notna())&(ref1.PhastCon.notna())])

## All
with_AA_tv_PP_sift_phastcon=len(ref1[(ref1.PhyloP.notna())&(ref1.SIFT_score.notna())&(ref1["Type"]=="V")&(ref1.AA.notna())&(ref1.PhastCon.notna())])

## Write summary
w=open(args.name+'reference_file_summary.txt','w')
w.write(
'#Summary from processing of ' + "test" +'\n'
'#Can be read as tsv with comment="#" and sep="\t"\n'
'Name\tNumber\t%of all\n'
'SNPs where an ancestral allele could be determined:\t' + str(with_AA) + '\t' + str(round(with_AA/len(ref1)*100,decimals=2)) + '%\n'
'SNPs with a transversion:\t' + str(tv) + '\t' + str(round(tv/len(ref1)*100,decimals=2)) + '%\n'
'Transversions with known ancestral allele:\t' + str(with_AA_tv)+ '\t' + str(round(with_AA_tv/len(ref1)*100,decimals=2)) + '%\n'
'SNPs with PhyloP score:\t' + str(with_PP) + '\t' + str(round(with_PP/len(ref1)*100,decimals=2)) + '%\n'
'SNPs with sift:\t' + str(with_sift) + '\t' + str(round(with_sift/len(ref1)*100,decimals=2)) + '%\n'
'SNPs with Phastcon score:\t' + str(round(with_phastcon,decimals=2)) + '\t' + str(round(with_phastcon/len(ref1)*100,decimals=2)) + '%\n'
'SNPs with all conservation scores:\t' +str(round(with_all_c,decimals=2)) + '\t' + str(round(with_all_c/len(ref1)*100,decimals=2)) + '%\n'
'Transversions with all information:\t' +str(round(with_AA_tv_PP_sift_phastcon,decimals=2)) + '\t' + str(round(with_AA_tv_PP_sift_phastcon/len(ref1)*100,decimals=2)) + '%\n'
)
w.close()