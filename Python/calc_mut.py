#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import argparse
import os
import time
import sys
import csv
import pandas as pd 
import numpy as np
import io
import subprocess
import sys
import time

start = time.time()


parser=argparse.ArgumentParser(description="Calculates mutational load")
parser.add_argument("-s","--sample_file",help="vcf file for bcftools")
parser.add_argument("-r","--reference_file",
help="What's your reference")
parser.add_argument("-n","--name", help="Name for produced files.",default="ny")
#argcomplete.autocomplete(parser)

args = parser.parse_args()


# Input
## Get all animals directly from stdout
#Get all animals directly from stdout
input_file=subprocess.run(['bcftools','query','-l','/Users/jqc305/Downloads/chr38_outgroup.vcf.gz'],stdout=subprocess.PIPE)
ids=input_file.stdout.decode().split("\n")
print("Canid IDs imported after: "+str(round(time.process_time(),2))+"s")

## Import the reference file
refa = pd.read_csv(args.reference_file, sep='\t',dtype = {'CHROM': object, 'POS': int, 'AA': object, 'DER': object, 'Type': object, 'PhyloP': float, 'SIFT_txt': object, 'SIFT_score': float, 'Consequence': object })
print("Reference imported after: "+str(round(time.process_time(),2))+"s")


# Run the pipeline
i=0
t=pd.DataFrame([])
while i < len(ids):
    #Get the animals directly from stdout one at a time to avoid ridiculus memory use
    input_file=subprocess.run(['bcftools','query','-f%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n',"/Users/jqc305/Downloads/chr38_outgroup.vcf.gz",'-H','-s',ids[i]],stdout=subprocess.PIPE)
    data = io.StringIO(input_file.stdout.decode())
    canids=pd.read_csv(data,sep="\t")
    canids.columns=canids.columns.str.lstrip(" # [1234567890]").str.replace(":GT","")
    # Combine the 2 dataframes
    canids_for_calc=canids.merge(refa, how = "left")

    print("DF combined after: "+str(round(time.process_time(),2))+"s")
    # Variables
    w_pp_ph = (canids_for_calc.PhastCon.notna())&(canids_for_calc.PhyloP.notna())
    w_pp_ph_sift = (canids_for_calc.PhastCon.notna())&(canids_for_calc.PhyloP.notna())&(canids_for_calc.SIFT_score)
    intergen = (canids_for_calc["Consequence"]=="intergenic_variant")
    gen = (canids_for_calc["Consequence"]!="intergenic_variant")
    anca = (canids_for_calc[ids[i]] == canids_for_calc["AA"] + "/" + canids_for_calc["AA"])
    
    ## Figure out if the GT is homozygous derived allel or if either PhyloP, sift or derived allele is unknown 
    der_and_tv = (canids_for_calc[ids[i]] == canids_for_calc["DER"] + "/" + canids_for_calc["DER"])

    # Conservation scores
    ## Across genome
    phylo_score_hom_tv = canids_for_calc[der_and_tv].PhyloP.sum()
    sift_score_hom_tv = canids_for_calc[der_and_tv].SIFT_score.sum()
    phastcon_score_hom_tv = canids_for_calc[der_and_tv].PhastCon.sum()

    ## Genic
    phylo_score_hom_tv_genic = canids_for_calc[der_and_tv&gen].PhyloP.sum()
    sift_score_hom_tv_genic = canids_for_calc[der_and_tv&gen].SIFT_score.sum()
    phastcon_score_hom_tv_genic = canids_for_calc[der_and_tv&gen].PhastCon.sum()
    
    ## Non-genic
    phylo_score_hom_tv_nongenic = canids_for_calc[der_and_tv&intergen].PhyloP.sum()
    sift_score_hom_tv_nongenic = canids_for_calc[der_and_tv&intergen].SIFT_score.sum()
    phastcon_score_hom_tv_nongenic = canids_for_calc[der_and_tv&intergen].PhastCon.sum()    
    
    # General results

    ## Genic and nongenic SNPs
    variants=len(canids_for_calc[der_and_tv])
    genic_variants=len(canids_for_calc[der_and_tv&gen])
    intergenic_variants=len(canids_for_calc[der_and_tv&intergen])


    ## Denominators
    ### In total
    hom_der=len(canids_for_calc[(canids_for_calc["Type"]=="V")&der_and_tv].index)
    hom_der_genic=len(canids_for_calc[(canids_for_calc["Type"]=="V")&der_and_tv&gen].index)
    hom_der_nongenic=len(canids_for_calc[(canids_for_calc["Type"]=="V")&(der_and_tv)&intergen].index)
    
    ### With PhyloP, and PhastCon
    hom_der_pp_ph=len(canids_for_calc[(canids_for_calc["Type"]=="V")&der_and_tv&w_pp_ph].index)
    hom_der_genic_pp_ph=len(canids_for_calc[(canids_for_calc["Type"]=="V")&der_and_tv&gen&w_pp_ph].index)
    hom_der_nongenic_pp_ph=len(canids_for_calc[(canids_for_calc["Type"]=="V")&der_and_tv&intergen&w_pp_ph].index)
    
    ### With SIFT, PhyloP, and PhastCon
    hom_der_s_pp_ph=len(canids_for_calc[(canids_for_calc["Type"]=="V")&der_and_tv&w_pp_ph_sift].index)
    hom_der_genic_s_pp_ph=len(canids_for_calc[(canids_for_calc["Type"]=="V")&der_and_tv&gen&w_pp_ph_sift].index)
    hom_der_nongenic_s_pp_ph=len(canids_for_calc[(canids_for_calc["Type"]=="V")&(der_and_tv)&intergen&w_pp_ph_sift].index)
    
    # Number of locations where gt=the ancestral allele
    ### In total
    hom_anc=len(canids_for_calc[anca].index)
    hom_anc_genic=len(canids_for_calc[anca&gen].index)
    hom_anc_nongenic=len(canids_for_calc[anca&intergen].index)
    
    ### With PhyloP, and PhastCon
    hom_anc_pp_ph=len(canids_for_calc[anca&w_pp_ph].index)
    hom_anc_genic_pp_ph=len(canids_for_calc[anca&gen&w_pp_ph].index)
    hom_anc_nongenic_pp_ph=len(canids_for_calc[anca&intergen&w_pp_ph].index)
    
    ### With SIFT, PhyloP, and PhastCon
    hom_anc_s_pp_ph=len(canids_for_calc[anca&w_pp_ph_sift].index)
    hom_anc_genic_s_pp_ph=len(canids_for_calc[anca&gen&w_pp_ph_sift].index)
    hom_anc_nongenic_s_pp_ph=len(canids_for_calc[anca&intergen&w_pp_ph_sift].index)
    
    ## Results
    ### All
    phylop_mutational_load=phylo_score_hom_tv/(hom_der_pp_ph+hom_anc_pp_ph)
    phast_mutational_load=phastcon_score_hom_tv/(hom_der_pp_ph+hom_anc_pp_ph)
    
    sift_mutational_load=sift_score_hom_tv/(hom_der_s_pp_ph+hom_anc_s_pp_ph)

    
    ### Genic
    phylop_mutational_load_genic=phylo_score_hom_tv_genic/(hom_der_genic_pp_ph+hom_anc_genic_pp_ph)
    phast_mutational_load_genic=phastcon_score_hom_tv_genic/(hom_der_genic_pp_ph+hom_anc_genic_pp_ph)

    sift_mutational_load_genic=sift_score_hom_tv_genic/(hom_der_nongenic_s_pp_ph+hom_anc_genic_s_pp_ph)

    
    ### Non-genic
    #print(phylo_score_hom_tv_nongenic)
    #print(hom_der_nongenic_pp_ph)
    phylop_mutational_load_nongenic=phylo_score_hom_tv_nongenic/(hom_der_nongenic_pp_ph+hom_anc_nongenic_pp_ph)
    phast_mutational_load_nongenic=phastcon_score_hom_tv_nongenic/(hom_der_nongenic_pp_ph+hom_anc_nongenic_pp_ph)
    if(hom_der_nongenic_s_pp_ph==0&hom_anc_nongenic_s_pp_ph==0):
        sift_mutational_load_nongenic=0
    else:
        sift_mutational_load_nongenic=sift_score_hom_tv_nongenic/(hom_der_nongenic_s_pp_ph+hom_anc_nongenic_s_pp_ph)

    
    #Make dataframe
    mutational_load = [[ids[i],phylop_mutational_load,sift_mutational_load,phast_mutational_load,
                        phylop_mutational_load_genic, sift_mutational_load_genic,phast_mutational_load_genic,
                        phylop_mutational_load_nongenic,sift_mutational_load_nongenic,phast_mutational_load_nongenic,
                        phylo_score_hom_tv,sift_score_hom_tv,phastcon_score_hom_tv,
                        phylo_score_hom_tv_genic,sift_score_hom_tv_genic,phastcon_score_hom_tv_genic,
                        phylo_score_hom_tv_nongenic,sift_score_hom_tv_nongenic,phastcon_score_hom_tv_nongenic,
                        hom_anc,hom_anc_genic,hom_anc_nongenic,
                        hom_der,hom_der_genic,hom_der_nongenic]
                      ]
     
    # Create the pandas DataFrame
    df = pd.DataFrame(mutational_load, columns = ["ID", 'PhyloP mutational load','Sift mutational load','PhastCon mutational load',
                                                  'Genic PhyloP mutational load','Genic Sift mutational load','Genic PhastCon mutational load',
                                                  'Nongenic PhyloP mutational load','Nongenic Sift mutational load','Nongenic PhastCon mutational load',
                                                  'Sum of PhyloP','Sum of Sift','Sum of PhastCon',
                                                  'Sum of PhyloP in genes','Sum of Sift in genes','Sum of PhastCon in genes',
                                                  'Sum of PhyloP outside genes','Sum of Sift outside genes','Sum of PhastCon outside genes',                        
                                                  'Ancestral alleles',"Genic Ancestral alleles",'Non-genic Ancestral alleles','Derived transversion','Genic derived transversions','Non-genic derived transversions'])
    t=pd.concat([t,df])

    print(str(ids[i])+" finished after: "+str(round(time.process_time(),2))+"s")
    i=i+1

t.to_csv(args.name+'mutational_load.tsv', sep = "\t", index = False,mode="w")


