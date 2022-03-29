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


parser=argparse.ArgumentParser(description="Calculates mutational load")
parser.add_argument("-s","--sample_file",help="tsv with CHROM, POS, TGT")
parser.add_argument("-r","--reference_file",
help="What's your reference")
parser.add_argument("-n","--name", help="Name for produced files.",default="ny")
#argcomplete.autocomplete(parser)

args = parser.parse_args()



#Import the specific canids genotype info
canid = pd.read_csv(args.sample_file, sep = '\t', names = ["CHROM","POS","TGT"])


#Import the reference file
refa = pd.read_csv(args.reference_file, sep='\t',dtype = {'CHROM': object, 'POS': int, 'AA': object, 'DER': object, 'Type': object, 'PhyloP': float, 'SIFT_txt': object, 'SIFT_score': float, 'Consequence': object })

# Combine the 2 dataframes
canid_for_calc=canid.merge(refa, how = "left")

# Figure out if the GT is homozygous derived allel or if either PhyloP, sift or derived allele is unknown 
der_and_tv = (canid_for_calc["TGT"] == canid_for_calc["DER"] + "/" + canid_for_calc["DER"])

## Across genome
phylo_score_hom_tv = canid_for_calc[der_and_tv].PhyloP.sum()
sift_score_hom_tv = canid_for_calc[der_and_tv].SIFT_score.sum()
#phastcon_score_hom_tv = canid_for_calc[der_and_tv].PhastCon.sum()

## Genic
phylo_score_hom_tv_genic = canid_for_calc[(der_and_tv)&(canid_for_calc["Consequence"]!="intergenic_variant")].PhyloP.sum()
sift_score_hom_tv_genic = canid_for_calc[(der_and_tv)&(canid_for_calc["Consequence"]!="intergenic_variant")].SIFT_score.sum()
#phastcon_score_hom_tv = canid_for_calc[(der_and_tv)&(canid_for_calc["Consequence"]!="intergenic_variant")].PhastCon.sum()

## Non-genic
phylo_score_hom_tv_nongenic = canid_for_calc[(der_and_tv)&(canid_for_calc["Consequence"]=="intergenic_variant")].PhyloP.sum()
sift_score_hom_tv_nongenic = canid_for_calc[(der_and_tv)&(canid_for_calc["Consequence"]=="intergenic_variant")].SIFT_score.sum()
#phastcon_score_hom_tv = canid_for_calc[(der_and_tv)&(canid_for_calc["Consequence"]=="intergenic_variant")].PhastCon.sum()

# Generals

## Genic and nongenic SNPs
intergenic_variants=len(canid_for_calc[(der_and_tv)&(canid_for_calc["Consequence"]=="intergenic_variant")])
genic_variants=len(canid_for_calc[(der_and_tv)&(canid_for_calc["Consequence"]!="intergenic_variant")])


## Number of locations where a transversion has happened and gt=the derived allele
hom_der=len(canid_for_calc[(canid_for_calc["Type"]=="V")&(canid_for_calc["TGT"] == canid_for_calc["DER"] + "/" + canid_for_calc["DER"])].index)
hom_der_genic=len(canid_for_calc[(canid_for_calc["Type"]=="V")&(canid_for_calc["TGT"] == canid_for_calc["DER"] + "/" + canid_for_calc["DER"])&(canid_for_calc["Consequence"]=="intergenic_variant")].index)
hom_der_nongenic=len(canid_for_calc[(canid_for_calc["Type"]=="V")&(canid_for_calc["TGT"] == canid_for_calc["DER"] + "/" + canid_for_calc["DER"])&(canid_for_calc["Consequence"]!="intergenic_variant")].index)

## Number of locations where gt=the ancestral allele
hom_anc=len(canid_for_calc[(canid_for_calc["Type"]=="V")&(canid_for_calc["TGT"] == canid_for_calc["AA"] + "/" + canid_for_calc["AA"])].index)
hom_anc_genic=len(canid_for_calc[(canid_for_calc["Type"]=="V")&(canid_for_calc["TGT"] == canid_for_calc["AA"] + "/" + canid_for_calc["AA"])&(canid_for_calc["Consequence"]=="intergenic_variant")].index)
hom_anc_nongenic=len(canid_for_calc[(canid_for_calc["Type"]=="V")&(canid_for_calc["TGT"] == canid_for_calc["AA"] + "/" + canid_for_calc["AA"])&(canid_for_calc["Consequence"]!="intergenic_variant")].index)

## Results
### All
phylop_mutational_load=phylo_score_score_hom_tv/(hom_der+hom_anc)
sift_mutational_load=sift_score_score_hom_tv/(hom_der+hom_anc)

### Genic
phylop_mutational_load_genic=phylo_score_hom_tv_genic/(hom_der_genic+hom_anc_genic)
sift_mutational_load_genic=sift_score_hom_tv_genic/(hom_der+hom_anc)

### Non-genic
phylop_mutational_load_nongenic=phylo_score_hom_tv_nongenic/(hom_der_nongenic+hom_anc_nongenic)
sift_mutational_load_nongenic=sift_score_hom_tv_nongenic/(hom_der_nongenic+hom_anc_nongenic)

#Summary
w=open(args.name+"mutational_load_summary.txt",'w')
w.write(
f'''#Summary from processing of test
#Can be read as tsv with comment="#" and sep="\\t"
## Background info
Number of non-genic variants:\t {intergenic_variants} 
Number of non-genic variants:\t {genic_variants} 

Number of derived alleles:\t {hom_der}
Number of genic derived alleles:\t {hom_der_genic}
Number of nongenic derived alleles:\t {hom_der_nongenic}
Number of ancestral alleles:\t {hom_anc}
Number of genic ancestral alleles:\t {hom_anc_genic}
Number of nongenic ancestral alleles:\t {hom_anc_genic}

## Results
### All positions
PhyloP mutational load:\t  {round(phylop_mutational_load,sigfigs=2)}  \n
Sift mutational load:\t  {round(sift_mutational_load,sigfigs=2)}  \n
    
### Genic positions
Genic PhyloP mutational load:\t  {round(phylop_mutational_load_genic,sigfigs=2)}
Genic Sift mutational load:\t  {round(sift_mutational_load_genic,sigfigs=2)}
    
### Nongenic positions
Nongenic PhyloP mutational load:\t  {round(phylop_mutational_load_nongenic,sigfigs=2)}
Nongenic Sift mutational load:\t  {round(sift_mutational_load_nongenic,sigfigs=2)}
'''
)
w.close()

