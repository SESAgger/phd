#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
import os
import sys

parser=argparse.ArgumentParser(description="Makes sbatch files based on a sample file (-s) and a template file (-t)")
parser.add_argument("-s","--sample_file",help="sample file, can be in any format",default="samples.txt")
parser.add_argument("-t","--template_file",help="template file, can be in any format",default="sbatch.txt")
parser.add_argument("-n","--name", help="Name for produced files.",default="ny")
argcomplete.autocomplete(parser)

args = parser.parse_args()

i=0
j=0
k=0
filer=[]


lfil=open(args.sample_file,'r')
liste=lfil.readlines()

tfil=open(args.template_file)
template=tfil.readlines()

while i < len(liste):
	sti=liste[i].replace('\n','')
	i+=1
	isample=str(sti).find('FIL1')
	navn=str(sti)[isample+1:]
	filnavn=navn+'_'+str(args.name)+'.sh'
	pth = ''
	#print(navn+"done")
	nyfil=open(filnavn,"w")
	for j in range(len(template)):
		nyline=template[j].replace('FIL1',navn)
		nyline=nyline.replace('TID','5:0:0')
		nyfil.write(nyline)
		j+=1
		#print(nyline)
	nyfil.close()
	filer.append(filnavn)
w=open('all_'+str(args.name)+'.sh','w')
w.write('''#!/bin/bash -l 

''')

print("Made:")
while k<len(filer):
	if str(args.name) in filer[k]:
		w.write('sbatch ' + filer[k] + '\n \n')
		print(filer[k])
	k+=1
w.close()

print("\n"+str(args.sample_file) +" is used as the sample file\n"
+ str(args.template_file) +" is used as the template file\n"
+"files will be named sample_name_"+str(args.name)+".sh")