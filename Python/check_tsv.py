#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
import os
import sys
import os.path 
parser=argparse.ArgumentParser(description="Check tsv for sarek")
parser.add_argument("tsv",help="file to check",default="samples.tsv")
parser.add_argument("step",help="step the file is to be used for\n mapping - fastq-mapping\n dup_mar_tab for check output from prepare recalibration",default="mapping")


argcomplete.autocomplete(parser)
args = parser.parse_args()

lfil=open(args.tsv,'r')
liste=lfil.readlines()

step=args.step


size=7
pattern1="_R1"
pattern2="_R2"

if("bam" in step):
    size=6
    pattern1="bam"
    pattern2="bai"

if("dup_mar_tab" in step):
    pattern1="bam"
    pattern2="bai"
    pattern3="recal.table"
    size=7 

i=0
total_errors=0

nyfil=open("tsv_log.log","w")

while i < len(liste):
    errors=0
    SUCCESS=str("")
    ERRORS=str("")
    start=str("") 

    tsv=liste[i].replace("\n","").split("\t") 
    
    #Columns      
    if(len(tsv)==size):
        format=str("Format ok\n\n")
    elif (len(tsv)>size):
        print("You have too many columns\n\n:")
        exit()
    elif(len(tsv)<size):
        print("You're missing a column\n\n")
        exit()
    else:
        print("Something's wrong with the file format")
        exit()
    
    navn=tsv[0]
    if("dup_mar_tab" in step):
        file1=tsv[-3]
        file2=tsv[-2]
        file3=tsv[-1]
        name1=file1.split("/")[-1].split(".")[0]   
        name2=file2.split("/")[-1].split(".")[0]   
        name3=file3.split("/")[-1].split(".")[0]        
    else:
        file1=tsv[-2]
        file2=tsv[-1]
        name1=file1.split("/")[-1].split(".")[0]   
        name2=file2.split("/")[-1].split(".")[0]     

    start=str("Checking " + navn + " on line " + str(i+1) +"\n") 

    #Names
    
    if(name1.replace(pattern1,"")==name2.replace(pattern2,"")):
        if("dup_mar_tab" in step and name1.replace(pattern1,"")==name3.replace(pattern3,"") ):
            SUCCESS=SUCCESS+str("   Names \n     "+ file1 + ", \n     " + file2 + " and \n     " + file3 + " match\n")
        elif("dup_mar_tab" in step and name1.replace(pattern1,"")!=name3.replace(pattern3,"") ):
            ERRORS=ERRORS+str("   dup_mar_tab table name: " + name3 + " does not match \n      " + name1 +"\n")
        else:
                SUCCESS=SUCCESS+str("   Names\n     " + file1 + " and \n     " + file2 + " match\n")
    else:
        errors+=1
        ERRORS=ERRORS+str("   File names: " + name1 + " " + name2 + " do not match \n")

    #file1    
    if(pattern1 in file1):
        SUCCESS=SUCCESS+str("   "+file1 + " contains " + pattern1 + "\n")
    else:
        errors+=1
        ERRORS=ERRORS+str("   "+ file1 + " does not contain " +pattern1 +"\n") 
    
    if(os.path.isfile(file1)==True): 
        SUCCESS=SUCCESS+str("   " + file1 + " exist\n")    
    else: 
        errors+=1
        ERRORS=ERRORS+str("   " + file1 + " does not exist \n")
    
    #file 2
    if(pattern2 in file2):
        SUCCESS=SUCCESS+str("   " +file2 + " contains " + pattern2 + "\n")
    else:
        errors+=1
        ERRORS=ERRORS+str("   " + file2 + " does not contain " + pattern2 + "\n") 
        
    if os.path.isfile(file2)==True:
        SUCCESS=SUCCESS+str("   " + file2 + " exist\n") 
    else:
        errors+=1
        ERRORS=ERRORS+str("   "  + file2 + " does not exist \n")
  
    #file3
    if("dup_mar_tab" in step):
        
 
        if(pattern3 in file3):
             SUCCESS=SUCCESS+str("   " + file3 + " contains " + pattern3 + "\n")
        elif(pattern3 not in file3):
            errors+=1
            ERRORS=ERRORS+str("   " + file3 + " does not contain " + pattern3 + "\n")
            
        if os.path.isfile(file3)==True:
            SUCCESS=SUCCESS+str("   " + file3 + " exist\n") 
        else:
            errors+=1
            ERRORS=ERRORS+str("   "  + file3 + " does not exist \n")
            
    
    ERRORS=ERRORS+str(navn + " had \t " + str(errors) + " error(s) on line: "+str(i+1)+"\n")
    logs=start+SUCCESS+ERRORS+"\n"
    #print(logs)
    if errors!=0:
        print(start+ERRORS)
    
    total_errors=total_errors+errors
    i+=1
    nyfil.write(logs)
resultat=str("Total errors: " + str(total_errors)+"\n")
print(resultat)
nyfil.close()

print("See tsv_log.log for details")
