[toc]

# 160 dogs

## Full dataset

### $\bold{F_{ST}}$ for all

```bash
PLINK v1.90b4.9 64-bit (13 Oct 2017)
Options in effect:
  --allow-extra-chr
  --allow-no-sex
  --bfile all_by_bcf_use_this
  --dog
  --double-id
  --fst case-control
  --geno 0.05
  --keep phenos/pheno_fcr_vs_eriks.txt
  --out fcr_vs_eriks_fst_g5percent.fst
  --pheno phenos/pheno_fcr_vs_eriks.txt
  --within phenos/pheno_fcr_vs_eriks.txt

Hostname: rackham2.uppmax.uu.se
Working directory: /crex/proj/uppstore2017228/KLT.05.GRUS/sophieflatr/newrun/onlyvcf/SNP
Start time: Wed Jul  8 10:03:16 2020

Random number seed: 1594195396
257386 MB RAM detected; reserving 128693 MB for main workspace.
57459622 variants loaded from .bim file.
309 dogs (0 males, 0 females, 309 ambiguous) loaded from .fam.
Ambiguous sex IDs written to fcr_vs_eriks_fst_g5percent.fst.nosex .
179 phenotype values present after --pheno.
--keep: 179 dogs remaining.
--within: 2 clusters loaded, covering a total of 179 dogs.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 179 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate in remaining samples is 0.274524.
50541144 variants removed due to missing genotype data (--geno).
6918478 variants and 179 dogs pass filters and QC.
Among remaining phenotypes, 19 are cases and 160 are controls.
Writing --fst report (2 populations) to fcr_vs_eriks_fst_g5percent.fst.fst ...
done.
6711884 markers with valid Fst estimates (104551 excluded).
Mean Fst estimate: 0.117067
Weighted Fst estimate: 0.156497

End time: Wed Jul  8 10:03:51 2020
```





## Pruned dataset

###  Prune FCR 

```bash
PLINK v1.90b4.9 64-bit (13 Oct 2017)
Options in effect:
  --allow-extra-chr
  --allow-no-sex
  --bfile /proj/uppstore2017228/KLT.05.GRUS/sophieflatr/newrun/onlyvcf/SNP/all_by_bcf_use_this
  --dog
  --double-id
  --geno 0.05
  --indep-pairwise 50 10 0.5
  --out fcrpruned_pendleton_r2_0.05_geno_0.05
  --pheno /proj/uppstore2017228/KLT.05.GRUS/sophieflatr/newrun/onlyvcf/SNP/phenos/pheno_fcr_vs_eriks.txt
  --remove /proj/uppstore2017228/KLT.05.GRUS/sophieflatr/newrun/onlyvcf/SNP/phenos/3dogs.txt

Hostname: rackham2.uppmax.uu.se
Working directory: /crex/proj/uppstore2017228/KLT.05.GRUS/sophieflatr/newrun/onlyvcf/SNP
Start time: Fri Jul 24 10:10:09 2020

Random number seed: 1595578209
257386 MB RAM detected; reserving 128693 MB for main workspace.
Allocated 72389 MB successfully, after larger attempt(s) failed.
57459622 variants loaded from .bim file.
309 dogs (0 males, 0 females, 309 ambiguous) loaded from .fam.
Ambiguous sex IDs written to fcrpruned_pendleton_r2_0.05_geno_0.05.nosex .
179 phenotype values present after --pheno.
--remove: 306 dogs remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 306 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate in remaining samples is 0.535039.
50886053 variants removed due to missing genotype data (--geno).
6573569 variants and 306 dogs pass filters and QC.
Among remaining phenotypes, 19 are cases and 160 are controls.  (127 phenotypes
are missing.)
Pruned 254902 variants from chromosome 1, leaving 54156.
Pruned 161054 variants from chromosome 2, leaving 41020.
Pruned 218485 variants from chromosome 3, leaving 46448.
Pruned 207288 variants from chromosome 4, leaving 46243.
Pruned 199788 variants from chromosome 5, leaving 51672.
Pruned 170215 variants from chromosome 6, leaving 42265.
Pruned 174101 variants from chromosome 7, leaving 41023.
Pruned 182259 variants from chromosome 8, leaving 47365.
Pruned 110894 variants from chromosome 9, leaving 32979.
Pruned 137612 variants from chromosome 10, leaving 33068.
Pruned 146532 variants from chromosome 11, leaving 34772.
Pruned 183033 variants from chromosome 12, leaving 40786.
Pruned 156827 variants from chromosome 13, leaving 36749.
Pruned 145485 variants from chromosome 14, leaving 35336.
Pruned 150219 variants from chromosome 15, leaving 34392.
Pruned 151089 variants from chromosome 16, leaving 39196.
Pruned 159343 variants from chromosome 17, leaving 36864.
Pruned 117487 variants from chromosome 18, leaving 29991.
Pruned 157455 variants from chromosome 19, leaving 32825.
Pruned 113965 variants from chromosome 20, leaving 27977.
Pruned 140605 variants from chromosome 21, leaving 38656.
Pruned 145780 variants from chromosome 22, leaving 34273.
Pruned 130286 variants from chromosome 23, leaving 33586.
Pruned 103742 variants from chromosome 24, leaving 28056.
Pruned 118249 variants from chromosome 25, leaving 30157.
Pruned 102118 variants from chromosome 26, leaving 33378.
Pruned 114297 variants from chromosome 27, leaving 30304.
Pruned 92120 variants from chromosome 28, leaving 24518.
Pruned 114580 variants from chromosome 29, leaving 30502.
Pruned 85715 variants from chromosome 30, leaving 26751.
Pruned 106837 variants from chromosome 31, leaving 27562.
Pruned 118256 variants from chromosome 32, leaving 27722.
Pruned 78446 variants from chromosome 33, leaving 21889.
Pruned 115235 variants from chromosome 34, leaving 28697.
Pruned 78660 variants from chromosome 35, leaving 26689.
Pruned 69258 variants from chromosome 36, leaving 23374.
Pruned 68401 variants from chromosome 37, leaving 23286.
Pruned 69895 variants from chromosome 38, leaving 25913.
Pruned 95020 variants from chromosome 39, leaving 27470.
Pruned 104 variants from chromosome 42, leaving 22.
Pruning complete.  5245637 of 6573569 variants removed.
Marker lists written to fcrpruned_pendleton_r2_0.05_geno_0.05.prune.in and
fcrpruned_pendleton_r2_0.05_geno_0.05.prune.out .

End time: Fri Jul 24 10:11:15 2020
```

### $\bold{F_{ST}}$ for pruned 

```bash
PLINK v1.90b4.9 64-bit (13 Oct 2017)
Options in effect:
  --allow-extra-chr
  --allow-no-sex
  --bfile /proj/uppstore2017228/KLT.05.GRUS/sophieflatr/newrun/onlyvcf/SNP/all_by_bcf_use_this
  --dog
  --double-id
  --extract fcrpruned_pendleton_r2_0.05_geno_0.05.prune.in
  --fst case-control
  --out fcr_e_pruned_0.5_pendleton
  --pheno /proj/uppstore2017228/KLT.05.GRUS/sophieflatr/newrun/onlyvcf/SNP/phenos/pheno_fcr_vs_eriks.txt
  --remove /proj/uppstore2017228/KLT.05.GRUS/sophieflatr/newrun/onlyvcf/SNP/phenos/3dogs.txt
  --within /proj/uppstore2017228/KLT.05.GRUS/sophieflatr/newrun/onlyvcf/SNP/phenos/pheno_fcr_vs_eriks.txt

Hostname: rackham2.uppmax.uu.se
Working directory: /crex/proj/uppstore2017228/KLT.05.GRUS/sophieflatr/newrun/onlyvcf/SNP
Start time: Fri Jul 24 10:22:10 2020

Random number seed: 1595578930
257386 MB RAM detected; reserving 128693 MB for main workspace.
Allocated 72389 MB successfully, after larger attempt(s) failed.
57459622 variants loaded from .bim file.
309 dogs (0 males, 0 females, 309 ambiguous) loaded from .fam.
Ambiguous sex IDs written to fcr_e_pruned_0.5_pendleton.nosex .
179 phenotype values present after --pheno.
--extract: 1327932 variants remaining.
--remove: 306 dogs remaining.
--within: 2 clusters loaded, covering a total of 179 dogs.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 306 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate in remaining samples is 0.98324.
1327932 variants and 306 dogs pass filters and QC.
Among remaining phenotypes, 19 are cases and 160 are controls.  (127 phenotypes
are missing.)
Writing --fst report (2 populations) to fcr_e_pruned_0.5_pendleton.fst ...
done.
1287278 markers with valid Fst estimates (13162 excluded).
Mean Fst estimate: 0.104456
Weighted Fst estimate: 0.145947

End time: Fri Jul 24 10:22:46 2020
```



