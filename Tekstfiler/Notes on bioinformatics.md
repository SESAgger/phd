# BCFTOOLS

```bash
bcftools annotate fil.vcf.gz --set-id +'%CHROM\_%POS' --output fil_renamed.vcf # Add snp-name to vcf
```







# VEP

```bash
vep --cache $VEP_CAHCE \
--species canis_familiaris -i fil.vcf.gz 	\
--flag-pick  \  #Flags the most likely consequence
--sift b \ #Annotates with SIFT both value and text
-o annotation_w_sift_pick 
```

