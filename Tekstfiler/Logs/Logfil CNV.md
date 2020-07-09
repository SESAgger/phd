Logfil CNV

```
gatk PreprocessIntervals \
	-R /proj/uppstore2017228/KLT.05.GRUS/DATA_merge/canFam3.1.fa  \
	--padding 0 \
	-imr OVERLAPPING_ONLY \
	-O canfam3.1.preprocessed.interval_list
gatk AnnotateIntervals \
	-L canfam3.1.preprocessed.interval_list \
	-R /proj/uppstore2017228/KLT.05.GRUS/DATA_merge/canFam3.1.fa \
	-imr OVERLAPPING_ONLY \
	-O canfam3.1.annotated.
```

