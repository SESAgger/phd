[toc]

# Pipeline

Selve pipeline run dvs. fra FastQC til GT

```bash
# T.o.m haplotypeCaller, chrashede derefter, så blev startet op igen med vcfonly
export NXF_LAUNCHER=$SNIC_TMP    # For Rackham
export NXF_TEMP=$SNIC_TMP        # For Rachham
export NXF_SINGULARITY_CACHEDIR=/proj/uppstore2017228/KLT.05.GRUS/sophieflatr/K9-WGS-Pipeline/singularity

nextflow -log runallagain.log run -resume /proj/uppstore2017228/KLT.05.GRUS/sophieflatr/K9-WGS-Pipeline/main.nf  -profile rackham \
        --fastqDir /home/sesagger/FCRfastqs/ \
        --reference /proj/uppstore2017228/KLT.05.GRUS/DATA_merge/canFam3.1.fa \
        --known /proj/uppstore2017228/KLT.05.GRUS/DATA_merge/dog_dbSNP.recode.vcf \
        --outdir /proj/uppstore2017228/KLT.05.GRUS/sophieflatr/newrun \
        --project snic2017-7-384
        
# Run som stod for genotypning og derfra        
export NXF_LAUNCHER=$SNIC_TMP    # For Rackham
export NXF_TEMP=$SNIC_TMP        # For Rachham
export NXF_SINGULARITY_CACHEDIR=/proj/uppstore2017228/KLT.05.GRUS/sophieflatr/K9-WGS-Pipeline/singularity

nextflow -log runallagain.log run /proj/uppstore2017228/KLT.05.GRUS/sophieflatr/K9-WGS-Pipeline/vcf-new.nf  -profile rackham1 \
	--vcfDir /proj/uppstore2017228/KLT.05.GRUS/sophieflatr/newrun/haplotypeCaller \
        --reference /proj/uppstore2017228/KLT.05.GRUS/DATA_merge/canFam3.1.fa \
        --known /proj/uppstore2017228/KLT.05.GRUS/DATA_merge/dog_dbSNP.recode.vcf \
        --outdir /proj/uppstore2017228/KLT.05.GRUS/sophieflatr/newrun/onlyvcfs \
        --project snic2017-7-384
```

## Make files

Jge kaf ifjk stave

### Annotate files

```bash
eriks
---------------------------------
module load bioinfo-tools htslib/1.10 bcftools/1.10
bcftools annotate 160DG.99.9.recalibrated.variants.eff.vcf.gz --set-id +'%CHROM\_%POS' --output erik_renamed.vcf

bgzip -c erik_renamed.vcf > erik_renamed.vcf.gz

tabix erik_renamed.vcf.gz

iDog
---------------------------------
bcftools annotate Filtred_Published.vcf.gz --set-id +'%CHROM\_%POS' --output iDog-filtered_renamed.vcf

bgzip -c iDog-filtered_renamed.vcf > iDog-filtered_renamed.vcf.gz

tabix iDog-filtered_renamed.vcf.gz
```

### Merging

```bash
module load bioinfo-tools bcftools/1.10 vcftools/0.1.15 htslib/1.10 plink/1.90b4.9 plink2/2.00-alpha-2-20190429

bcftools merge mine_hunde_renamed.vcf.gz ../../../eriks/erik_renamed.vcf.gz ../iDog/iDog-filtered_renamed.vcf.gz --output all_by_bcfttools.vcf

bgzip -c all_by_bcfttools.vcf > all_by_bcfttools.vcf.gz

tabix all_by_bcfttools.vcf.gz
```

## Analyses

### Calculation of pooled heterogozity

```bash
module load bioinfo-tools bcftools/1.10 vcftools/0.1.15 htslib/1.10 plink/1.90b4.9 plink2/2.00-alpha-2-20190429
plink \
  --allow-extra-chr \
  --allow-no-sex  \
  --bfile all_by_bcf_use_this \
  --dog \
  --double-id \
  --freq counts \
  --keep phenos/pheno_hsfcr_vs_fcr.txt  \
  --out fcr_frq_count_all \
  --pheno phenos/pheno_hsfcr_vs_fcr.txt
  
grep -v "NaN" fcr_frq_count_all.frq.counts
```

```R
#Libraries
# Rcpp_1.0.3           compiler_3.6.0       pillar_1.4.3  
# prettyunits_1.1.1    progress_1.2.2       bitops_1.0-6  
# tools_3.6.0          digest_0.6.23        zeallot_0.1.0 
# bit_1.1-15.2         RSQLite_2.2.0        memoise_1.1.0 
# lifecycle_0.1.0      tibble_2.1.3         gtable_0.3.0  
# pkgconfig_2.0.3      rlang_0.4.0          DBI_1.1.0     
# curl_4.3             parallel_3.6.0       httr_1.4.1    
# hms_0.5.1            IRanges_2.18.3       S4Vectors_0.22
# vctrs_0.2.0          gtools_3.8.1         stats4_3.6.0  
# bit64_0.9-7          grid_3.6.0           tidyselect_0.2
# Biobase_2.44.0       glue_1.3.1           calibrate_1.7.
# R6_2.4.1             AnnotationDbi_1.46.1 XML_3.99-0.3  
# purrr_0.3.2          blob_1.2.0           magrittr_1.5  
# BiocGenerics_0.30.0  backports_1.1.5      MASS_7.3-51.4 
# assertthat_0.2.1     colorspace_1.4-1     stringi_1.4.5 
# RCurl_1.98-1.1       crayon_1.3.4         tidyr_1.0.2 
# biomaRt_2.40.5 			 forcats_0.4.0  			gdata_2.18.0      
# stringr_1.4.0        dplyr_0.8.4    			qqman_0.1.4   

hentfil("fcr_freq_counts_no_allmiss.tsv",df)
df$Location<-str_split_fixed(df$SNP,"_",2)[,2]
hp_calc<-function (df,bin_size,bin_step)
{ 
    df3<-c()
    i<-0
    j<-1
    while (j<39){
        a<-subset(df,CHR==j)
        if(dim(a)[1]!=0){
            i<-0
            while (i < max(a$Location)){
                bin_start<-i
                bin_end<-i+bin_size
                b<-subset(a,a$Location>=bin_start&a$Location<bin_end)
                C1_sum<-sum(b$C1,na.rm=TRUE)
                C2_sum<-sum(b$C2,na.rm=TRUE)
                G0_sum<-sum(b$G0,na.rm=TRUE)
                df2<-cbind(j,bin_start,bin_end,C1_sum,C2_sum,G0_sum)
                i<-i+bin_step
                if(C1_sum!=0|C2_sum!=0|G0_sum!=0){
                    df3<-bind_rows(df3,as.data.frame(df2))
                }
            }
        }     
        j<-j+1
    }
    df3$hp<-2*df3$C1_sum*df3$C2_sum/(df3$C1_sum+df3$C2_sum)^2
    df3$Z_hp<-scale(df3$hp,center=TRUE,scale=TRUE)
    assign(paste(deparse(substitute(df)), "summed", sep = "_"), 
        df3, envir = .GlobalEnv)
}
```



### Calculate $\mathrm{F_{ST}}$ 

#### Windowed

##### Make the file

```bash
module load bioinfo-tools bcftools/1.10 vcftools/0.1.15 htslib/1.10 plink/1.90b4.9 plink2/2.00-alpha-2-20190429

vcftools \
    --gzvcf all_by_bcfttools.vcf.gz \
    --fst-window-size {100000/200000} \
    --fst-window-step {50000/100000} \
    --weir-fst-pop fcr_for_fst.txt \
    --weir-fst-pop {all,iDog,eriks}_for_fst.txt \
    --keep fcr_for_fst.txt \
    --keep {all,iDog,eriks}_for_fst.txt \
    --out fcr_{all,iDog,eriks}_fst_w_{100000,200000}_s_{50000,100000}kb
```

##### Databehandling

$Weighted-F_{ST}$ is used, unless otherwise specified.

Fejltjek:

```R 
str(df)
summary(df)
#new_df <- something(df)
dim(newdf)==dim(df)
# Summarize
  df %>% dplyr::select(-Feature) %>% subset(filter options) %>% distinct(Location,.keep_all=TRUE) #Check length=length(count(Location)) - to see if any locations have been dropped
MAF=ifelse(as.character(fcr_e_g5$Allele)==as.character(fcr_e_g5$A1),fcr_e_g5$MAF_A,1-fcr_e_g5$MAF_A)#To take into accountt thatt plink somettimes switches the min og maj
a<-c("CHR","SNP","Gene","Location","POS","NMISS","FST","Z_FST","#Uploaded_variation","Allele","Consequence","Amino_acids","Codons","Extra","external_gene_name","A1","MAF_A","MAF_U","r_MAF_A")#Relevante kolonner
```

```R
#Make tables:
df %>% dplyr::select(-Feature) %>% subset(filter options) %>% distinct(Location,.keep_all=TRUE)%>% dplyr::select(Consequence ,external_gene_name,Gene,Location,FST , Z_FST,MAF_A,MAF_U,Codons,Extra,NCHROBS_A,NCHROBS_U)

dplyr::select(Consequence ,external_gene_name,Gene,Location,A1,A2,NCHROBS_A,FST , Z_FST,MAF_A,Codons,Extra)
```



##### Figures

```R
gg.manhattan_chr
function (df1, threshold, hlight, col, ylims, title) 
{
    if (is.data.frame(df1) == FALSE) {
        print("Input must be a data.frame")
    }
    else if (is.data.table(df1) == TRUE) {
        print("Input must be data.frame, not data.table")
    }
    else {
        if ((missing(threshold))) {
            threshold <- 5
        }
        if ((missing(col))) {
            col <- c("#5D82BB", "#3B64A5", "#1E4F9E", "#103B7E", 
                "#082B64")
        }
        if ((missing(title))) {
            title <- "Unknown"
        }
        df <- subset(df1, Z_FST > 0)
        df.tmp <- df %>% group_by(CHR) %>% summarise(chr_len = max(POS)) %>% 
            mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>% 
            dplyr::select(-chr_len) %>% left_join(df, ., by = c(CHR = "CHR")) %>% 
            arrange(CHR, POS) %>% mutate(BPcum = POS + tot) %>% 
            mutate(is_annotate = ifelse(Z_FST > threshold, "yes", 
                "no")) %>% mutate(is_annotate1 = ifelse(Z_FST > 
            5.5, "yes", "no"))
        axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center = (max(BPcum) + 
            min(BPcum))/2)
        ggplot(df.tmp, aes(x = BPcum, y = Z_FST)) + geom_point(aes(color = as.factor(CHR)), 
            alpha = 0.8, size = 2) + scale_color_manual(values = rep(col, 
            38)) + scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) + 
            scale_y_continuous(expand = c(0, 0), limits = ylims) + 
            ggtitle(paste0(title)) + labs(x = "Chromosome", y = expression(ZF["ST"])) + 
            geom_hline(yintercept = 5, color = "orange", linetype = "dashed") + 
            geom_hline(yintercept = 5.5, color = "orange") + 
            geom_point(data = subset(df.tmp, is_annotate == "yes"), 
                color = "orange", size = 2) + theme_bw(base_size = 22) + 
            theme(plot.title = element_text(hjust = 0.5), legend.position = "none", 
                panel.border = element_blank(), panel.grid.major.x = element_blank(), 
                panel.grid.minor.x = element_blank())
    }
}
```



#### All

```bash
module load bioinfo-tools bcftools/1.10 vcftools/0.1.15 htslib/1.10 plink/1.90b4.9 plink2/2.00-alpha-2-20190429
plink \
  --allow-extra-chr \
  --allow-no-sex \
  --bfile all_by_bcf_use_this \
  --dog \
  --double-id \
  --fst case-control \
  --geno 0.05 \
  --keep phenos/pheno_fcr_vs_{both,iDog,eriks}.txt \
  --out fcr_vs_{both,iDog,eriks}_fst_g5percent.fst \
  --pheno phenos/pheno_fcr_vs_{both,iDog,eriks}.txt \
  --within phenos/pheno_fcr_vs_{both,iDog,eriks}.txt
```

### Annotering af SNPs

```bash

```

## Ostranders

Dogs that in file that could be found under bioproject, was included