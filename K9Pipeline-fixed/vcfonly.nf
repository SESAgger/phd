 
if (params.help) {
    exit 0, usageMessage()
}

/*
find ud af at sætte anden tid på fastqc relativt til bwa. Spild af booket tid.
--maf 0.000001 is removed from hardfilters_{snp,indels}
*/

checkInputParams()

reference         = file(params.reference)
refdir            = file(reference.getParent())
referenceBwaIndex = file("${reference}.{amb,ann,bwt,pac,sa}")
referenceFaIndex  = file("${reference}.fai")
referenceDict     = file("${refdir}/${reference.getBaseName()}.dict")
known             = file(params.known)

chromosomes = params.chromosomes.split(',').collect { it.trim() }

fastqFiles = Channel.fromFilePairs(params.fastqDir + '/*R{1,2}*.f*q.gz') // /**/ indsat for at medtage subfoldere
fastqFiles.into { fastq_qc; fastq_bwa }

bdir = params.bamDir + '/*.bam'
bamFiles = Channel.fromPath(bdir)
bamFiles
  .map { [it.baseName, it, infer_bam_index_from_bam(it)] }
  .set { readyBamFiles }


vcfFiles = Channel.fromPath(params.vcfDir+'/*vcf.gz') // /**/ indsat for at medtage subfoldere
vcfFiles
  .map { [it.baseName, it, infer_vcf_index_from_vcf(it)] }
  .set { compress_haplocalled }


infoMessage()
	    /* COMBINE GVCF
	    This needs to collect 150/200 samples from samples.list
	    Unsure how to get those
	    Script also runs one per chromosome it seems, could use each
	    */

        //compress_haplocalled.toList().set { collect_haplovcfs }

        compress_haplocalled.toList().transpose().toList().set { collect_haplovcfs }

	    gVCFCombine_ch = Channel.create()
	    genotyping = Channel.create()
	    if ( params.combineByChromosome ) {
	        collect_haplovcfs.set { gVCFCombine_ch }
	        genotyping.close()
	    }
	    else {
	        collect_haplovcfs.set { genotyping }
	        gVCFCombine_ch.close()
	    }


	    process gVCFCombine {
	        tag "$chrom"

	        input:
	            set file(reference), file(refindex), file(refdict) from Channel.value([reference, referenceFaIndex, referenceDict])
	            set val(key), file(vcfs), file(ix_vcfs) from gVCFCombine_ch
	            each chrom from chromosomes
	        output:
	            set val(chrom), file("${chrom}.vcf") into combined

	        when params.combineByChromosome

	        script:
	        """
	        gatk CombineGVCFs \
	            -V ${vcfs.join(' -V ')} \
	            -R $reference \
	            -O  ${chrom}.vcf -L $chrom
	        """
	    }


	    process bgZipCombinedGVCF {
	        tag "$chrom"

	        input:
	            set val(chrom), file(combined_gvcf) from combined
	        output:
	            set val(chrom), file("${combined_gvcf}.gz"), file("*.gz.tbi") into compressed_comb_gvcf


	        """
	        bgzip $combined_gvcf
	        tabix ${combined_gvcf}.gz
	        """
	    }


	    genotyping.mix( compressed_comb_gvcf )
	              .map { it[0] instanceof List ? ['all', it[1], it[2]] : it }
	              .set { genotyping }


	    process genotype {
	        tag "$key"

	        input:
	            set val(key), file(vcfs), file(ix_vcfs) from genotyping
	            set file(reference), file(refindex), file(refdict) from Channel.value([reference, referenceFaIndex, referenceDict])
	        output:
	            set val(key), file("${key}_genotyping.vcf.gz"), file("${key}_genotyping.vcf.gz.tbi") into genotyped


	        script:
	        """
	        gatk GenotypeGVCFs \
	            -R $reference \
	            -V ${vcfs.join(' -V ')} \
	            -O  ${key}_genotyping.vcf.gz

	        """
	    }


	    genotyped.toList().transpose().toList().set { comb_input }


	    process combineChrVCFs {
	        input:
	            set val(keys), file(vcf), file(idx) from comb_input
	        output:
	            set val('all'), file('all.vcf.gz'), file('all.vcf.gz.tbi') into hardfilters

	        publishDir "${params.outdir}/genotype", mode: 'copy'


	        script:
	        """
	        gatk GatherVcfs \
	            -R $reference \
	            -I ${vcf.join(' -I ')} \
	            -O  all.vcf.gz 

	        tabix all.vcf.gz
	        """
	    }



	    hardfilters.into { hardfilters_snp; hardfilters_indel }


	    process hardfilters_snp {
	        input:
	            set val(key), file(vcf), file(index) from hardfilters_snp
	            set file(reference), file(refindex), file(refdict) from Channel.value([reference, referenceFaIndex, referenceDict])
	        output:
	            set file('*SNP*vcf'), file('*filtered_snps*vcf')

	        publishDir "${params.outdir}/genotype", mode: 'copy'


	        script:
	        """
	        gatk SelectVariants -R $reference \
	            -V $vcf -select-type SNP -O  ${key}_raw_SNP_1.vcf
	        gatk VariantFiltration -R $reference -V ${key}_raw_SNP_1.vcf \
	            --filter-expression "QD < 2.0" --filter-name "QD_less_than_2_filter" \
	            --filter-expression "FS > 60.0" --filter-name "FS_greater_than_60_filter" \
	            --filter-expression "SOR > 3.0" --filter-name "SOR_greater_than_3_filter" \
	            --filter-expression "MQ < 40.0" --filter-name "MQ_less_than_40_filter" \
	            --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum_less_than_-12.5_filter" \
	            --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum_less_than_-8_filter" \
	            -O ${key}_filtered_snps_1.vcf

	        vcftools --vcf ${key}_filtered_snps_1.vcf --keep-filtered PASS --out  ${key}_pass_SNP_1 --remove-filtered-geno-all  \
	            --remove-filtered-all --recode --recode-INFO-all --max-maf 0.99992
	         """
	    } 


	    process hardfilters_indel {
	        input:
	            set val(key), file(vcf), file(index) from hardfilters_indel
	            set file(reference), file(refindex), file(refdict) from Channel.value([reference, referenceFaIndex, referenceDict])
	        output:
	            set file('*INDEL*vcf'), file('*filtered_indels*vcf')

	        publishDir "${params.outdir}/genotype", mode: 'copy'


	        script:
	        """
	        gatk SelectVariants -R $reference \
	            -V $vcf -select-type INDEL -O  ${key}_raw_INDEL_1.vcf

	        gatk VariantFiltration -R $reference -V ${key}_raw_INDEL_1.vcf \
	            --filter-expression "QD < 2.0" --filter-name "QD_less_than_2_filter" \
	            --filter-expression "FS > 200.0" --filter-name "FS_greater_than_60_filter" \
	            --filter-expression "SOR > 10.0" --filter-name "SOR_greater_than_3_filter" \
	            --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum_less_than_-8_filter" \
	            -O  ${key}_filtered_indels_1.vcf

	        vcftools --vcf ${key}_filtered_indels_1.vcf  --keep-filtered PASS --out  ${key}_pass_INDEL_1 --remove-filtered-geno-all  \
	            --remove-filtered-all --recode --recode-INFO-all --max-maf 0.99992
	        """
	}


////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS                                                                  //
////////////////////////////////////////////////////////////////////////////////


def checkInputParams() {
    // Check required parameters and display error messages
    boolean fatal_error = false
    if ( ! params.fastqDir && ! params.bamDir && ! params.vcfDir) {
        log.warn("You need to provide a fastqDir (--fastqDir) or a bamDir (--bamDir) or a vcfDir (--vcfDir)")
        fatal_error = true
    }
    if ( ! params.reference ) {
        log.warn("You need to provide a genome reference (--reference)")
        fatal_error = true
    }
    if ( ! params.known ) {
        log.warn("You need to provide a file of known sites (--known)")
        fatal_error = true
    }
    if ( params.onlyMap && !params.fastqDir ) {
        log.warn("You need to specify --fastqDir when doing the mapping (--onlyMap)")
        fatal_error = true
    }
    if ( params.onlyMap && params.bamDir ) {
        log.warn("You need to specify --fastqDir not --bamDir when doing only doing the mapping (--onlyMap)")
        fatal_error = true
    }

    if (fatal_error) {
        log.warn "See --help for more information"
        exit 1
    }
}

def usageMessage() {
    log.info """\
    Usage:
        nextflow run main.nf --fastqDir <directory>
    Options:
        --help
           Print this help message
        --fastqDir <Dir>
           Directory containing fastq samples (.fq.gz)
        --bamDir <Dir>
           Instead of --fastqDir, directory containing bam and indexed bam sample files (.bam, .bam.bai)
        --reference <file>
           Genome reference file (has to be indexed)
        --known <file>
           File with known sites for quality calibration.
        --O dir <dir>
           Directory for output files
        --onlyMap
           Only run the mapping steps
        --project
           Slurm project to run with
        --QCdone 
        	Without this flag, only fastQC ± Adapterremoval is done
        --AR
        	Will remove adapters as specified in file
    """
}

def infoMessage() {
    log.info """\
*** K9 WGS Pipeline ***
Configuration environment:
    Out directory:   $params.outdir
    Fastq directory: $params.fastqDir
    Bam directory:   $params.bamDir
    Reference:       $params.reference
    OnlyMap?:        $params.onlyMap
    """
}

def infer_bam_index_from_bam(f) {
    // If the ".bam.bai" file does not exist, try ".bai" without ".bam"
    return infer_filepath(f, /$/, '.bai')
        ?: infer_filepath(f, /.bam$/, '.bai')
        ?: filepath_from(f, /$/, '.bai') // Default filename if none exist
}

def infer_vcf_index_from_vcf(f) {
    // If the ".vcf.tbi" file does not exist, try ".tbi" without "vcf.tbi"
    return infer_filepath(f, /$/, '.tbi')
        ?: infer_filepath(f, /.tbi$/, '.tbi')
        ?: filepath_from(f, /$/, '.tbi') // Default filename if none exist
}

def filepath_from(from, match, replace) {
    path = file( from.toString().replaceAll(match, replace) )
    return path
}

def infer_filepath(from, match, replace) {
    path = file( from.toString().replaceAll(match, replace) )
    if (path.exists()) {
        return path
    }
    return false
}
