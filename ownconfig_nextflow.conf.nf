profiles {
  own {
    params {
      config_profile_description = 'Testing & development profile for UPPMAX, provided by nf-core/configs.'
      // Max resources to be requested by a devel job
      max_memory = 120.GB
      max_time = 5.h
      markdup_java_options="-Xmx6g -Xms4g"
      singleCPUMem = null
      totalMemory = 120.GB
      max_cpus = 1
    }
    process.queue = 'core'
    process {
      withName:CNVkit{
        time = 
      }
      withName:Output_documentation{
       time = 10.min
      }
      withName:BamQC {
        time =5.h
      }
      withName:MapReads{
        time = 16.h
      }
      withName:BuildIntervals{
        time = 20.m
      }
      withName:CreateIntervalBeds{
        time = 10.m
      }
      withName:MergeBamMapped{
        time = 2.h
      }
      withName:MarkDuplicates{
        time = 10.h
        memory = null
        cpus = 10

      }
      withName:BaseRecalibrator{
        time = 1.h
      }
      withName:ApplyBQSR{
        time = 30.m
      }
      withName:get_software_versions{
        time = 5.m
        cpus = 1
      }
    }
  }
}


//Max BamQC time= 3 hours
//Max MapReads = 14 hours
//Max BuildIntervals = 11 min
//Max CreateIntervalBeds = 5 min
//Max MergeBam_Mapped = 41 min
//Max MarkDuplicates = 5 hours
//Max BaseRecalibrator = 25 min


nextflow run /proj/uppstore2017228/KLT.05.GRUS/sophie_sarekstuff/sarek/main.nf --project snic2017-7-384 -profile uppmax,own --genome CanFam4 --step mapping --input /home/sesagger/KLT5/sophie_sarekstuff/test --intervals /proj/uppstore2017228/KLT.05.GRUS/TI-2654/run/results/reference_genome/cf4.b6.14.fa.bed -c ownconfig.conf


//withName:FreebayesSingle, {time = x.h*task.attempt)}}
//withName:VEP{time = x.h*task.attempt)}}
//withName:Mpileup{time = x.h*task.attempt)}}
//withName:MergeMutect2Stats{time = x.h*task.attempt)}}
//withName:UMIMapBamFile{time = x.h*task.attempt)}}
//withName:IndexBamRecal{time = x.h*task.attempt)}}
//withName:Ascat{time = x.h*task.attempt)}}
//withName:AlleleCounter{time = x.h*task.attempt)}}
//withName:Vcftools{time = x.h*task.attempt)}}
//withName:Sentieon_BQSR{time = x.h*task.attempt)}}
//withName:BuildGermlineResourceIndex{time = x.h*task.attempt)}}
//withName:TrimGalore{time = x.h*task.attempt)}}
//withName:Snpeff{time = x.h*task.attempt)}}
//withName:BuildKnownIndelsIndex{time = x.h*task.attempt)}}
//withName:GroupReadsByUmi{time = x.h*task.attempt)}}
//withName:StrelkaSingle{time = x.h*task.attempt)}}
//withName:Sentieon_TNscope{time = x.h*task.attempt)}}
//withName:HaplotypeCaller{time = x.h*task.attempt)}}
//withName:Sentieon_DNAseq{time = x.h*task.attempt)}}
//withName:StrelkaBP{time = x.h*task.attempt)}}
//withName:Output_documentation{time = x.h*task.attempt)}}
//withName:ControlFREEC{time = x.h*task.attempt)}}
//withName:PileupSummariesForMutect2{time = x.h*task.attempt)}}
//withName:ControlFREECSingle{time = x.h*task.attempt)}}
//withName:MantaSingle{time = x.h*task.attempt)}}
//withName:BuildPonIndex{time = x.h*task.attempt)}}
//withName:ControlFreecViz{time = x.h*task.attempt)}}
//withName:CompressVCFvep{time = x.h*task.attempt)}}
//withName:BuildFastaFai{time = x.h*task.attempt)}}
//withName:CalculateContamination{time = x.h*task.attempt)}}
//withName:VEPmerge{time = x.h*task.attempt)}}
//withName:IndexBamFile{time = x.h*task.attempt)}}
//withName:BcftoolsStats{time = x.h*task.attempt)}}
//withName:Strelka{time = x.h*task.attempt)}}
//withName:MergePileupSummaries{time = x.h*task.attempt)}}
//withName:UMIFastqToBAM{time = x.h*task.attempt)}}
//withName:IndexBamMergedForSentieon{time = x.h*task.attempt)}}
//withName:BuildDict{time = x.h*task.attempt)}}
//withName:FastQCFQ{time = x.h*task.attempt)}}
//withName:Manta{time = x.h*task.attempt)}}
//withName:GenotypeGVCFs{time = x.h*task.attempt)}}
//withName:CompressVCFsnpEff{time = x.h*task.attempt)}}
//withName:MergeMpileup{time = x.h*task.attempt)}}
//withName:GatherBQSRReports{time = x.h*task.attempt)}}
//withName:BuildDbsnpIndex{time = x.h*task.attempt)}}
//withName:TIDDIT{time = x.h*task.attempt)}}
//withName:FastQCBAM{time = x.h*task.attempt)}}
//withName:MSIsensor_msi{time = x.h*task.attempt)}}
//withName:CNVkit{time = x.h*task.attempt)}}
//withName:MultiQC{time = x.h*task.attempt)}}
//withName:MergeBamRecal{time = x.h*task.attempt)}}
//withName:SamtoolsStats{time = x.h*task.attempt)}}
//withName:BuildBWAindexes{time = x.h*task.attempt)}}
//withName:FreeBayes{time = x.h*task.attempt)}}
//withName:CompressSentieonVCF{time = x.h*task.attempt)}}
//withName:ConvertAlleleCounts{time = x.h*task.attempt)}}
//withName:ControlFreecVizSingle{time = x.h*task.attempt)}}
//withName:CallMolecularConsensusReads{time = x.h*task.attempt)}}
//withName:Sentieon_MapReads{time = x.h*task.attempt)}}
//withName:Mutect2{time = x.h*task.attempt)}}
//withName:ConcatVCF{time = x.h*task.attempt)}}
//withName:Sentieon_DNAscope{time = x.h*task.attempt)}}
//withName:FilterMutect2Calls{time = x.h*task.attempt)}}
//withName:MSIsensor_scan{time = x.h*task.attempt)}}
//withName:Sentieon_Dedup{time = x.h*task.attempt)}}
//withName:ConcatVCF_Mutec{time = x.h*task.attempt)}}