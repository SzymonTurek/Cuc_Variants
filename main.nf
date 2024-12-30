#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastq_files  = file("$baseDir/samples.csv")
params.reference_genome = file("$baseDir/reference_genome/pb_b10_ill1.fasta")
params.reference_genome_dir = "$baseDir/reference_genome"
params.fastqc_outdir="$baseDir/fastqc_output"
params.BWA_outdir="$baseDir/BWA_output"
params.Bowtie2_outdir="$baseDir/Bowtie2_output"
params.BCFTOOLS_variants_outdir="$baseDir/BCFTOOLS_variants_outdir"
params.BCFTOOLS_BOWTIE2_variants_outdir="$baseDir/BCFTOOLS_variants/BCFTOOLS_BOWTIE2_variants_outdir"
params.BCFTOOLS_BWA_variants_outdir="$baseDir/BCFTOOLS_variants/BCFTOOLS_BWA_variants_outdir"
params.FREEBAYES_BOWTIE2_variants_outdir="$baseDir/FREEBAYES_variants/FREEBAYES_BOWTIE2_variants_outdir"
params.FREEBAYES_BWA_variants_outdir="$baseDir/FREEBAYES_variants/FREEBAYES_BWA_variants_outdir"
params.FREEBAYES_variants_outdir="$baseDir/FREEBAYES_variants_outdir"

process FastQC {

    tag "FastQC on $sample_ID"
    publishDir params.fastqc_outdir, mode: 'copy'

    input:
    tuple val(sampleId), path(read1), path(read2) from sample_ch

    output:
    path "*.{html,zip}" into QC_Report

    script:
    """
    fastqc -t 2 -q "${read1}" "${read2}"
    """
}


process READ_FASTQS {
    input:
    tuple val(sampleId), path(read1), path(read2)

    script:
    """
    echo Sample_ID: $sampleId Read1: $read1 Read2: $read2
    """
}


process FASTQC {
    tag "FastQC on $sample_id"
    publishDir params.fastqc_outdir, mode: 'copy'


    input:
    tuple val(sample_id), path(read1), path(read2)
        

    output:
    path "fastqc_${sample_id}"

    script:
    """
    mkdir fastqc_${sample_id}
    fastqc -o fastqc_${sample_id} -f fastq -q ${read1} -q ${read2} 
    """


}


process MULTIQC {

    publishDir params.fastqc_outdir, mode: 'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'


    script:
    """
    multiqc .
    """

}

process BWA_INDEX {
   
    publishDir params.reference_genome_dir, mode: 'copy'
    
    input:

    path reference_genome

    output:
    path '*'

    script:
    """
    bwa index ${reference_genome}
    """
}


process BWA_ALIGN {
    publishDir params.BWA_outdir, mode: 'copy'

    label 'big_mem'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path reference_genome
    path reference_genome_dir

    output:
    path "${sample_id}_bwa.bam" 
    path "${sample_id}_bwa_sorted.bam"  , emit: sorted_bam 
    path "${sample_id}_mean_read_depth.log"
    path "${sample_id}_breadth_of_coverage.log" 
    path "${sample_id}_flagstats.log"
    val  "${sample_id}" , emit: sample_id 

    script:
    """
    bwa mem -t 12 ${reference_genome} ${read1} ${read2} | samtools view -bS  > ${sample_id}_bwa.bam
    samtools sort ${sample_id}_bwa.bam -o ${sample_id}_bwa_sorted.bam
    samtools depth -a ${sample_id}_bwa_sorted.bam | awk '{c++;s+=\$3}END{print s/c}' > ${sample_id}_mean_read_depth.log
    samtools depth -a ${sample_id}_bwa_sorted.bam | awk '{c++; if(\$3>0) total+=1}END{print (total/c)*100}' > ${sample_id}_breadth_of_coverage.log
    samtools flagstat ${sample_id}_bwa_sorted.bam > ${sample_id}_flagstats.log

    """
}


process BOWTIE_INDEX {
   
    publishDir params.reference_genome_dir, mode: 'copy'
    
    input:

    path reference_genome

    output:
    path '*'

    script:
    """
    bowtie2-build ${reference_genome} bowtie_index    
    """
}



process BOWTIE2_ALIGN {
    
    publishDir params.Bowtie2_outdir, mode: 'copy'
    
    label 'big_mem'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    path(reference_genome_dir)

    output:
    path "${sample_id}_bowtie.bam" 
    path "${sample_id}_bowtie_sorted.bam" , emit: sorted_bam 
    path "${sample_id}_mean_read_depth.log" 
    path "${sample_id}_breadth_of_coverage.log" 
    path "${sample_id}_flagstats.log" 
    val  "${sample_id}" , emit: sample_id 

    script:
    """
    bowtie2 -x bowtie_index -p 12 -1 ${read1} -2 ${read2}  | samtools view -bS  > ${sample_id}_bowtie.bam
    samtools sort ${sample_id}_bowtie.bam -o ${sample_id}_bowtie_sorted.bam
    samtools depth -a ${sample_id}_bowtie_sorted.bam | awk '{c++;s+=\$3}END{print s/c}' > ${sample_id}_mean_read_depth.log
    samtools depth -a ${sample_id}_bowtie_sorted.bam | awk '{c++; if(\$3>0) total+=1}END{print (total/c)*100}' > ${sample_id}_breadth_of_coverage.log
    samtools flagstat ${sample_id}_bowtie_sorted.bam > ${sample_id}_flagstats.log

    """
}

process BCFTOOLS_CALL_VARIANTS_BOWTIE2 {
    publishDir params.BCFTOOLS_BOWTIE2_variants_outdir, mode: 'copy'
    
    input:
    path(sorted_bam)
    path(ref)
    val(sample_id)


    output:
    path("${sample_id}_bcftools_variants.vcf"), emit: bcftols_vcf 
    val("${sample_id}") , emit: sample_id
    
    script:
    
    """
    bcftools mpileup -Ou -f ${ref} ${sorted_bam} | bcftools call -mv -Ov -o ${sample_id}_bcftools_variants.vcf

    """
}


process BCFTOOLS_CALL_VARIANTS_BWA {
    publishDir params.BCFTOOLS_BWA_variants_outdir, mode: 'copy'
    
    input:
    path(sorted_bam)
    path(ref)
    val(sample_id)


    output:
    path("${sample_id}_bcftools_variants.vcf"), emit: bcftols_vcf 
    val("${sample_id}") , emit: sample_id
    
    script:
    
    """
    bcftools mpileup -Ou -f ${ref} ${sorted_bam} | bcftools call -mv -Ov -o ${sample_id}_bcftools_variants.vcf

    """
}



process FREEBAYES_CALL_VARIANTS_BWA {
    publishDir params.FREEBAYES_BWA_variants_outdir, mode: 'copy'
    
    input:
    path(sorted_bam)
    path(ref)
    val(sample_id)


    output:
    path("${sample_id}_freebayes_variants.vcf") , emit: freebayes_vcf 
    val("${sample_id}") , emit: sample_id
    
    script:
    
    """
    freebayes -f ${ref} ${sorted_bam} > ${sample_id}_freebayes_variants.vcf
    """
}


process FREEBAYES_CALL_VARIANTS_BOWTIE2 {
    publishDir params.FREEBAYES_BOWTIE2_variants_outdir, mode: 'copy'
    
    input:
    path(sorted_bam)
    path(ref)
    val(sample_id)


    output:
    path("${sample_id}_freebayes_variants.vcf") , emit: freebayes_vcf 
    val("${sample_id}") , emit: sample_id
    
    script:
    
    """
    freebayes -f ${ref} ${sorted_bam} > ${sample_id}_freebayes_variants.vcf
    """
}


process MARK_DUPLICATES_BOWTIE2 {
    publishDir params.Bowtie2_outdir, mode: 'copy'
    
    input:
    path(sorted_bam)
    val(sample_id)


    output:
    path("${sample_id}_marked_duplicates.bam")
    path("${sample_id}_marked_dup_metrics.txt")
    path("${sample_id}_marked_dup_flagstats.log")

    script:
    
    """
    gatk AddOrReplaceReadGroups -I ${sorted_bam} -O ${sample_id}_marked_duplicates_head.bam --RGLB library --RGPL illumina --RGPU barcode --RGSM name
    gatk MarkDuplicatesSpark -I ${sample_id}_marked_duplicates_head.bam -O ${sample_id}_marked_duplicates.bam -M ${sample_id}_marked_dup_metrics.txt
    samtools flagstat ${sample_id}_marked_duplicates.bam > ${sample_id}_marked_dup_flagstats.log
    """
}


process MARK_DUPLICATES_BWA {
    publishDir params.BWA_outdir, mode: 'copy'
    
    input:
    path(sorted_bam)
    val(sample_id)


    output:
    path("${sample_id}_marked_duplicates.bam")
    path("${sample_id}_marked_dup_metrics.txt")
    path("${sample_id}_marked_dup_flagstats.log")

    script:
    
    """
    gatk AddOrReplaceReadGroups -I ${sorted_bam} -O ${sample_id}_marked_duplicates_head.bam --RGLB library --RGPL illumina --RGPU barcode --RGSM name
    gatk MarkDuplicatesSpark -I ${sample_id}_marked_duplicates_head.bam -O ${sample_id}_marked_duplicates.bam -M ${sample_id}_marked_dup_metrics.txt
    samtools flagstat ${sample_id}_marked_duplicates.bam > ${sample_id}_marked_dup_flagstats.log
    """
}


process FILTER_VARIANTS_BWA_FREEBAYES {

    publishDir params.FREEBAYES_BWA_variants_outdir, mode: 'copy'
    
    input:
    path(vcf_file)
    path(reference_genome)
    val(sample_id)
    path(fai_idx)
    path(fasta_dict)


    output:
    path("${sample_id}_filtered.vcf")

    script:
    
    """
    gatk VariantFiltration -R ${reference_genome} -V ${vcf_file} -O ${sample_id}_filtered.vcf --filter-name "my_filter1"  --filter-expression "MQ<20" --filter-name "my_filter2" --filter-expression "DP<5"

    """
}

process FILTER_VARIANTS_BOWTIE_FREEBAYES {

    publishDir params.FREEBAYES_BOWTIE2_variants_outdir, mode: 'copy'
    
    input:
    path(vcf_file)
    path(reference_genome)
    val(sample_id)
    path(fai_idx)
    path(fasta_dict)


    output:
    path("${sample_id}_filtered.vcf")

    script:
    
    """
    gatk VariantFiltration -R ${reference_genome} -V ${vcf_file} -O ${sample_id}_filtered.vcf --filter-name "my_filter1"  --filter-expression "MQ<20" --filter-name "my_filter2" --filter-expression "DP<5"

    """
}

process FILTER_VARIANTS_BWA_BCFTOOLS {

    publishDir params.BCFTOOLS_BWA_variants_outdir, mode: 'copy'
    
    input:
    path(vcf_file)
    path(reference_genome)
    val(sample_id)
    path(fai_idx)
    path(fasta_dict)


    output:
    path("${sample_id}_filtered.vcf")

    script:
    
    """
    gatk VariantFiltration -R ${reference_genome} -V ${vcf_file} -O ${sample_id}_filtered.vcf --filter-name "my_filter1"  --filter-expression "MQ<20" --filter-name "my_filter2" --filter-expression "DP<5"

    """
}

process FILTER_VARIANTS_BOWTIE_BCFTOOLS {

    publishDir params.BCFTOOLS_BOWTIE2_variants_outdir , mode: 'copy'
    
    input:
    path(vcf_file)
    path(reference_genome)
    val(sample_id)
    path(fai_idx)
    path(fasta_dict)


    output:
    path("${sample_id}_filtered.vcf")

    script:
    
    """
    gatk VariantFiltration -R ${reference_genome} -V ${vcf_file} -O ${sample_id}_filtered.vcf --filter-name "my_filter1"  --filter-expression "MQ<20" --filter-name "my_filter2" --filter-expression "DP<5"

    """
}


process INDEX_FASTA {
    publishDir params.reference_genome_dir, mode: 'copy'
    
    input:

    path reference_genome

    output:
    path '*.fai'

    script:
    """
    samtools faidx ${reference_genome}
    """
}
   

process CREATE_SEQUENCE_DICTIONARY {
    publishDir params.reference_genome_dir, mode: 'copy'
    
    input:

    path reference_genome

    output:
    path '*.dict'

    script:
    """
    gatk CreateSequenceDictionary -R ${reference_genome}
    """
}



workflow {

    samples_ch = Channel.fromPath( params.fastq_files )
    .splitCsv( header: true, sep: ',' )
    .map { row -> tuple( row.sample_id, file(row.fastq1), file(row.fastq2) ) }



    //samples_ch = Channel.fromPath(params.fastq_files)
    //    .splitCsv(header: true)
    //    .map { row -> tuple(row.sample_id, file(row.read1)) }

    
    //collected_samples_ch = Channel.fromPath(params.fastq_files)
    //    .splitCsv(header: true)
    //    .map { row -> tuple(file(row.read1)) }
    //    .collect()
    //    .view()


    READ_FASTQS(samples_ch)
    fastqc_ch = FASTQC(samples_ch)

    MULTIQC(fastqc_ch.collect())
    fasta_index_ch = INDEX_FASTA(params.reference_genome)
    fasta_dict_ch = CREATE_SEQUENCE_DICTIONARY(params.reference_genome)
    
    BWA_INDEX_ch = BWA_INDEX(params.reference_genome)
    BWA_ALIGN(samples_ch, params.reference_genome, BWA_INDEX.out)
    
    BOWTIE_INDEX(params.reference_genome)
    BOWTIE2_ALIGN(samples_ch, BOWTIE_INDEX.out)
    
    MARK_DUPLICATES_BOWTIE2(BOWTIE2_ALIGN.out.sorted_bam, BOWTIE2_ALIGN.out.sample_id)
    MARK_DUPLICATES_BWA(BWA_ALIGN.out.sorted_bam, BWA_ALIGN.out.sample_id)
    
    BCFTOOLS_CALL_VARIANTS_BOWTIE2(BOWTIE2_ALIGN.out.sorted_bam, params.reference_genome, BOWTIE2_ALIGN.out.sample_id)
    BCFTOOLS_CALL_VARIANTS_BWA(BWA_ALIGN.out.sorted_bam, params.reference_genome, BWA_ALIGN.out.sample_id)
    
    FREEBAYES_CALL_VARIANTS_BWA(BOWTIE2_ALIGN.out.sorted_bam, params.reference_genome, BOWTIE2_ALIGN.out.sample_id)
    FREEBAYES_CALL_VARIANTS_BOWTIE2(BWA_ALIGN.out.sorted_bam, params.reference_genome, BWA_ALIGN.out.sample_id)
    
    FILTER_VARIANTS_BWA_FREEBAYES( FREEBAYES_CALL_VARIANTS_BWA.out.freebayes_vcf, params.reference_genome, FREEBAYES_CALL_VARIANTS_BWA.out.sample_id, fasta_index_ch, fasta_dict_ch )
    FILTER_VARIANTS_BOWTIE_FREEBAYES(FREEBAYES_CALL_VARIANTS_BOWTIE2.out.freebayes_vcf, params.reference_genome, FREEBAYES_CALL_VARIANTS_BOWTIE2.out.sample_id, fasta_index_ch, fasta_dict_ch )
    FILTER_VARIANTS_BWA_BCFTOOLS(BCFTOOLS_CALL_VARIANTS_BWA.out.bcftols_vcf, params.reference_genome, BCFTOOLS_CALL_VARIANTS_BWA.out.sample_id, fasta_index_ch, fasta_dict_ch )
    FILTER_VARIANTS_BOWTIE_BCFTOOLS(BCFTOOLS_CALL_VARIANTS_BOWTIE2.out.bcftols_vcf, params.reference_genome, BCFTOOLS_CALL_VARIANTS_BOWTIE2.out.sample_id, fasta_index_ch, fasta_dict_ch )

}
