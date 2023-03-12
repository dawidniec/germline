/*
 * CHANNELS
 */

ch_reference                    = Channel.fromPath( params.reference, checkIfExists: true )
ch_capture_targets_bed          = Channel.fromPath( params.capture_targets_bed, checkIfExists: true )
ch_primary_targets_bed          = Channel.fromPath( params.primary_tatrgets_bed, checkIfExists: true )
ch_bait_interval_list           = Channel.fromPath( params.bait_interval_list, checkIfExists: true )
ch_gold_standard_vcf            = Channel.fromPath( params.gold_standard_vcf, checkIfExists: true )
ch_high_confidence_regions_bed  = Channel.fromPath( params.high_confidence_regions_bed, checkIfExists: true )
ch_input_fastq                  = Channel.fromPath( params.input, checkIfExists: true )
    .map{ file -> 
            def key = file.simpleName.toString().tokenize("_")[0]
            return tuple(key, file) 
        }
    .groupTuple()

/*
 * PROCESSES
 */
 
process CAT_FASTQ {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.fastq.gz")

    script:
    """
    find . -name "*_R1_*.fastq.gz" | sort | xargs cat
    find . -name "*_R2_*.fastq.gz" | sort | xargs cat
    """
}

process FASTQC{
    container 'staphb/fastqc:0.11.9'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.zip"

    script:
    """
    fastqc --threads 2 --nogroup --extract ${reads}
    """
}

process SEQTK_SUBSAMPLE {
    container 'staphb/seqtk:1.3'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R?_subset.fastq.gz")

    script:
    """
    seqtk sample -s 12345 ${reads[0]} 30000000 | gzip - > ${sample_id}_R1_subset.fastq.gz
    seqtk sample -s 12345 ${reads[1]} 30000000 | gzip - > ${sample_id}_R2_subset.fastq.gz
    """
}

process FASTP {
    container 'biocontainers/fastp:v0.20.1_cv1'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R?_trimmed.fastq.gz")
    path "*.html"
    path "*.log"

    script:
    """
    fastp \
        -i ${reads[0]} \
        -I ${reads[1]} \
        -o ${sample_id}_R1_trimmed.fastq.gz \
        -O ${sample_id}_R2_trimmed.fastq.gz \
        -q 20 \
        -u 20 \
        -3 -W 5 \
        -x -g \
        -l 50 \
        -c 
        -h ${sample_id}_fastp.html \
        &> ${sample_id}_fastp.log
    """
}

process BWAMEM {
    container 'staphb/bwa:0.7.17'
    cpus 2

    input:
    tuple val(sample_id), path(reads)
    path reference

    output:
    tuple val(sample_id), stdout

    script:
    """
    bwa mem \
        -R '@RG\tID:${sample_id}\tDS:KAPA_TE\tPL:ILLUMINA\tLB:${sample_id}\tSM:${sample_id}' \
        ${reference} \
        -t 2 \
        -M ${reads}
    """
}

process GATK_CONVERT_FIX_SORT_BAM {
    container 'broadinstitute/gatk:4.2.0.0'

    input:
    tuple val(sample_id), path(sam_file)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam")

    script:
    """
    gatk SamFormatConverter \
        -I ${sam_file} \
        -O ${sample_id}_bwa.bam

    gatk FixMateInformation \
        -I ${sample_id}_bwa.bam \
        -O ${sample_id}_fixmate.bam

    gatk SortSam \
        -I ${sample_id}_fixmate.bam \
        -O ${sample_id}_sorted.bam \
        -SO coordinate \
        --CREATE_INDEX true
    """
}

process GATK_MAPPING_METRICS {
    container 'broadinstitute/gatk:4.2.0.0'

    input:
    tuple val(sample_id), path(bam_file)
    path reference

    output:
    path "${sample_id}_alignment_metrics.txt"

    script:
    """
    gatk CollectAlignmentSummaryMetrics \
        --INPUT ${bam_file} \
        --OUTPUT ${sample_id}_alignment_metrics.txt \
        --REFERENCE_SEQUENCE ${reference} \
        --METRIC_ACCUMULATION_LEVEL ALL_READS \
        --VALIDATION_STRINGENCY LENIENT
    """
}

process SAMTOOLS_MAPPING_RATES {
    container 'staphb/samtools:1.12'
    input:
    tuple val(sample_id), path(bam_file)

    output:
    path "${sample_id}_flagstat.txt"

    script:
    """
    samtools flagstat ${bam_file} > ${sample_id}_flagstat.txt
    """
}

process GATK_MARK_DUPLICATES {
    container 'broadinstitute/gatk:4.2.0.0'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("*.bam"), path("*.bai")

    script:
    """
    gatk MarkDuplicates \
        -I ${bam_file} \
        -O ${sample_id}_mdups.bam \
        --METRICS_FILE ${sample_id}_markduplicates_metrics.txt \
        --VALIDATION_STRINGENCY LENIENT \
        --REMOVE_DUPLICATES false \
        --ASSUME_SORTED true \
        --CREATE_INDEX true
    """
}

process GATK_INSERT_SIZE_DISTRIBUTION {
    container 'broadinstitute/gatk:4.2.0.0'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    path "*.txt"
    path "*.pdf"

    script:
    """
    gatk CollectInsertSizeMetrics \
        -I ${bam_file} \
        -O ${sample_id}_insert_size_metrics_sorted.txt
        -H ${sample_id}_insert_size_plot_sorted.pdf \
        --VALIDATION_STRINGENCY LENIENT \

    """
}

process GATK_COUNT_ON_TARGET_READS {
    container 'broadinstitute/gatk:4.2.0.0'

    input:
    tuple val(sample_id), path(bam_file)
    path reference
    path capture_target_bed_file

    output:
    path "ontarget_reads_sorted.txt"

    script:
    """
    gatk CountReads \
        -R ${reference} \
        -I ${bam_file} \
        -L ${capture_target_bed_file} \
        --read-filter MappedReadFilter \
        --read-filter NotSecondaryAlignmentReadFilter > ontarget_reads_sorted.txt
    """
}

process GATK_COLLECT_HS_METRICS {
    container 'broadinstitute/gatk:4.2.0.0'

    input:
    tuple val(sample_id), path(bam_file), path(bai_file)
    path reference
    path bait_interval_list

    output:
    path "*.txt"

    script:
    """
    gatk CollectHsMetrics \
        --INPUT ${bam_file} \
        --OUTPUT ${sample_id}_hs_metrics_sorted.txt \
        --BAIT_INTERVALS ${bait_interval_list} \
        --BAIT_SET_NAME DESIGN \
        --TARGET_INTERVALS ${bait_interval_list} \
        --METRIC_ACCUMULATION_LEVEL ALL_READS \
        --REFERENCE_SEQUENCE ${reference} \
        --VALIDATION_STRINGENCY LENIENT \
        --COVERAGE_CAP 100000 \
        --PER_BASE_COVERAGE ${sample_id}_per_base_coverage_sorted.txt
    """
}

process DEEPVARIANT {
    container 'google/deepvariant:0.9.0'

    input:
    tuple val(sample_id), path(bam_file), path(bai_file)
    path reference
    path capture_targets_bed_file
    path primary_targets_bed_file

    output:
    tuple val(sample_id), path("${sample_id}_deepvariant.vcf.{gz,gz.tbi}"), emit: vcf
    path "capture_primary_targets_union.bed", emit: bed
    path "*"

    script:
    """
    cat ${capture_targets_bed_file} ${primary_targets_bed_file} \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - > capture_primary_targets_union.bed

    /opt/deepvariant/bin/run_deepvariant \
        --model_type=WES \
        --ref="${reference}" \
        --reads="${bam_file}" \
        --regions="DESIGN_capture_primary_targets_union.bed" \
        --output_vcf="${sample_id}_deepvariant.vcf.gz" \
        --output_gvcf="${sample_id}_deepvariant.gvcf.gz" \
        --num_shards=$(nproc)
    """
}

process HAPPY {
    container 'pkrusche/hap.py:v0.3.9'
    input:
    tuple val(sample_id), path(vcf_file)
    path reference
    path gold_standard_vcf
    path confident_bed
    path capture_primary_target_bed

    output:
    path "*happy*"

    script:
    """
     /opt/hap.py/bin/hap.py \
        ${gold_standard_vcf} \
        ${vcf_file} \
        -o "${sample_id}_happy" \
        --gender auto \
        --threads 1 \
        -f ${confident_bed} \
        -T ${capture_primary_target_bed} \
        -r ${reference}
    """
}

/*
 * WORKFLOW
 */

workflow {
    CAT_FASTQ( ch_input_fastq )
    ch_fastq = CAT_FASTQ.out
    
    FASTQC( ch_fastq )
    
    SEQTK_SUBSAMPLE( ch_fastq )
    ch_subset_fastq = SEQTK_SUBSAMPLE.out

    FASTP( ch_subset_fastq )
    ch_trimmed_fastq = FASTP.out

    BWAMEM( ch_trimmed_fastq, ch_reference )
    ch_bwa = BWAMEM.out
    
    GATK_CONVERT_FIX_SORT_BAM( ch_bwa )
    ch_sorted_bam = GATK_CONVERT_FIX_SORT_BAM.out
    
    GATK_MAPPING_METRICS( ch_sorted_bam, ch_reference )
    
    SAMTOOLS_MAPPING_RATES( ch_sorted_bam )
    
    GATK_MARK_DUPLICATES( ch_sorted_bam )
    ch_mdups_bam = GATK_MARK_DUPLICATES.out
    
    GATK_INSERT_SIZE_DISTRIBUTION( ch_mdups_bam )
    
    GATK_COUNT_ON_TARGET_READS( 
        ch_mdups_bam, 
        ch_reference, 
        ch_capture_targets_bed 
    )
    
    GATK_COLLECT_HS_METRICS( 
        container 'broadinstitute/gatk:4.2.0.0'

        ch_mdups_bam, 
        ch_reference, 
        ch_bait_interval_list 
    )
    
    DEEPVARIANT( 
        ch_mdups_bam, 
        ch_reference, 
        ch_capture_targets_bed, 
        ch_primary_targets_bed
    )
    ch_deepvariant_vcf = DEEPVARIANT.out.vcf
    ch_deepvariant_bed = DEEPVARIANT.out.bed

    if (params.happy) {
        HAPPY( 
            ch_deepvariant_vcf,
            ch_reference,
            ch_gold_standard_vcf,
            ch_high_confidence_regions_bed,
            ch_deepvariant_bed
        )
    }
}