process FASTQC{
    input:
    tuple val(id), path(reads)

    output:
    path("sample_R1_fastqc.zip")
    path("sample_R2_fastqc.zip")
    tuple val(id), path("sample_R1.fastq.gz"), path("sample_R2.fastq.gz")

    script:
    """
    gunzip -c sample_R1.fastq.gz > sample_R1.fastq
    gunzip -c sample_R2.fastq.gz > sample_R2.fastq

    fastqc --threads 2 --nogroup --extract sample_R1.fastq sample_R2.fastq
    """
}

process SEQTK {

    output:
    tuple val(id), path("sample_R1_subset.fastq"), path("sample_R1_subset.fastq")

    script:
    """
    seqtk sample -s 12345 sample_R1.fastq 30000000 | gzip - > sample_R1_subset.fastq
    seqtk sample -s 12345 sample_R2.fastq 30000000 | gzip - > sample_R2_subset.fastq
    """
}

process FASTP {
    script:
    """
    fastp \
        -i sample_R1_subset.fastq \
        -I sample_R2_subset.fastq \
        -o sample_R1_subset_trimmed.fastq \
        -O sample_R2_subset_trimmed.fastq \
        -q 20 \
        -u 20 \
        -3 -W 5 \
        -x -g \
        -l 50 \
        -c 
        -h sample_fastp.html \
        -w
    """
}

process BWAMEM {
    input:
    path reference
    tuple val(id), path(reads)

    script:
    """
    bwa mem \
        -R '@RG\tID:1\tDS:KAPA_TE\tPL:ILLUMINA\tLB:SAMPLE\tSM:SAMPLE' \
        ${reference} \
        -t 2 \
        -M ${reads}
    """
}