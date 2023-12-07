#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process GENERATE_READ_COUNT {
    publishDir "${params.out_dir}/readcount", mode: 'copy'

    maxForks "4"
    input:
    tuple val(sample_id), path(bam)
    output:
    tuple val(sample_id), path("${sample_id}.rc"), emit:result

    script:
    """
    bam_file=${bam}
    if [ -L ${bam} ]; then
        bam_file="\$(readlink -f ${bam})"
    fi

    ln -s \${bam_file}.bai ${bam}.bai


    /opt/readCounter -w 100000 -q 30 \
        -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
        ${bam} > ${sample_id}.rc

    """
}


process SINGLE_SAMPLE_GC {
    publishDir "${params.out_dir}/gc", mode: 'copy'

    maxForks "4"
    input:
    tuple val(sample_id), path(rc)
    output:
    tuple path("${sample_id}.alldata.txt"), emit:result

    script:
    """
    python /opt/sample_feature.py -i ${rc} -model /opt/model/model_Duke.txt -o ${sample_id}.alldata.txt
    """
}


process MERGE_ALL_SAMPLE {
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path(gc)
    output:
    path("yj_bin_merge_${params.outbin_alias}.txt")

    script:
    """
    python /opt/merge_feature.py -i ${params.out_dir}/gc -inbin ${params.inbin} -outbin ${params.outbin} -o       yj_bin_merge_${params.outbin_alias}.txt
    """

}
