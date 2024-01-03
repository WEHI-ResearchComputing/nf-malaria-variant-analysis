process BWA {

    label 'BWAAlign'
    
    tag "${sampleId}"
    
    input:
    tuple val(sampleId), val(groupId),val(fastqbase),val(ref)
    
  
    output:
    tuple val(sampleId), path("*.sam") 

    script:
    """
    echo "BWA Alignment"
    bwa mem -t ${task.cpus} -o ${sampleId}.sam \
            -R "@RG\\tID:${sampleId}\\tSM:${sampleId}\\tPL:ILLUMINA" \
            ${ref} \
            '${params.input_seq_path}/${fastqbase}_R1_001.fastq.gz' '${params.input_seq_path}/${fastqbase}_R2_001.fastq.gz' 
    """
}

process Index {
    label 'BWAAlign'
    tag "${sampleId}"
    publishDir "${params.outdir}/align", mode: 'copy'
    
    input:
    tuple val(sampleId), path(samfile) 
    
  
    output:
    val(sampleId), emit: sampleid 
    path("*_nodup.bam"), emit: bamnodup
    path("*_s.bam"), emit: bams
    path("*.metrics")

    script:
    """
    echo "SAM Indexing"
    samtools view -b ${samfile}  |  \
    samtools sort -o ${sampleId}_s.bam -O bam -@ 4 -
    samtools index ${sampleId}_s.bam
    MarkDuplicates I=${sampleId}_s.bam O=${sampleId}_nodup.bam M=${sampleId}_duplic.metrics \
                REMOVE_DUPLICATES=TRUE
    samtools index ${sampleId}_nodup.bam


    """
}

process Merge{
    label 'Samtools'
    tag "${mergedparent}"
    publishDir "${params.outdir}/align", mode: 'copy'
    
    input:
    tuple  val(groupId),val(mergedparent), val(parentbamlist), val(ref)
    path(bamnodup)
    
    output:
    path("${mergedparent}")

    script:
    """
    echo "SAM Merge"
    samtools merge ${mergedparent} ${parentbamlist}
    """
}