process Bwa {

    label 'BWAAlign'
    
    tag "${sampleId}"
    
    input:
    tuple val(groupId),val(sampleId), val(fastqbase),val(ref)
    
  
    output:
    tuple val(sampleId),val(groupId), path("*.sam") 

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
    tuple val(sampleId),val(groupId), path(samfile) 
    
  
    output:
    tuple val(groupId), val(sampleId), path("*_nodup.bam"), emit: bamnodup
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
    tuple  val(groupId),val(ref),val(refpath),val(prefix),val(bsref), val(parentId), val(parentbamlist)
    tuple val(groupId), val(sampleId), path(bamnodup)
    
    output:
    tuple val(groupId),path("${parentId}.bam")

    script:
    """
    echo "SAM Merge"
    samtools merge -f ${parentId}.bam ${parentbamlist}
    """
}

process WriteBamLists {
    label 'Write'

    publishDir "${params.outdir}", mode: 'copy'
    input:
    path(sample_key)

    output:
    path("*_bams.txt")
    
    script:
    """
    {
        read
        while IFS=\$'\\t' read -r groupId sampleId fastqbase ref; do
            echo \$sampleId"_nodup.bam" >> \$groupId"_bams.txt"
        done 
    }< ${sample_key}
     
    """
}