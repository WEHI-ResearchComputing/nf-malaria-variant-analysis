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
    val sampleId
    tuple val(groupId), path("*_nodup.bam"), emit: bamnodup
    tuple val(groupId), path("*_nodup*bai"), emit: bai
    path("*_nodup.bam"), emit: bam_4qc
    path("*_s.bam"), emit: s_bams
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
    tag "${groupId}-${parentId}"
    publishDir "${params.outdir}/align", mode: 'copy'
    
    input:
    tuple  val(groupId),val(ref),val(refpath),val(prefix),val(bsref), val(parentId), val(parentbamlist),
            path(parentbams)
    
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
    path(group_key)

    output:
    path("*_bams.txt")
    
    script:
    """
    {
        read
        while IFS=\$'\\t' read -r groupId	ref	parentId	parentbamlist; do
            IFS=',' read -r -a array <<< "\$parentbamlist"
            # Iterate over the array
            for item in "\${array[@]}"; do
                echo "\$item">> \$groupId"_bams.txt"
            done
        done 
    }< ${group_key}
    {
        read
        while IFS=\$'\\t' read -r groupId sampleId fastqbase ref; do
            echo \$sampleId"_nodup.bam" >> \$groupId"_bams.txt"
        done 
    }< ${sample_key}
    
     
    """
}