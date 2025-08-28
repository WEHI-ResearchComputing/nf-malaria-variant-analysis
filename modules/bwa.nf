process LaneMerge {
    label 'MergeFastqFiles'
    tag "${fastqbase}"

    publishDir "${params.outdir}/merged_seqdata", mode: 'copy'

    input:
    tuple val(fastqbase), path(R1fastqs), path(R2fastqs)

    output:
    tuple val(fastqbase), path("*.fastq")

    script:
    """
    zcat ${R1fastqs.join(' ')} > ${fastqbase}_R1.fastq
    zcat ${R2fastqs.join(' ')} > ${fastqbase}_R2.fastq
    """
}

process Bwa {

    label 'BWAAlign'
    
    tag "${sampleId}"
    
    input:
    tuple val(fastqbase), val(sampleId),
            val(groupId), val(ref),
            path(refpath), path(fastqs)
    
  
    output:
    tuple val(sampleId),val(groupId), path("*.sam") 

    script:
    """
    bwa mem -t ${task.cpus} -o ${sampleId}.sam \
            -R "@RG\\tID:${sampleId}\\tSM:${sampleId}\\tPL:ILLUMINA" \
            ${ref} \
            '${fastqs[0]}' '${fastqs[1]}' 
    """
}

process Index {
    label 'IndexDedup'
    tag "${sampleId}"
    publishDir "${params.outdir}/align", mode: 'copy'
    
    input:
    tuple val(sampleId),val(groupId), path(samfile) 
  
    output:
    val sampleId
    tuple val(groupId), path("*_nodup.bam"), emit: bamnodup
    tuple val(groupId), path("*_nodup*bai"), emit: bai
    tuple val(sampleId), path("*_nodup.bam"), path("*_nodup*bai"), emit: bysampleid

    script:
    """
    samtools view -b ${samfile}  |  \
    samtools sort -o ${sampleId}_s.bam -O bam -@ 4 -
    samtools index ${sampleId}_s.bam
    picard MarkDuplicates I=${sampleId}_s.bam O=${sampleId}_nodup.bam M=${sampleId}_duplic.metrics \
                REMOVE_DUPLICATES=TRUE
    samtools index ${sampleId}_nodup.bam


    """
}

process Merge{
    label 'Samtools'
    tag "${parentId}"
    publishDir "${params.outdir}/align", mode: 'copy'
    
    input:
    tuple  val(parentId), path(parentbams), val(parentbamlist)
            
    
    output:
    tuple val(parentId),path("${parentId}.bam")

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
    path(inputfile)

    output:
    path("*_bams.txt"), emit:bams

    script:
    """
    sed -i -e '\$a\\' ${inputfile}
    declare -A assocArray

    {
        read
        while IFS=\$'\\t' read -r groupId sampleId fastqbase ref parentId; do
            assocArray[\$groupId]=\${parentId//[\$'\\t\\r\\n ']}
        done 
    }< 	${inputfile}

    {
        read
        while IFS=\$'\\t' read -r groupId sampleId fastqbase ref parentId; do
            for gid in \${!assocArray[@]}
            do
            if [[ "\${assocArray[\${gid}]}" == "\$groupId" ]]; then
                    echo \$sampleId"_nodup.bam" >> \$gid"_bams.txt"
                fi
                
            done
        done
    }< ${inputfile} 
    {
        read
        while IFS=\$'\\t' read -r groupId sampleId fastqbase ref parentId; do
            echo \$sampleId"_nodup.bam" >> \$groupId"_bams.txt"
        done 
    }< 	${inputfile}

    """
}