process Bwa {

    label 'BWAAlign'
    
    tag "${sampleId}"
    
    input:
    tuple val(sampleId),val(groupId), val(fastqbase),val(ref)
    
  
    output:
    tuple val(sampleId),val(groupId), path("*.sam") 

    script:
    """
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
    path("*.metrics")

    script:
    """
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
            echo \$sampleId"_nodup.bam" >> \$groupId"_bams.txt"
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
    """
}