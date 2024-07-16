process Bwa {

    label 'BWAAlign'
    
    tag "${sampleId}"
    
    input:
    tuple val(fastqbase), val(sampleId),val(groupId), val(ref),
                path(fastqs)
    
  
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


