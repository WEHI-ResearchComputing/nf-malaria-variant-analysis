

process BCF{
    label 'Bcftools'
    tag "${groupId}"   
    publishDir "${params.outdir}/variants/snvs_indels", mode: 'copy'

    input:
    
    tuple  val(groupId),val(ref),val(prefix)
    path(bamlist)
    path(bamnodup)
    

    output:
    path("${groupId}${prefix}.bcf")
    path("${groupId}${prefix}.vcf")

    script:
    """
    bcftools mpileup -Ou --max-depth 800 --threads ${task.cpus}  \
    -f ${ref}  \
    --bam-list ${bamlist}  | \
    bcftools call --ploidy 1 --threads ${task.cpus} -mv -Ob  \
    --output ${groupId}${prefix}.bcf -
    
    bcftools view -Ov ${groupId}${prefix}.bcf > ${groupId}${prefix}.vcf
    """
}

process WriteBamList {

    publishDir "${params.outdir}/", mode: 'copy'

    input:
    path(list)

    output:
    path("bam_list.txt")
    
    script:
    """
    for file in \$(ls ${list}); do
        echo \${file} >> bam_list.txt
    done    
    """
}