

process Bcf{
    label 'Bcftools'
    tag "${groupId}"   
    publishDir "${params.outdir}/variants/snvs_indels", mode: 'copy'

    input:
    tuple  val(groupId),val(ref),val(refpath),val(prefix),val(bsref), val(parentId), val(parentbamlist)
    path(bamlist)
    path(bamnodup)

    
    output:
    path("${groupId}${prefix}.bcf")
    path("${groupId}${prefix}.vcf")


    script:
    """
    bcftools mpileup -Ou --max-depth 800 --threads ${task.cpus}  \
    -f ${ref}.fasta  \
    --bam-list ${bamlist}  | \
    bcftools call --ploidy 1 --threads ${task.cpus} -mv -Ob  \
    --output ${groupId}${prefix}.bcf -
    
    bcftools view -Ov ${groupId}${prefix}.bcf > ${groupId}${prefix}.vcf
    """
}

process Gridss{
    label 'Gridss'    
    tag "${groupId}"   
    publishDir "${params.outdir}/variants/SVs", mode: 'copy'

    input:
    tuple  val(groupId),val(ref),val(refpath),val(prefix),val(bsref), val(parentId), val(parentbamlist)
    path(bamlist)
    path(bamnodup)
    val(bams)

    output:
    tuple val(groupId), path("${groupId}.bam"), emit: bam
    tuple val(groupId), path("${groupId}.vcf"), emit: vcf

    script:
    """
    gridss --reference ${ref}.fasta  \
    --jar ${params.gridss_jar_path} --assembly ${groupId}.bam  \
    --workingdir . --threads 8 --skipsoftcliprealignment \
    --output ${groupId}.vcf  ${bams}
    
    
    """
}

process SomaticFilter{
    label 'Gridss'  
    tag "${groupId}"   
    publishDir "${params.outdir}/variants/SVs", mode: 'copy'

    input:
    tuple  val(groupId),val(ref),val(refpath),val(prefix),val(bsref), val(parentId), val(parentbamlist), path(vcf)
    path(bamlist)


    output:
    path("*.vcf*")
    path("output.txt")
    stdout

    script:
    
    """
    parentcount=\$(echo ${parentbamlist} | awk -F' ' '{print NF}')
    samplecount=\$(wc -l < ${groupId}_bams.txt)
    tumourordinals=\$(seq -s \' \' \$(expr \$parentcount + 1) \$samplecount)
    gridss_somatic_filter \
        --input ${groupId}.vcf \
        --fulloutput ${groupId}_high_and_low_confidence_somatic.vcf.bgz \
        --scriptdir  \$(dirname \$(which gridss_somatic_filter))\
        --ref ${bsref}  \
        --normalordinal 1  --tumourordinal \$tumourordinals
    echo "Writing ${groupId}_high_and_imprecise.vcf"
    bgzip -dc ${groupId}_high_and_low_confidence_somatic.vcf.bgz | awk  \
    \'/^#/ || \$7 ~ /^PASS\$/ || \$7 ~ /^imprecise\$/' >  \
    ${groupId}_high_and_imprecise.vcf

    echo "Output:"\$parentcount, \$samplecount, \$tumourordinals > output.txt
    
    """
}
process RCopyNum {
    label 'Rbcf'  
    label 'Rscript'  

    tag "${groupId}"   
    publishDir "${params.outdir}/variants/copynumfiles", mode: 'copy'
    
    input:
    path(bam)
    tuple  val(groupId),val(ref),val(refpath),val(prefix),val(bsref), val(parentId), val(parentbamlist)
    path(samplekeyfile)
    path(groupkeyfile)
    

    output:
    path("*.rds")

    script:
    """
    Rscript ${projectDir}/bin/malDrugR/copynumQDNAseqParents.R \
        --samplegroup ${groupId} \
        --samplekeyfile ${samplekeyfile} \
        --groupkeyfile ${groupkeyfile} \
        --bin_in_kbases ${params.bin_in_kbases}
    """
}

process FilterBcfVcf {
    label 'Rbcf'    
    tag "${groupId}"   
    publishDir "${params.outdir}/variants/SVs", mode: 'copy'
    
    input:
    tuple  val(groupId),val(ref),val(refpath),val(prefix),val(bsref), val(parentId), val(parentbamlist)
    path(groupkeyfile)
    path(vcfs)

    output:
    path "*.tsv"
    path "*.vcf"

    script:
    """
    Rscript ${projectDir}/bin/malDrugR/filterBcfVcf.R \
        --samplegroup $samplegroup \
        --vardir ./ \
        --groupkeyfile $groupkeyfile \
        --QUALcrit ${params.qualcrit} \
        --critsamplecount ${params.critsamplecount}
    """
}

process MajorityFilter {
    label 'Rbcf'    
    tag "${groupId}"   
    publishDir "${params.outdir}/variants/SVs", mode: 'copy'
    script:
    """
    """
}

process Plot {
    label 'Rbcf'    
    tag "${groupId}"   
    publishDir "${params.outdir}/variants/SVs", mode: 'copy'
    script:
    """
    """
}


