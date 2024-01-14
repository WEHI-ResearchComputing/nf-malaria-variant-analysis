

process Bcf{
    label 'Bcftools'
    tag "${groupId}"   
    publishDir "${params.outdir}/variants/snvs_indels", mode: 'copy'

    input:
    tuple  val(groupId),val(ref),val(refpath),val(prefix),val(bsref), val(parentId), val(parentbamlist),
            path(bams),
            path(bamlist),
            path(parentbams)
    

    
    output:
    path("${groupId}${prefix}.vcf")


    script:
    """
    bcftools mpileup -Ou --max-depth 800 --threads ${task.cpus}  \
    -f ${ref}.fasta  \
    --bam-list ${bamlist}  | \
    bcftools call --ploidy 1 --threads ${task.cpus} -mv -Ov  \
    --output ${groupId}${prefix}.vcf -
    """
}

process Gridss{
    label 'Gridss' 
    label 'RGridss'   
    tag "${groupId}"   
    publishDir "${params.outdir}/variants/SVs", mode: 'copy'

    input:
    tuple  val(groupId),val(ref),val(refpath),val(prefix),val(bsref), val(parentId), val(parentbamlist),
            path(bams), 
            val(bamfilenames),
            path(parentbams)

    output:
    tuple val(groupId), path("${groupId}.bam"), emit: bam
    tuple val(groupId), path("${groupId}.vcf"), emit: vcf

    script:
    """
    gridss --reference ${ref}.fasta  \
    --jar ${params.gridss_jar_path} --assembly ${groupId}.bam  \
    --workingdir . --threads ${task.cpus} --skipsoftcliprealignment \
    --output ${groupId}.vcf  ${bamfilenames}
    
    
    """
}

process SomaticFilter{
    label 'Rfilter' 

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
    echo \$CONDA_PREFIX
    parentcount=\$(echo ${parentbamlist} | awk -F' ' '{print NF}')
    samplecount=\$(wc -l < ${groupId}_bams.txt)
    tumourordinals=\$(seq -s \' \' \$(expr \$parentcount + 1) \$samplecount)
    echo "Installing R packages."
    Rscript ${projectDir}/bin/installR.R
    echo "All R packages installed. Start Somatic filter."
    Rscript ${projectDir}/bin/gridss_assets/gridss_somatic_filter.R \
        --input ${groupId}.vcf \
        --fulloutput ${groupId}_high_and_low_confidence_somatic.vcf.bgz \
        --scriptdir ${projectDir}/bin/gridss_assets/  ##\$(dirname \$(which gridss_somatic_filter))\
        --ref ${bsref}  \
        --normalordinal 1  --tumourordinal \$tumourordinals
    bgzip -dc ${groupId}_high_and_low_confidence_somatic.vcf.bgz | awk  \
    \'/^#/ || \$7 ~ /^PASS\$/ || \$7 ~ /^imprecise\$/' >  \
    ${groupId}_high_and_imprecise.vcf

    echo "Output:"\$parentcount, \$samplecount, \$tumourordinals > output.txt
    
    """
}
process RCopyNum {
    label 'Gridss' 
    label 'Rfilter' 

    tag "${groupId}"   
    publishDir "${params.outdir}/variants/copynumfiles", mode: 'copy'
    
    input:
    val(dummy)
    // tuple  val(groupId),val(ref),val(refpath),val(prefix),val(bsref), val(parentId), val(parentbamlist),
    //         path(bams), 
    //         val(bamfilenames),
    //         path(parentbams)

    output:
    stdout
    //path("*.rds"), optional:true

    script:
    """
    echo ${dummy}
    Rscript ${projectDir}/bin/installR.R
    
    """
    // Rscript malDrugR/copynumQDNAseqParents.R \
    //     --samplegroup ${groupId} \
    //     --parentId ${parentId} \
    //     --bams ${bamfilenames} \
    //     --bin_in_kbases ${params.bin_in_kbases}
    // Rscript ${projectDir}/bin/Rlibs/forgeBSplasmo52.R

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


