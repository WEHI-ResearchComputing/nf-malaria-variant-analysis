

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
    tuple  val(groupId), path("${groupId}.vcf")


    script:
    """
    bcftools mpileup -Ou --max-depth 800 --threads ${task.cpus}  \
    -f ${ref}.fasta  \
    --bam-list ${bamlist}  | \
    bcftools call --ploidy 1 --threads ${task.cpus} -mv -Ov  \
    --output "${groupId}.vcf" -
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

process InstallR{
    label 'Rfilter' 

    output:
    val("1"), emit: dummy 
    
    script:
    
    """
    echo \$CONDA_PREFIX
    echo "Installing R packages."
    Rscript ${projectDir}/bin/installR.R    
    """
}

process SomaticFilter{
    label 'Rfilter' 

    tag "${groupId}"   
    publishDir "${params.outdir}/variants/SVs", mode: 'copy'
    
    input:
    tuple  val(groupId),val(ref),val(refpath),val(prefix),val(bsref), val(parentId), val(parentbamlist), path(vcf)
    path(bamlist)
    val(dummy)

    output:
    tuple val(groupId),path("${groupId}_high_and_imprecise.vcf"), emit:vcf
    path("output.txt")
    
    script:
    
    """
    echo \$CONDA_PREFIX
    parentcount=\$(echo ${parentbamlist} | awk -F' ' '{print NF}')
    samplecount=\$(wc -l < ${groupId}_bams.txt)
    tumourordinals=\$(seq -s \' \' \$(expr \$parentcount + 1) \$samplecount)
   
    echo "Start Somatic filter."
    Rscript ${projectDir}/bin/gridss_assets/gridss_somatic_filter.R \
        --input ${groupId}.vcf \
        --fulloutput ${groupId}_high_and_low_confidence_somatic.vcf.bgz \
        --scriptdir ${projectDir}/bin/gridss_assets/  ##\$(dirname \$(which gridss_somatic_filter))\
        --ref ${bsref}  \
        --normalordinal 1  --tumourordinal \$tumourordinals
    module load gzip
    gzip -dc ${groupId}_high_and_low_confidence_somatic.vcf.bgz.bgz | awk  \
    \'/^#/ || \$7 ~ /^PASS\$/ || \$7 ~ /^imprecise\$/' >  \
    ${groupId}_high_and_imprecise.vcf

    echo "Output:"\$parentcount, \$samplecount, \$tumourordinals > output.txt
    
    """
}
process RCopyNum {
    label 'Rfilter' 

    tag "${groupId}"   
    publishDir "${params.outdir}/variants/copynumfiles", mode: 'copy'
    
    input:
   
    tuple  val(groupId),val(ref),val(refpath),val(prefix),val(bsref), val(parentId), val(parentbamlist),
            path(bams), 
            val(bamfilenames),
            path(parentbams)
    val(dummy)

    output:
    tuple val(groupId), path("*.rds")

    script:
    """
    Rscript ${projectDir}/bin/malDrugR/copynumQDNAseqParents_mod.R \
        --samplegroup ${groupId} \
        --parentId ${parentId} \
        --bams "${bamfilenames}" \
        --bin_in_kbases ${params.bin_in_kbases} \
        --refDir ${refpath}\
        --bsref ${bsref}
    
    """

}

process FilterBcfVcf {
    label 'Rfilter'  
    tag "${groupId}"   
    publishDir "${params.outdir}/variants/SVs", mode: 'copy'

    input:
    tuple  val(groupId),val(ref),val(refpath),val(prefix),val(bsref), val(parentId), val(parentbamlist), 
            path(vcf),
            path(bams),
            path(bai),
            path (parentbams),
            path (parentindices)
    val(dummy)

    output:
    path "*.tsv", emit: tsv 
    path "*plus*.Qcrit*.vcf", emit: vcf

    script:
    """
    
    Rscript ${projectDir}/bin/malDrugR/filterBcfVcf_mod.R \
        --samplegroup ${groupId} \
        --refpath ${refpath} \
        --refstrain ${prefix} \
        --QUALcrit ${params.qualcrit} \
        --parentlist "${parentbamlist}"\
        --critsamplecount ${params.critsamplecount} 
        
    """
}

process MajorityFilter {
    label 'Rfilter'    
    tag "${groupId}"   
    publishDir "${params.outdir}/variants/SVs", mode: 'copy'


    input:
    tuple  val(groupId),val(ref),val(refpath),val(prefix),val(bsref), val(parentId), val(parentbamlist), 
        path(vcf)
    val(dummy)
    
    output:
    tuple val(groupId), path("${groupId}_all_*.vcf"), emit: vcf
    tuple val(groupId), path("${groupId}*.tsv"), emit: txt 
 

    script:
    """
    Rscript ${projectDir}/bin/malDrugR/gridss_majorityfilt_mod.R \
        --scriptdir ${projectDir}/bin/gridss_assets/ \
        --samplegroup ${groupId} \
        --parentlist "${parentbamlist}" \
        --critsamplecount ${params.critsamplecount} 
        
    """
}

process RPlot {
    label 'Rfilter'    
    tag "${groupId}"   
    
    publishDir "${params.outdir}/variants/copynumPlots", mode: 'copy'
    
    input:
    tuple  val(groupId),val(ref),val(refpath),val(prefix),val(bsref), val(parentId), val(parentbamlist), 
        path(rds)

    output:
    tuple val(groupId), path("*.pdf")

    script:
    """
    Rscript ${projectDir}/bin/malDrugR/copynumPlots_mod.R \
        --samplegroup ${groupId} \
        --parentId ${parentId} \
        --strain ${prefix} \
        --refDir ${refpath} \
        --bin_in_kbases ${params.bin_in_kbases} 
    """
}


