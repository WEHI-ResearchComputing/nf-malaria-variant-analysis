

process Bcf{
    label 'Bcftools'
    tag "${groupId}"   
    publishDir "${params.outdir}/variants/snvs_indels", mode: 'copy'
    cache true
    input:
    tuple  val(parentId),val(groupId),val(ref),
            path(bamlist),
            path(bams),path(parentbams)
    
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
    cache true  
    publishDir "${params.outdir}/variants/SVs", mode: 'copy'

    input:
    tuple  val(parentId),val(groupId),val(ref),
            val(bamfilenames), 
            path(bams),path(parentbams)

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
    Rscript --vanilla ${projectDir}/Rtools/installR.R    
    """
}

process SomaticFilter{
    label 'Rfilter' 
    cache true
    tag "${groupId}"   
    publishDir "${params.outdir}/variants/SVs", mode: 'copy'
    
    input:
    tuple val(groupId),val(bsref),val(parentbamlist), 
            path(bamlist), 
            path(vcf),
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
    Rscript --vanilla ${projectDir}/Rtools/gridss_assets/gridss_somatic_filter.R \
        --input ${groupId}.vcf \
        --fulloutput ${groupId}_high_and_low_confidence_somatic.vcf.bgz \
        --scriptdir ${projectDir}/Rtools/gridss_assets/  ##\$(dirname \$(which gridss_somatic_filter))\
        --ref ${bsref}  \
        --normalordinal 1  --tumourordinal \$tumourordinals
    module load gzip
    gzip -dc ${groupId}_high_and_low_confidence_somatic.vcf.bgz.bgz | awk  \
    \'/^#/ || \$7 ~ /^PASS\$/ || \$7 ~ /^imprecise\$/' >  \
    ${groupId}_high_and_imprecise.vcf

    echo "GRIDSS somatic filter used 1st line of ${groupId}_bams.txt as Normal sample" > output.txt
    echo "and lines \$tumourordinals as 'tumour' samples." >> output.txt
    
    """
}
process RCopyNum {
    label 'Rfilter'         
    cache true

    tag "${groupId}"   
    publishDir "${params.outdir}/variants/copynumfiles", mode: 'copy'
    
    input:
   
    tuple  val(groupId),val(parentId),val(refpath),val(bsref),
            path(mergedparent), 
            val(bamfilenames), 
            path(bams), 
            val(dummy),val(bins)

    output:
    tuple val(groupId), path("*.rds")

    script:
    """
    Rscript --vanilla ${projectDir}/Rtools/malDrugR/copynumQDNAseqParents_mod.R \
        --samplegroup ${groupId} \
        --parentId ${parentId} \
        --bams "${bamfilenames}" \
        --bin_in_kbases ${bins} \
        --refDir ${refpath}\
        --bsref ${bsref}
    
    """

}

process FilterBcfVcf {
    label 'Rfilter'  
    tag "${groupId}"   
    publishDir "${params.outdir}/variants/snvs_indels", mode: 'copy'

    input:
    tuple  val(groupId),val(refpath),val(prefix),
            path(parentbai),path(parentbam), val(parentbamlist), path(vcf), 
            path(bams), path(bai), val(dummy)

    output:
    path "*.tsv", emit: tsv 
    path "*plus*.Qcrit*.vcf", emit: vcf

    script:
    """
    Rscript --vanilla ${projectDir}/Rtools/malDrugR/filterBcfVcf_mod.R \
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
    tuple  val(groupId),val(parentbamlist), path(vcf), val(dummy)
    
    output:
    tuple val(groupId), path("${groupId}_all_*.vcf"), emit: vcf
    tuple val(groupId), path("${groupId}*.tsv"), emit: txt 
    script:
    """
    Rscript --vanilla ${projectDir}/Rtools/malDrugR/gridss_majorityfilt_mod.R \
        --scriptdir ${projectDir}/Rtools/gridss_assets/ \
        --samplegroup ${groupId} \
        --parentlist "${parentbamlist}" \
        --critsamplecount ${params.critsamplecount} 
        
    """
}

process RPlotFull {
    label 'Rfilter'    
    tag "${groupId}"   
    
    publishDir "${params.outdir}/variants/copynumPlots", mode: 'copy'
    
    input:
    tuple  val(groupId),val(parentId),
            path(rds)

    output:
    tuple val(groupId), path("*.pdf")

    script:
    """
    Rscript --vanilla ${projectDir}/Rtools/malDrugR/copynumPlotsFull.R \
        --samplegroup ${groupId} \
        --parentId ${parentId} \
        --bin_in_kbases ${params.bin_CNfull}  \
        --lowerbound_plot ${params.lowerbound_fullCN}  \
        --upperbound_plot ${params.upperbound_fullCN}
    """
}
process RPlotROI {
    label 'Rfilter'    
    tag "${groupId}"   
    
    publishDir "${params.outdir}/variants/copynumPlots", mode: 'copy'
    
    input:
    tuple  val(groupId),val(parentId),
            path(rds)

    output:
    tuple val(groupId), path("*.pdf")

    script:
    """
    Rscript --vanilla ${projectDir}/Rtools/malDrugR/copynumPlotsROI.R \
        --samplegroup ${groupId} \
        --parentId ${parentId} \
        --bin_in_kbases ${params.bin_CNroi} \
        --chrOI ${params.chr_CNroi} \
        --startROI ${params.start_CNroi} \
        --endROI ${params.end_CNroi}
    """
}


