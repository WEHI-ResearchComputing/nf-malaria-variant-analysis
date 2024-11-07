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
    tuple  val(groupId), path("${groupId}.snvs_indels.vcf")


    script:
    """
    bcftools mpileup -Ou --max-depth 800 --threads ${task.cpus}  \
    -a AD,DP -f ${ref}.fasta  \
    --bam-list ${bamlist}  | \
    bcftools call --ploidy 1 --threads ${task.cpus} -mv -Ov  \
    --output "${groupId}.snvs_indels.vcf" -
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
    tuple val(groupId), path("${groupId}.GRIDSS.bam"), emit: bam
    tuple val(groupId), path("${groupId}.GRIDSS.vcf"), emit: vcf

    script:
    """
    gridss --reference ${ref}.fasta  \
    --jar ${params.gridss_jar_path} --assembly ${groupId}.GRIDSS.bam  \
    --workingdir . --threads ${task.cpus} --skipsoftcliprealignment \
    --output ${groupId}.GRIDSS.vcf  ${bamfilenames}
    
    
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
            path(vcf)
            
    output:
    tuple val(groupId),path("${groupId}.SV_high_and_low_confidence_somatic.vcf.bgz"), emit:vcf
    
    script:
    
    """
    
    parentcount=\$(echo ${parentbamlist} | awk -F' ' '{print NF}')
    samplecount=\$(wc -l < ${groupId}_bams.txt)
    tumourordinals=\$(seq -s \' \' \$(expr \$parentcount + 1) \$samplecount)
   
    Rscript --vanilla ${projectDir}/Rtools/gridss_assets/gridss_somatic_filter.R \
        --input ${groupId}.GRIDSS.vcf \
        --fulloutput ${groupId}.SV_high_and_low_confidence_somatic.vcf \
        --scriptdir ${projectDir}/Rtools/gridss_assets/  ##\$(dirname \$(which gridss_somatic_filter))\
        --ref ${bsref}  \
        --normalordinal 1  --tumourordinal \$tumourordinals
        
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
            val(bins)

    output:
    tuple val(groupId), path("*.rds"), path("*.csv")

    script:
    """
    Rscript --vanilla ${projectDir}/Rtools/malDrugR/copynumQDNAseq.R \
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
            path(bams), path(bai)

    output:
    path "*.csv", emit: csv 
    path "*plus*.Qcrit*.vcf", emit: vcf

    script:
    """
    Rscript --vanilla ${projectDir}/Rtools/malDrugR/filterBcfVcf.R \
        --samplegroup ${groupId} \
        --refpath ${refpath} \
        --refstrain ${prefix} \
        --QUALcrit ${params.BCFqualcrit} \
        --parentlist "${parentbamlist}"\
        --critsamplecount ${params.critsamplecount} 
        
    """
}

process FilterGridssSV {
    label 'Rfilter'    
    tag "${groupId}"   
    publishDir "${params.outdir}/variants/SVs", mode: 'copy'


    input:
    tuple  val(groupId),val(parentbamlist), path(vcf)
    
    output:
    tuple val(groupId), path("${groupId}*.vcf"), emit: vcf
    tuple val(groupId), path("${groupId}*.csv"), emit: csv 
    script:
    """
    Rscript --vanilla ${projectDir}/Rtools/malDrugR/gridss_majorityfilt.R \
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


