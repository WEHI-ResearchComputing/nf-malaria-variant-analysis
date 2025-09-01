process Bcf{
    label 'Bcftools'
    tag "${groupId}"   
    publishDir "${params.outdir}/variants/snvs_indels", mode: 'copy'
    cache true
    input:
    tuple  val(parentId),val(groupId),path(ref),
            path(bamlist),
            path(bams),path(parentbams)
    
    output:
    tuple  val(groupId), path("${groupId}.snvs_indels.vcf")


    script:
    """
    bcftools mpileup -Ou --max-depth 800 --threads ${task.cpus}  \
    -a AD,DP -f ${ref}  \
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
    tuple  val(parentId),val(groupId),path(ref),
            val(bamfilenames), 
            path(bams),path(parentbams), path(jarfile)

    output:
    tuple val(groupId), path("${groupId}.GRIDSS.bam"), emit: bam
    tuple val(groupId), path("${groupId}.GRIDSS.vcf"), emit: vcf

    script:
    """
    gridss --reference ${ref} \
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
            path(vcf), path(script)
            
    output:
    tuple val(groupId),path("${groupId}.SV_high_and_low_confidence_somatic.vcf.bgz"), emit:vcf
    
    script:
    
    """
    
    parentcount=\$(echo ${parentbamlist} | awk -F' ' '{print NF}')
    samplecount=\$(wc -l < ${groupId}_bams.txt)
    tumourordinals=\$(seq -s \' \' \$(expr \$parentcount + 1) \$samplecount)
   
    Rscript --vanilla ${script} \
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
   
    tuple  val(groupId),val(parentId),path(refpath),val(bsref),
            path(mergedparent), 
            path(bams), 
            val(bins),path(script)

    output:
    tuple val(groupId), path("*.rds"), path("*.csv")

    script:
    """
    Rscript --vanilla ${script} \
        --samplegroup ${groupId} \
        --parentId ${parentId} \
        --bams "${bams}" \
        --parentbam ${mergedparent} \
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
    tuple  val(groupId),path(refpath),val(prefix),
            path(parentbai),path(parentbam), val(parentbamlist), path(vcf), 
            path(bams), path(bai),path(script)

    output:
    path "*.csv", emit: csv 
    path "*plus*.Qcrit*.vcf", emit: vcf

    script:
    """
    Rscript --vanilla ${script} \
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
    tuple  val(groupId),val(parentbamlist), path(vcf),
        path(script),
        path(scriptdir)

    output:
    tuple val(groupId), path("${groupId}*.vcf"), emit: vcf
    tuple val(groupId), path("${groupId}*.csv"), emit: csv 
    script:
    """
    Rscript --vanilla ${projectDir}/Rtools/malDrugR/gridss_majorityfilt.R \
        --scriptdir ${scriptdir} \
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
            path(rds),path(script)


    output:
    tuple val(groupId), path("*copynums*.pdf")

    script:
    """
    Rscript --vanilla ${script} \
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
            path(rds),path(script)

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
process genesROI {
    label 'Rfilter'    
    tag "${groupId}"   
    
    publishDir "${params.outdir}/variants/copynumPlots", mode: 'copy'
    
    input:
    tuple  path(refpath),val(prefix),path(script)

    output:
    path("*.csv")

    script:
    """
    Rscript --vanilla ${projectDir}/Rtools/malDrugR/genesROI.R \
        --refpath ${refpath} \
        --refstrain ${prefix} \
        --region "${params.genesRegion}"
    """
}


