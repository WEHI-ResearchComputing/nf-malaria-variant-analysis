process FastQC {

    label 'Fastqc' // assumes resources are defined in a label
    
    publishDir "${params.outdir}/QC", mode: 'copy'
    cache true
    input:
    path(bams)
    
    
    output:
    path("*.html"), emit: html
    path("*.zip") , emit: zip
    

    script:
    """
    fastqc -t ${task.cpus} --noextract ${bams} -o ./
    """
}

process MosDepth {
    label 'Mosdepth'
    tag "${sampleId}"
    publishDir "${params.outdir}/QC", mode: 'copy'
    cache true

    input:
    tuple  val(sampleId), path(bam), path(bai)

    output:
    path("*.global.dist.txt"), emit: txt


    script:
    """
    mosdepth -t ${task.cpus} --no-per-base --mapq 10 ${bam.baseName} ${bam} 
    """
}

process FlagStats {
    label 'Flagstats'

    tag "${sampleId}"
    publishDir "${params.outdir}/QC", mode: 'copy'
    cache true

    input:
    tuple val(sampleId), path(bam), path(bai)

    output:
    path("*.flagstats"), emit: txt


    script:
    """
    samtools flagstats ${bam} > ${bam.baseName}.flagstats
    """
}

process MultiQC {
    label 'Multiqc'
    cache true
    publishDir "${params.outdir}/multiqc/", mode: 'copy'

    input:
    path  multiqc_files, stageAs: "?/*" 
    path(multiqc_config)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots

   
    script:
    """

    multiqc --config ${params.multiqc_config} . 
    """
}
