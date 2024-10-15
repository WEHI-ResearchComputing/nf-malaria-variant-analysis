process FastQC {
    label 'Fastqc' // assumes resources are defined in a label
    tag "${groupId}"
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

process Mosdepth {
    label 'Mosdepth'
    tag "${groupId}"
    publishDir "${params.outdir}/QC", mode: 'copy'
    cache true
    input:
    path(bams)
    
    output:
    path("*.global.dist.txt"), emit: coverage


    script:
    """
    mosdepth -t ${task.cpus} --no-per-base --mapq 10 $groupId ${bams} 
    """
}

process MultiQC {
    label 'Multiqc'
    cache true
    publishDir "${params.outdir}/multiqc/", mode: 'copy'

    input:
    path  multiqc_files, stageAs: "?/*" 

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots

   
    script:
    """

    multiqc --config ${params.multiqc_config} . 
    """
}
