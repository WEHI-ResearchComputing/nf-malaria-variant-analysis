process FastQC {
    label 'Fastqc' // assumes resources are defined in a label
    tag "${groupId}"
    publishDir "${params.outdir}/QC", mode: 'copy'

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

process MultiQC {
    label 'Multiqc'
    
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
