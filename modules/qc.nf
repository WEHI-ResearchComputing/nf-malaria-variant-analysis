process FastQC {
    label 'Fastqc' // assumes resources are defined in a label

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

process BioBloomCategorizer {
    
    label 'Biobloom'

    tag "${sampleId}"
    publishDir "${params.outdir}/BioBloom", mode: 'copy'
    
    conda "${projectDir}/envs/biobloom.yaml"

    input:
    tuple val(sampleId), path(reads)
    path(biobloom_filters)
    path(biobloom_txts)

    output:
    path "*.tsv", emit: biobloom

    script:
    def read1, read2
    if (params.bb_first_n_reads > 0) {
        if (reads[0].toString().endsWith(".gz")) {
            read1 = "zcat ${reads[0]} | head -n ${4 * params.bb_first_n_reads} > temp_1.fq"
            read2 = "zcat ${reads[1]} | head -n ${4 * params.bb_first_n_reads} > temp_2.fq"
        } else {
            read1 = "head -n ${4 * params.bb_first_n_reads} ${reads[0]} > temp_1.fq"
            read2 = "head -n ${4 * params.bb_first_n_reads} ${reads[1]} > temp_2.fq"
        }
    } else {
        read1 = "ln -s ${reads[0]} temp_1.fq"
        read2 = "ln -s ${reads[1]} temp_2.fq"
    }

    """
    # extract first n reads
    ${read1}
    ${read2}

    biobloomcategorizer -t ${task.cpus} -e -p ${sampleId} -f "${biobloom_filters}" temp_1.fq temp_2.fq
    """
}

process BioBloomMaker {
    
    label 'Biobloom'

    tag "${ref_id}"
    publishDir "${params.outdir}/BioBloom", mode: 'copy'
    
    conda "${projectDir}/envs/biobloom.yaml"

    input:
    tuple val(ref_id), path(ref)

    output:
    path("*.bf"), emit: bb_filter
    path("*.txt"), emit: bb_txt
    
    script:
    """
    biobloommaker -t ${task.cpus} -p ${ref_id} ${ref}
    """
}

process SamtoolsStats {
    label 'SamtoolsStats'

    tag "${sampleId}"

    publishDir "${params.outdir}/samtools/", mode: 'copy'

    input:
    tuple val(sampleId), path(bamFile)

    output:
    tuple val(sampleId), path("*.stats"), emit: stats
    tuple val(sampleId), path("*.flagstat"), emit: flagstat

    script:
    """
    
    samtools stats \\
        --threads ${task.cpus} \\
        ${bamFile} > ${sampleId}.stats
    samtools flagstat \\
            --threads ${task.cpus} \\
            ${bamFile} > ${sampleId}.flagstat
    """
}

process Qualimap {
    label 'Qualimap' 
    
    tag "${sampleId}"
    publishDir "${params.outdir}/qualimap/", mode: 'copy'

    conda   "${projectDir}/envs/qualimap.yaml"

    input:
    tuple val(sampleId), path(bamFile)

    output:
    tuple val(sampleId), path("${sampleId}"), emit: results
    

    script:
    """
    unset DISPLAY
    mkdir -p tmp
    export _JAVA_OPTIONS="-Djava.io.tmpdir=./tmp"
    export JAVA_OPTS="-Djava.awt.headless=true -Xmx${task.memory.toMega()}M" 

    qualimap bamqc \\
        --java-mem-size=${task.memory.toMega()}M \\
        -bam ${bamFile} \\
        -outdir ${sampleId} \\
        -outformat HTML \\
        -nt ${task.cpus} 
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
