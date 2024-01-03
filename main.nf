println "*****************************************************"
println "*  Nextflow Malaria Variant analysis pipeline       *"
println "*  A Nextflow wrapper pipeline                      *"
println "*  Written by WEHI Research Computing Platform      *"
println "*  research.computing@wehi.edu.au                   *"
println "*                                                   *"
println "*****************************************************"
println " Required Pipeline parameters                        "
println "-----------------------------------------------------"
println "Input Sequence Path: $params.input_seq_path          "
println "Group Key file     : $params.input_group_key_file    "
println "Sample Key file    : $params.input_sample_key_file   "
println "Output directory   : $params.outdir                  " 
println "*****************************************************"

// include modules
include {  BWA  } from './modules/bwa.nf'
include { Index } from './modules/bwa.nf'
include { Merge } from './modules/bwa.nf'
include { BCF } from './modules/bcf.nf'
include { WriteBamList } from './modules/bcf.nf'
include{ FastQC } from './modules/qc.nf'
include{ MultiQC } from './modules/qc.nf'

process foo {
    debug true
    input:
    tuple val(sampleId), val(groupId),val(fastqbase),val(ref)
    output:
    stdout

    script:
    """
    echo $sampleId , $groupId,$fastqbase,$ref
    sleep 10
    """
}

workflow {
    
    Channel.fromPath(params.input_sample_key_file,checkIfExists:true)
    .ifEmpty{
        error("""
        No samples could be found in sample key file! Please check your sample key directory path
        is correct. 
        """)
    }
    .splitCsv(header:true,sep:'\t')
    .map { row -> 
        //groupId	sampleId	fastqbase	ref
        if (row.ref =="3D7") ref=params.ref3D7_path+"/PlasmoDB-52_Pfalciparum3D7_Genome"
        else if (row.ref =="Dd2") ref=params.refDd2_path_path+"/PlasmoDB-57_PfalciparumDd2_Genome"

        return tuple(row.sampleId, row.groupId, row.fastqbase,ref) 
    }
    .set{samplekey_ch}
    
    Channel.fromPath(params.input_group_key_file,checkIfExists:true)
    .ifEmpty{
        error("""
        No samples could be found in group key file! Please check your sample key directory path
        is correct. 
        """)
    }
    .splitCsv(header:true,sep:'\t')
    .map { row -> 
        //groupId	ref parentId	parentbamlist
        def parentbamlist = row.parentbamlist.replace(',', ' ')
        def mergedparent=row.parentId+".bam"
        if (row.ref =="3D7") ref=params.ref3D7_path+"/PlasmoDB-52_Pfalciparum3D7_Genome"
        else if (row.ref =="Dd2") ref=params.refDd2_path_path+"/PlasmoDB-57_PfalciparumDd2_Genome"

        return tuple(row.groupId,mergedparent, parentbamlist, ref) 
    }
    .set{groupkey_ch}

    Channel.fromPath(params.input_group_key_file,checkIfExists:true)
    .ifEmpty{
        error("""
        No samples could be found in group key file! Please check your sample key directory path
        is correct. 
        """)
    }
    .splitCsv(header:true,sep:'\t')
    .map { row -> 
        //groupId	ref parentId	parentbamlist
        def ref=""
        def prefix=""
        if (row.ref =="3D7") {
            ref=params.ref3D7_path+"/PlasmoDB-52_Pfalciparum3D7_Genome.fasta"
            prefix="_3D7ref"
        }        
        else if (row.ref =="Dd2") {
            ref=params.refDd2_path_path+"/PlasmoDB-57_PfalciparumDd2_Genome.fasta"
            prefix="_Dd2ref"
        }
        
        return tuple(row.groupId, ref,prefix) 
    }
    .set{groupkey_bcf_ch}   
    
    sam_ch=BWA(samplekey_ch)
    bam_ch=Index(sam_ch)


    merged_ch=Merge(groupkey_ch,bam_ch.bamnodup.collect())
    bamlist_ch=WriteBamList(bam_ch.bamnodup.collect())
    bcf_ch=BCF(groupkey_bcf_ch,bamlist_ch,bam_ch.bamnodup.collect())
    
    MultiQC(FastQC(bam_ch.bamnodup.collect()).zip.collect{ it[1] }.ifEmpty([]))  
}
