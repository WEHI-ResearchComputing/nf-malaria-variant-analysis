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
include {  Bwa } from './modules/bwa.nf'
include { WriteBamLists } from './modules/bwa.nf'
include { Index } from './modules/bwa.nf'
include { Merge } from './modules/bwa.nf'
include { Bcf } from './modules/stvariant.nf'
include { Gridss } from './modules/stvariant.nf'
include { SomaticFilter } from './modules/stvariant.nf'
include { RCopyNum } from './modules/stvariant.nf'
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
    
    //----------------Alignment-----------------------------------------
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
        def ref=""
        if (row.ref =="3D7") ref=params.ref3D7_path+"/PlasmoDB-52_Pfalciparum3D7_Genome"
        else if (row.ref =="Dd2") ref=params.refDd2_path_path+"/PlasmoDB-57_PfalciparumDd2_Genome"
        return tuple(row.groupId,row.sampleId, row.fastqbase,ref) 
    }
    .set{samplekey_ch}
    bamlist_ch=WriteBamLists(Channel.fromPath(params.input_sample_key_file,checkIfExists:true))
    sam_ch=Bwa(samplekey_ch)
    bam_ch=Index(sam_ch)
    //-----------------------------------------------------------------

    //----------------Merge&List---------------------------------------
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
        def ref=""
        def refpath=""
        def prefix=""
        def bsref=""
        if (row.ref =="3D7") {
            refpath=params.ref3D7_path
            ref=params.ref3D7_path+"/PlasmoDB-52_Pfalciparum3D7_Genome"
            prefix="_3D7ref"
            bsref="BSgenome.Pfalciparum3D7.PlasmoDB.52"
        }        
        else if (row.ref =="Dd2") {
            refpath=params.refDd2_path_path
            ref=params.refDd2_path_path+"/PlasmoDB-57_PfalciparumDd2_Genome"
            prefix="_Dd2ref"
            bsref="BSgenome.PfalciparumDd2.PlasmoDB.57"
        }
        return tuple(row.groupId,ref,refpath,prefix,bsref, row.parentId, parentbamlist) 
    }
    .set{groupkey_ch}
    merged_ch=Merge(groupkey_ch,bam_ch.nodup_bams.collect())
    //-----------------------------------------------------------------

    //----------------QC tools------------------------------------------
    MultiQC(FastQC(bam_ch.nodup_bams.collect()).zip.collect().ifEmpty([]))  
    //-----------------------------------------------------------------

    //----------------BCF tools----------------------------------------
    bcf_ch=Bcf(groupkey_ch,bamlist_ch.collect(),bam_ch.nodup_bams.collect())
    //-----------------------------------------------------------------

    //----------------------Gridss------------------------------------- 
    sv_ch=Gridss(groupkey_ch,bamlist_ch.collect(),
            bam_ch.bamnodup.collect(),
            bam_ch.bamnodup.collect().map { list ->
                                        list.join(' ')
                                        })
    
    //combined_ch=groupkey_ch.join(sv_ch.vcf, by: 0)
    //SomaticFilter(combined_ch,bamlist_ch.collect())
    //-----------------------------------------------------------------
    //----------------------CopyNum------------------------------------- 
    //RCopyNum(bam_ch,merged_ch.collect(),groupkey_ch,Channel.fromPath(params.input_sample_key_file),Channel.fromPath(params.input_group_key_file))
    //-----------------------------------------------------------------
}
