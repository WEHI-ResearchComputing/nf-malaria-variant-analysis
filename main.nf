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
include { InstallR } from './modules/stvariant.nf'
include { FilterBcfVcf } from './modules/stvariant.nf'
include { MajorityFilter } from './modules/stvariant.nf'
include { RPlot } from './modules/stvariant.nf'
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
    bamlist_ch=WriteBamLists(Channel.fromPath(params.input_sample_key_file,checkIfExists:true),
                    Channel.fromPath(params.input_group_key_file,checkIfExists:true))
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
            bsref="BSgenome.Pfalciparum3D7.PlasmoDB.52"
        }        
        else if (row.ref =="Dd2") {
            refpath=params.refDd2_path_path
            ref=params.refDd2_path_path+"/PlasmoDB-57_PfalciparumDd2_Genome"
            bsref="BSgenome.PfalciparumDd2.PlasmoDB.57"
        }
        return tuple(row.groupId,ref,refpath,row.ref,bsref, row.parentId, parentbamlist) 
    }
    .set{groupkey_ch}

    //Explain: .bamnodup.groupTuple() group tuples by the first value(groupid) 
    //Explain: then combine with groupkey_ch based on groupid
    allbams_grouped_ch=groupkey_ch.join(bam_ch.bamnodup.groupTuple(), by: 0,remainder:true)
    bam_bygroup_ch=allbams_grouped_ch.map{tuple-> if (tuple[1]!=null)
                            return tuple
                        }
    parentbam_ch=allbams_grouped_ch.map{tuple-> if (tuple[1]==null)
        return tuple[2]
    }
    merged_ch=Merge(groupkey_ch.combine(parentbam_ch.toList().collect()))
    //-----------------------------------------------------------------

    //----------------QC tools------------------------------------------
    MultiQC(FastQC(bam_ch.bam_4qc.collect()).zip.collect().ifEmpty([]))  
    //-----------------------------------------------------------------

    //----------------BCF tools----------------------------------------
    bcf_ch=Bcf(bam_bygroup_ch.join(bamlist_ch.flatten()
                                .map{filepath ->
                                    def groupid = filepath.baseName.split('_bam')[0]  // Splits the filename and takes the first part
                                    return tuple(groupid, filepath)  // Returns a tuple of the groupid and the original filepath
                                },by:0).combine(parentbam_ch.collect().toList()))
    //-----------------------------------------------------------------

    //----------------------Gridss------------------------------------- 
    
    gridss_combined_ch=bam_bygroup_ch.join(bamlist_ch.flatten().map{filepath ->
        def filename = filepath.baseName  // Extracts the filename without the path and extension
        def groupid = filename.split('_bam')[0]  // Splits the filename and takes the first part
        def fileLines = filepath.readLines() // Reads the file lines into a list
        def fileContents = fileLines.join(' ') // Joins the lines with a space
        return tuple(groupid, fileContents)  // Returns a tuple of the groupid and the original filepath
    },by:0)
    sv_ch=Gridss(gridss_combined_ch.combine(parentbam_ch.collect().toList()))
    
    combined_sv_ch=groupkey_ch.join(sv_ch.vcf, by: 0)
    dummy_ch=InstallR()
    smfilter_ch=SomaticFilter(combined_sv_ch,bamlist_ch.collect(),dummy_ch)
    //-------------------------------------------------------------------
    //----------------------CopyNum--------------------------------------
    copynum_ch=RCopyNum(gridss_combined_ch.join(merged_ch,by:0),
        dummy_ch)
    //-------------------------------------------------------------------
    //----------------------filter BCF----------------------------------- 
    parent_index_ch=groupkey_ch
                    .join(bcf_ch, by:0)
                    .join(bam_ch.bai.groupTuple(),by:0,remainder:true)
                    .map{
                        tuple-> if (tuple[1]==null)
                                return tuple[2]
                    }.collect().toList()               

    FilterBcfVcf(
                groupkey_ch.join(bcf_ch, by:0)
                .join(bam_ch.bamnodup.groupTuple(),by:0)
                .join(bam_ch.bai.groupTuple(),by:0)
                .combine(parentbam_ch.collect().toList())
                .combine(parent_index_ch),
                dummy_ch
            )
 
    //-------------------------------------------------------------------
    //----------------------Majority filter----------------------------------- 
    MajorityFilter(groupkey_ch
                    .join(smfilter_ch.vcf,by:0),dummy_ch
                  )

    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    //----------------------Plot-----------------------------------------
    RPlot(groupkey_ch.join(copynum_ch,by:0))
    //-------------------------------------------------------------------
}
