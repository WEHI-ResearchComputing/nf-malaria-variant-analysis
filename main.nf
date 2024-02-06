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
println "Sample Group file  : $params.input_file              "
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
    //----------------Imput Preparation-----------------------------------------
    Channel.fromPath(params.input_file,checkIfExists:true)
    .ifEmpty{
        error("""
        No groups could be found in group key file! Please check your directory path
        is correct. 
        """)
    }
    .splitCsv(header:true,sep:'\t')
    .set{input_ch} // Emits a row [groupId	sampleId	fastqbase	ref	parentId]
    
    // .map { row -> 
    //     //groupId	sampleId	fastqbase	ref	parentId	parentbamlist
    //     def parentbamlist = row.parentbamlist.replace(',', ' ')
    //     def mergedparent=row.parentId+".bam"
    //     def ref=""
    //     def refpath=""
    //     def prefix=""
    //     def bsref=""
    //     if (row.ref =="3D7") {
    //         refpath=params.ref3D7_path
    //         ref=params.ref3D7_path+"/PlasmoDB-52_Pfalciparum3D7_Genome"
    //         bsref="BSgenome.Pfalciparum3D7.PlasmoDB.52"
    //     }        
    //     else if (row.ref =="Dd2") {
    //         refpath=params.refDd2_path_path
    //         ref=params.refDd2_path_path+"/PlasmoDB-57_PfalciparumDd2_Genome"
    //         bsref="BSgenome.PfalciparumDd2.PlasmoDB.57"
    //     }
    //     return tuple(row.groupId,row.sampleId, row.fastqbase, row.parentId, ref, refpath, row.ref,bsref, parentbamlist) 
    // }.set{input_ch}

    // input_ch.view()
    //-----------------------------------------------------------------

    //----------------Alignment-----------------------------------------
    bamlist_ch=WriteBamLists(Channel.fromPath(params.input_file,checkIfExists:true))
    bamlist_ch=bamlist_ch.flatten()
                        .map{
                                row -> tuple(row.baseName.split("_bams")[0],row)
                            } // Emits groupID, bamslist.txt
    input_ch.map{row -> 
                    ref=""
                    if (row.ref =="3D7") {
                        ref=params.ref3D7_path+"/PlasmoDB-52_Pfalciparum3D7_Genome"
                    }        
                    else if (row.ref =="Dd2") {
                        ref=params.refDd2_path_path+"/PlasmoDB-57_PfalciparumDd2_Genome"
                    }
                    
                    return tuple(row.sampleId,row.groupId,row.fastqbase,ref)
                }.set{bwa_input_ch}// Emit tuple val(sampleId),val(groupId), val(fastqbase),val(ref)

    sam_ch=Bwa(bwa_input_ch)
    bam_ch=Index(sam_ch)
    //-----------------------------------------------------------------

    //----------------Merge&List---------------------------------------
    input_ch.map{row -> 
                return row.parentId
            }.unique().combine(bam_ch.bamnodup,by:0)
            .map{row -> tuple(row[0],row[1], row[1].baseName+".bam")
            }.groupTuple()
            .map{row -> 
                    tuple(row[0],row[1], row[2].join(" "))
            }.set{parent_ch} // Emit tuple val(parentId),path(bams), val(parentlist)
    
    input_ch.map{row -> 
                return row.parentId
            }.unique().combine(bam_ch.bai,by:0)
            .groupTuple()
            .join(parent_ch)
            .set{parent_all_ch} // Emits tuple val(parentId),path(bai),path(bams), val(parentlist)
   
    merged_ch=Merge(parent_ch) // Emits tuple val(parentId),path(parentId.bam)
    //-----------------------------------------------------------------

    //----------------QC tools------------------------------------------
    MultiQC(FastQC(bam_ch.bamnodup.map{row->row[1]}.unique({it.baseName}).collect()).zip.collect().ifEmpty([]))  
    //-----------------------------------------------------------------

    //----------------BCF tools----------------------------------------
    // BCF Input Channel Emits parentId,groupId,ref,bamlist,bams,parentbams
    input_ch.map{row -> 
                    def ref=""
                    if (row.ref =="3D7") {
                        ref=params.ref3D7_path+"/PlasmoDB-52_Pfalciparum3D7_Genome"
                    }        
                    else if (row.ref =="Dd2") {
                        ref=params.refDd2_path_path+"/PlasmoDB-57_PfalciparumDd2_Genome"
                    }
                    return tuple(row.groupId,row.parentId,ref)
                }.unique()
                .join(bamlist_ch)
                .combine(bam_ch.bamnodup,by:0)
                .groupTuple(by:[0,1,2,3])
                .combine(parent_ch.map{row->tuple(row[1],row[0])},by:1).set{bcf_input_ch} // Emits val(parentId),val(groupId),val(ref), path(bamlist), path(bams),path(parentbams)
    bcf_ch=Bcf(bcf_input_ch)
    //---------------------------------------------------------------------

    //----------------------Gridss------------------------------------- 
    input_ch.map{row -> 
                    ref=""
                    if (row.ref =="3D7") {
                        ref=params.ref3D7_path+"/PlasmoDB-52_Pfalciparum3D7_Genome"
                    }        
                    else if (row.ref =="Dd2") {
                        ref=params.refDd2_path_path+"/PlasmoDB-57_PfalciparumDd2_Genome"
                    }
                    
                    return tuple(row.groupId,row.parentId,ref)
                }.unique()
                .join(bamlist_ch.map{gid,filepath ->
                        def fileLines = filepath.readLines() // Reads the file lines into a list
                        def fileContents = fileLines.join(' ') // Joins the lines with a space
                        return tuple(gid, fileContents)  // Returns a tuple of the groupid and the original filepath
                })
                .combine(bam_ch.bamnodup,by:0)
                .groupTuple(by:[0,1,2,3])
                .combine(parent_ch.map{row->tuple(row[1],row[0])},by:1).set{gridss_input_ch}// Emits val(parentId),val(groupId),val(ref), val(bamlistcontent), path(bams),path(parentbams)

    gridss_ch=Gridss(gridss_input_ch)
    
    dummy_ch=InstallR()
    input_ch.map{row -> 
                    bsref=""
                    if (row.ref =="3D7") {
                        bsref="BSgenome.Pfalciparum3D7.PlasmoDB.52"
                    }        
                    else if (row.ref =="Dd2") {
                        bsref="BSgenome.PfalciparumDd2.PlasmoDB.57"
                    }
                    
                    return tuple(row.parentId, row.groupId,bsref )
                }.unique()
                .combine(parent_ch,by:0).map{row -> tuple(row[1],row[2],row[5])}
                .join(bamlist_ch).join(gridss_ch.vcf).combine(dummy_ch)
                .set{sv_input_ch} //Emit val(groupId),val(bsref),val(parentbamlist), path(bamlist), path(vcf), val(dummy)

    sfilter_ch=SomaticFilter(sv_input_ch)
    //-------------------------------------------------------------------
    //----------------------CopyNum--------------------------------------
    input_ch.map{row -> 
                    refpath=""
                    bsref=""
                    if (row.ref =="3D7") {
                        refpath=params.ref3D7_path
                        bsref="BSgenome.Pfalciparum3D7.PlasmoDB.52"
                    }        
                    else if (row.ref =="Dd2") {
                        refpath=params.refDd2_path_path
                        bsref="BSgenome.PfalciparumDd2.PlasmoDB.57"
                    }
                    return tuple(row.parentId, row.groupId,refpath,bsref )
                }.unique()
                .combine(merged_ch,by:0)
                .map{row -> if (row[0]) tuple(row[1],row[0],row[2],row[3],row[4])}
                .join(bamlist_ch.map{gid,filepath ->
                         def fileLines = filepath.readLines() // Reads the file lines into a list
                         def fileContents = fileLines.join(' ') // Joins the lines with a space
                         return tuple(gid, fileContents)  // Returns a tuple of the groupid and the original filepath
                })
                .combine(bam_ch.bamnodup,by:0)
                .groupTuple(by:[0,1,2,3,4,5])
                .combine(dummy_ch)
                .set{copynum_input_ch}//Emit val(groupId),val(parentId),val(refpath),val(bsref),path(mergedparent), path(bamlistcontent), path(bams), val(dummy)
    
    copynum_ch=RCopyNum(copynum_input_ch)
    //-------------------------------------------------------------------
    //----------------------filter BCF----------------------------------- 
    input_ch.map{row -> 
                    refpath=""
                    if (row.ref =="3D7") {
                        refpath=params.ref3D7_path
                    }        
                    else if (row.ref =="Dd2") {
                        refpath=params.refDd2_path
                    }
                    return tuple(row.parentId,row.groupId,refpath,row.ref)
                }.unique()
                .combine(parent_all_ch,by:0)
                .map{row-> tuple(row[1],row[2],row[3],row[4],row[5],row[6])}
                .join(bcf_ch, by:0)
                .join(bam_ch.bamnodup.groupTuple(),by:0)
                .join(bam_ch.bai.groupTuple(),by:0)
                .combine(dummy_ch)
                .set{fbcf_ch}   // Emit val(groupId),val(refpath),val(prefix),
                                // path(parentbai),path(parentbam), val(parentbamlist), 
                                // path(vcf), 
                                // path(bams), path(bai), val(dummy)
    FilterBcfVcf(fbcf_ch )            
    //-------------------------------------------------------------------
    //----------------------Majority filter----------------------------------- 
    input_ch.map{row ->     
                    return tuple(row.groupId,row.parentId)
                }.unique()
                .combine(parent_ch.map{row->tuple(row[2],row[0])},by:1)
                .map{row->tuple(row[1],row[2])}
                .join(sfilter_ch.vcf)
                .combine(dummy_ch)
                .set{mjf_ch} // Emits val(groupId),val(parentlist), path(vcf), val(dummy)

    MajorityFilter(mjf_ch)
    //-------------------------------------------------------------------
    //----------------------Plot-----------------------------------------
    // Input Channel Emits val(groupId),val(parentId),val(refpath),val(prefix),path(rds)
                
    RPlot(input_ch.map{row ->     
                    refpath=""
                    bsref=""
                    if (row.ref =="3D7") {
                        refpath=params.ref3D7_path
                        bsref="BSgenome.Pfalciparum3D7.PlasmoDB.52"
                    }        
                    else if (row.ref =="Dd2") {
                        refpath=params.refDd2_path_path
                        bsref="BSgenome.PfalciparumDd2.PlasmoDB.57"
                    }
                    return tuple(row.groupId,row.parentId,refpath,bsref)
                }.unique()
                .join(copynum_ch,by:0))
    //-------------------------------------------------------------------
}
