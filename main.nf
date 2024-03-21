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
include { RPlotFull } from './modules/stvariant.nf'
include { RPlotROI } from './modules/stvariant.nf'

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
    //----------------Input Preparation-----------------------------------------
    Channel.fromPath(params.input_file,checkIfExists:true)
    .ifEmpty{
        error("""
        No groups could be found in group key file! Please check your directory path
        is correct. 
        """)
    }
    .splitCsv(header:true,sep:'\t')
    .map{row ->
        ref=""
        refpath=""
        bsref=""
        if (row.ref =="3D7") {
            ref=params.ref3D7_path+"/PlasmoDB-52_Pfalciparum3D7_Genome"
            refpath=params.ref3D7_path
            bsref="BSgenome.Pfalciparum3D7.PlasmoDB.52"
        }        
        else if (row.ref =="Dd2") {
            ref=params.refDd2_path+"/PlasmoDB-57_PfalciparumDd2_Genome"
            refpath=params.refDd2_path
            bsref="BSgenome.PfalciparumDd2.PlasmoDB.57"
        }
        return tuple(row.groupId,row.sampleId,row.fastqbase,row.ref, row.parentId,ref,refpath,bsref)
    }
    .set{input_ch} // Emits 0->groupId,	1->sampleId, 2->fastqbase,	
                   // 3->ref_prefix, 4->parentId, 5->ref (path+name)
                   // 6->refpath , 7-> bsref
    //----------------Alignment-----------------------------------------
    bamlist_ch=WriteBamLists(Channel.fromPath(params.input_file,checkIfExists:true))
    bamlist_ch=bamlist_ch.flatten()
                        .map{
                                row -> tuple(row.baseName.split("_bams")[0],row)
                            } // Emits groupID, bamslist.txt
    Channel.fromFilePairs("${params.input_seq_path}/*_{,R}{1,2}*.{fq,fastq}{,.gz}", size: 2 )
            .ifEmpty {
                    error("""
                    No samples could be found! Please check whether your input directory
                    is correct, and that your samples match typical fastq paired end naming
                    convention(s).
                    """)
            }.map{ row-> tuple(row[0].split(/_R[0-9]{1}/)[0],row[1])}
            .set {fastq_input_channel}
    
    input_ch.map{row -> tuple(row[2],row[1],row[0],row[5])}
                .join(fastq_input_channel,by:0)
                .set{bwa_input_ch}// Emits tuple val(sampleId),val(groupId), val(fastqbase),val(ref),path(fastqs)
   
    sam_ch=Bwa(bwa_input_ch)
    bam_ch=Index(sam_ch)
    //----------------Merge&List---------------------------------------
    input_ch.map{row -> row[4]}
            .unique()
            .combine(bam_ch.bamnodup,by:0)
            .map{row -> tuple(row[0],row[1], row[1].baseName+".bam")
            }.groupTuple()
            .map{row -> 
                    tuple(row[0],row[1], row[2].join(" "))
            }
            .ifEmpty {
                    error("""
                    No parent samples found.
                    """)
            }
            .set{parent_ch} // Emits tuple val(parentId),path(bams), val(parentlist)
    
    input_ch.map{row -> row[4]}
            .unique().combine(bam_ch.bai,by:0)
            .groupTuple()
            .join(parent_ch)
            .ifEmpty {
                    error("""
                    No parent samples found.
                    """)
            }
            .set{parent_all_ch} // Emits tuple val(parentId),path(bai),path(bams), val(parentlist)
    
    merged_ch=Merge(parent_ch) // Emits tuple val(parentId),path(parentId.bam)
    //----------------QC tools------------------------------------------
    MultiQC(FastQC(bam_ch.bamnodup.map{row->row[1]}.unique({it.baseName}).collect()).zip.collect().ifEmpty([]))  
    //----------------BCF tools----------------------------------------
    //BCF Input Channel Emits parentId,groupId,ref,bamlist,bams,parentbams
    input_ch.map{row -> tuple(row[0],row[4],row[5])}
            .unique()
            .join(bamlist_ch)
            .combine(bam_ch.bamnodup,by:0)
            .groupTuple(by:[0,1,2,3])
            .combine(parent_ch.map{row->tuple(row[1],row[0])},by:1).set{bcf_input_ch} // Emits val(parentId),val(groupId),val(ref), path(bamlist), path(bams),path(parentbams)
    bcf_ch=Bcf(bcf_input_ch)
    //----------------------Gridss------------------------------------- 
    input_ch.map{row -> tuple(row[0],row[4],row[5])}
            .unique()
            .join(bamlist_ch.map{gid,filepath ->
                    def fileLines = filepath.readLines() 
                    def fileContents = fileLines.join(' ') 
                    return tuple(gid, fileContents)  
            })
            .combine(bam_ch.bamnodup,by:0)
            .groupTuple(by:[0,1,2,3])
            .combine(parent_ch.map{row->tuple(row[1],row[0])},by:1)
            .ifEmpty {
                error("""
                Input to gridss is empty.
                """)
            }
            .set{gridss_input_ch}// Emits val(parentId),val(groupId),val(ref), val(bamlistcontent), path(bams),path(parentbams)
    gridss_ch=Gridss(gridss_input_ch)
    
    input_ch.map{row -> tuple(row[4], row[0],row[7] )}
            .unique()
            .combine(parent_ch,by:0).map{row -> tuple(row[1],row[2],row[5])}
            .join(bamlist_ch).join(gridss_ch.vcf)
            .set{sv_input_ch} //Emit val(groupId),val(bsref),val(parentbamlist), path(bamlist), path(vcf)

    sfilter_ch=SomaticFilter(sv_input_ch)
    //----------------------CopyNum--------------------------------------
    input_ch.map{row -> tuple(row[4], row[0],row[6],row[7] )}
            .unique()
            .combine(merged_ch,by:0)
            .map{row -> if (row[0]) tuple(row[1],row[0],row[2],row[3],row[4])}
            .join(bamlist_ch.map{gid,filepath ->
                def fileLines = filepath.readLines() 
                def fileContents = fileLines.join(' ') 
                return tuple(gid, fileContents)  
            })
            .combine(bam_ch.bamnodup,by:0)
            .groupTuple(by:[0,1,2,3,4,5])
            .ifEmpty {
                error("""
                Input to CopyNum Analysis is empty.
                """)
            }
            .set{copynum_input_ch}//Emits val(groupId),val(parentId),val(refpath),val(bsref),path(mergedparent), path(bamlistcontent), path(bams), val(dummy)
    
    copynum_ch=RCopyNum(copynum_input_ch
                            .combine(Channel.fromList( [params.bin_CNroi, params.bin_CNfull] ))
                        )
    //----------------------filter BCF----------------------------------- 
    input_ch.map{row -> tuple(row[4],row[0],row[6],row[3])}
            .unique()
            .combine(parent_all_ch,by:0)
            .map{row-> tuple(row[1],row[2],row[3],row[4],row[5],row[6])}
            .join(bcf_ch, by:0)
            .join(bam_ch.bamnodup.groupTuple(),by:0)
            .join(bam_ch.bai.groupTuple(),by:0)
            .ifEmpty {
                error("""
                Input to filter is empty.
                """)
            }
            .set{fbcf_ch}   // Emits val(groupId),val(refpath),val(prefix),
                            // path(parentbai),path(parentbam), val(parentbamlist), 
                            // path(vcf), 
                            // path(bams), path(bai), val(dummy)
    FilterBcfVcf(fbcf_ch )            
    //----------------------Majority filter----------------------------------- 
    input_ch.map{row -> tuple(row[0],row[4])}
            .unique()
            .combine(parent_ch.map{row->tuple(row[2],row[0])},by:1)
            .map{row->tuple(row[1],row[2])}
            .join(sfilter_ch.vcf)
            .set{mjf_ch} // Emits val(groupId),val(parentlist), path(vcf), val(dummy)
    MajorityFilter(mjf_ch)
    //----------------------Plot-----------------------------------------
    // Input Channel Emits val(groupId),val(parentId),val(refpath),val(prefix),path(rds)
    input_ch.map{row -> tuple(row[0],row[4])}
            .unique()
            .join(copynum_ch.groupTuple().map{row-> tuple(row[0], row[1].flatten())},by:0).set{plot_ch}
    RPlotFull(plot_ch)
    RPlotROI(plot_ch)
    //-------------------------------------------------------------------
}




