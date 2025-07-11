// include modules

include {
        LaneMerge ;
        Bwa ;
        WriteBamLists ;
        Index ;
        Merge
} from './modules/bwa.nf'

include {
        Bcf ;
        Gridss ;
        SomaticFilter ;
        RCopyNum ;
        FilterBcfVcf ;
        FilterGridssSV ;
        RPlotFull ;
        RPlotROI
} from './modules/stvariant.nf'


include {
        FastQC ;
        MosDepth ;
        FlagStats ;
        MultiQC
} from './modules/qc.nf'

include { validateParameters ; paramsSummaryLog ; samplesheetToList } from 'plugin/nf-schema'

def validateSampleIdContainsGroupId(List row) {
        if (!row[1].toString().contains(row[0].toString())) {
                error("Validation failed: sampleId '${row[1]}' does not contain groupId '${row[0]}'")
        }
        return row
}
workflow {
        println("*****************************************************")
        println("*  Nextflow Malaria Variant analysis pipeline       *")
        println("*  A Nextflow wrapper pipeline                      *")
        println("*  Written by WEHI Research Computing Platform      *")
        println("*  research.computing@wehi.edu.au                   *")
        println("*                                                   *")
        println("*****************************************************")
        println(" Required Pipeline parameters                        ")
        println("-----------------------------------------------------")
        println("Input Sequence Path: ${params.input_seq_path}          ")
        println("Sample Group file  : ${params.input_file}              ")
        println("Output directory   : ${params.outdir}                  ")
        println("*****************************************************")

        validateParameters()
        // Create a new channel of metadata from a sample sheet passed to the pipeline through the --input parameter
        Channel.fromList(samplesheetToList(params.input_file, "input_schema.json"))
                .map { validateSampleIdContainsGroupId(it) }

        //----------------Input Preparation-----------------------------------------
        Channel.fromPath(params.input_file, checkIfExists: true)
                .ifEmpty {
                        error(
                                """
        No groups could be found in group key file! Please check your directory path
        is correct. 
        """
                        )
                }
                .splitCsv(header: true, sep: '\t')
                .map { row ->
                        def ref = ""
                        def refpath = ""
                        def bsref = ""
                        if (row.ref == "3D7") {
                                ref = params.ref3D7_path + "/PlasmoDB-52_Pfalciparum3D7_Genome"
                                refpath = params.ref3D7_path
                                bsref = "BSgenome.Pfalciparum3D7.PlasmoDB.52"
                        }
                        else if (row.ref == "Dd2") {
                                ref = params.refDd2_path + "/PlasmoDB-57_PfalciparumDd2_Genome"
                                refpath = params.refDd2_path
                                bsref = "BSgenome.PfalciparumDd2.PlasmoDB.57"
                        }
                        else if (row.ref == "Supp") {
                                ref = params.refSupp_path + "/PlasmoDB-52_Pfalciparum3D7_Genome_supplemented"
                                refpath = params.refSupp_path
                                bsref = "BSgenome.PfalciparumNF54iGP"
                        }
                        return tuple(row.groupId, row.sampleId, row.fastqbase, row.ref, row.parentId, ref, refpath, bsref)
                }
                .set { input_ch }
                // Emits 0->groupId,	1->sampleId, 2->fastqbase,	
                // 3->ref_prefix, 4->parentId, 5->ref (path+name)
                // 6->refpath , 7-> bsref

        //----------------Alignment-----------------------------------------
        bamlist_ch = WriteBamLists(Channel.fromPath(params.input_file, checkIfExists: true))
        bamlist_ch = bamlist_ch
                .flatten()
                .map { row ->
                        tuple(row.baseName.split("_bams")[0], row)
                }
        // Emits groupID, bamslist.txt

        if (params.merge_lanes) {
                files_ch = input_ch
                        .map { row ->
                                def pattern = "${params.input_seq_path}/${row[1]}*_R1*"
                                def files = file(pattern)

                                // Handle case where no files match
                                if (files.isEmpty()) {
                                        log.warn("No files found matching pattern: ${pattern}")
                                        return []
                                }

                                return tuple(row[2], files)
                        }
                        .map { row ->
                                def pattern = "${params.input_seq_path}/${row[0]}*_R2*"
                                def files = file(pattern)

                                // Handle case where no files match
                                if (files.isEmpty()) {
                                        log.warn("No files found matching pattern: ${pattern}")
                                        return []
                                }

                                return tuple(row[0], row[1], files)
                        }
                // Emits fastqbase, R1 files, R2 files

                LaneMerge(files_ch).set { mergedfastqfiles }
                // Emits 0->fastqbase, 2->[fastqs]
                input_ch
                        .map { row -> tuple(row[2], row[1], row[0], row[5], row[6]) }
                        // Emits 0->sampleId, 2->[fastqs]
                        .join(mergedfastqfiles, by: 0)
                        .set { bwa_input_ch }
                        // Emits tuple val(fastqbase), val(sampleId),val(groupId),val(ref),path(refpath),path(fastqs)
        }
        else {
                Channel.fromFilePairs("${params.input_seq_path}/*_{,R}{1,2}*.{fq,fastq}{,.gz}", size: 2)
                        .ifEmpty {
                                error(
                                        """
                    No samples could be found! Please check whether your input directory
                    is correct, and that your samples match typical fastq paired end naming
                    convention(s).
                    """
                                )
                        }
                        .map { row -> tuple(row[0].split(/_R[0-9]{1}/)[0], row[1]) }
                        .set { fastq_input_channel }

                input_ch
                        .map { row -> tuple(row[2], row[1], row[0], row[5], row[6]) }
                        .join(fastq_input_channel, by: 0)
                        .set { bwa_input_ch }
                        // Emits tuple val(fastqbase), val(sampleId),val(groupId),val(ref),path(refpath),path(fastqs)
        }
        sam_ch = Bwa(bwa_input_ch)
        bam_ch = Index(sam_ch)
        //----------------Merge&List---------------------------------------
        input_ch
                .map { row -> row[4] }
                .unique()
                .combine(bam_ch.bamnodup, by: 0)
                .map { row ->
                        tuple(row[0], row[1], row[1].baseName + ".bam")
                }
                .groupTuple()
                .map { row ->
                        tuple(row[0], row[1], row[2].join(" "))
                }
                .ifEmpty {
                        error(
                                """
                    No parent samples found.
                    """
                        )
                }
                .set { parent_ch }
        // Emits tuple val(parentId),path(bams), val(parentlist)

        input_ch
                .map { row -> row[4] }
                .unique()
                .combine(bam_ch.bai, by: 0)
                .groupTuple()
                .join(parent_ch)
                .ifEmpty {
                        error(
                                """
                    No parent samples found.
                    """)
            }
            .set{parent_all_ch} // Emits tuple val(parentId),path(bai),path(bams), val(parentlist)
    
    merged_ch=Merge(parent_ch) // Emits tuple val(parentId),path(parentId.bam)
    //----------------QC tools------------------------------------------
    //MultiQC integrates files from fastqc and mosdepth, found in input dir, into single report
    fastqc_ch=FastQC(bam_ch.bamnodup.map{row->row[1]}.unique({it.baseName}).collect())
    
    // MosDepth Input Channel include val(sampleId), path(bam), path(bai)
    mosdepth_ch=MosDepth(bam_ch.bysampleid) // join bams with their index based on groupId 
    // FlagStats Input Channel include val(sampleId), path(bam)
    flagstat_ch=FlagStats(bam_ch.bysampleid)  
    
    MultiQC( fastqc_ch.zip.mix(mosdepth_ch, flagstat_ch).collect().ifEmpty([]),Channel.fromPath(params.multiqc_config) )  
    //----------------BCF tools----------------------------------------
    //BCF Input Channel Emits parentId,groupId,ref,bamlist,bams,parentbams
    input_ch.map{row -> tuple(row[0],row[4],row[5]+".fasta")}
            .unique()
            .join(bamlist_ch)
            .combine(bam_ch.bamnodup,by:0)
            .groupTuple(by:[0,1,2,3])
            .combine(parent_ch.map{row->tuple(row[1],row[0])},by:1).set{bcf_input_ch}
               // Emits val(parentId),val(groupId),path(ref), path(bamlist), path(bams),path(parentbams)
    bcf_ch=Bcf(bcf_input_ch)
    //----------------------Gridss------------------------------------- 
    input_ch.map{row -> tuple(row[0],row[4],row[5]+".fasta")}
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
            .set{gridss_input_ch}// Emits val(parentId),val(groupId),path(ref), val(bamlistcontent), path(bams),path(parentbams)
    gridss_ch=Gridss(gridss_input_ch.combine(Channel.fromPath(params.gridss_jar_path)))
    
    input_ch.map{row -> tuple(row[4], row[0],row[7] )}
            .unique()
            .combine(parent_ch,by:0).map{row -> tuple(row[1],row[2],row[5])}
            .join(bamlist_ch).join(gridss_ch.vcf)
            .set{sv_input_ch} //Emit val(groupId),val(bsref),val(parentbamlist), path(bamlist), path(vcf)

    sfilter_ch=SomaticFilter(sv_input_ch
                                .combine(Channel.fromPath("${projectDir}/Rtools/gridss_assets/gridss_somatic_filter.R")))
    //----------------------CopyNum--------------------------------------
    input_ch.map{row -> tuple(row[4], row[0],row[6],row[7] )}
            .unique()
            .combine(merged_ch,by:0)
            .map{row -> if (row[0]) tuple(row[1],row[0],row[2],row[3],row[4])}
            .combine(bam_ch.bamnodup,by:0)
            .groupTuple(by:[0,1,2,3,4])
            .ifEmpty {
                error("""
                Input to CopyNum Analysis is empty.
                """)
            }
            .set{copynum_input_ch}
              //Emits val(groupId),val(parentId),path(refpath),val(bsref),path(mergedparent), path(bamlistcontent), path(bams), val(dummy)
    
    copynum_ch=RCopyNum(copynum_input_ch
                            .combine(Channel.fromList( [params.bin_CNroi, params.bin_CNfull] ))
                            .combine(Channel.fromPath("${projectDir}/Rtools/malDrugR/copynumQDNAseq.R"))
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
            .set{fbcf_ch}   // Emits val(groupId),path(refpath),val(prefix),
                            // path(parentbai),path(parentbam), val(parentbamlist), 
                            // path(vcf), 
                            // path(bams), path(bai), val(dummy)
    FilterBcfVcf(fbcf_ch
                        .combine(Channel.fromPath("${projectDir}/Rtools/malDrugR/filterBcfVcf.R"))     
    )       
    //----------------------GRIDSS filter----------------------------------- 
    input_ch.map{row -> tuple(row[0],row[4])}
            .unique()
            .combine(parent_ch.map{row->tuple(row[2],row[0])},by:1)
            .map{row->tuple(row[1],row[2])}
            .join(sfilter_ch.vcf)
            .set{gf_ch} // Emits val(groupId),val(parentlist), path(vcf), val(dummy)
    FilterGridssSV(gf_ch
                        .combine(Channel.fromPath("${projectDir}/Rtools/malDrugR/gridss_majorityfilt.R"))
                        .combine(Channel.fromPath("${projectDir}/Rtools/gridss_assets/"))
        )
    //----------------------Plot CopyNums-----------------------------------------
    // Input Channel Emits val(groupId),val(parentId),path(rds)
    input_ch.map{row -> tuple(row[0],row[4])}
            .unique()
            .join(copynum_ch.groupTuple().map{row-> tuple(row[0], row[1].flatten())},by:0).set{plot_ch}

    RPlotFull(plot_ch.combine(Channel.fromPath("${projectDir}/Rtools/malDrugR/copynumPlotsFull.R")))
    RPlotROI(plot_ch.combine(Channel.fromPath("${projectDir}/Rtools/malDrugR/copynumPlotsROI.R")))
    
    //----------------------List genes in a region of interest, if specified------
    // Input Channel emits path(refpath), val(ref_prefix)
    input_ch.map{row -> tuple(row[6],row[3])}.unique().set{genesROI_input_ch}

    if (params.genesRegion != "") 
      genesROI(genesROI_input_ch.combine(Channel.fromPath("${projectDir}/Rtools/malDrugR/genesROI.R")))
    //-------------------------------------------------------------------
}