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


workflow {
    //----------------Input Preparation-----------------------------------------
    Channel.fromPath(params.input_file,checkIfExists:true)
    .splitCsv(header:true,sep:'\t')
    .map{row ->
        ref=""
        refpath=""
        if (row.ref =="3D7") {
            ref=params.ref3D7_path+"/PlasmoDB-52_Pfalciparum3D7_Genome"
            refpath=params.ref3D7_path
        }        
        return tuple(row.groupId,row.sampleId,row.fastqbase,row.ref,row.parentId,ref,refpath)
    }
    .set{input_ch} // Emits 0->groupId,	1->sampleId, 2->fastqbase,	
                   // 3->ref_prefix, 4->parentId, 5->ref (path+name)
                   // 6->refpath 

    //----------------Alignment-----------------------------------------
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

}




