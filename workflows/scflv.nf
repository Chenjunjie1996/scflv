/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                } from '../modules/local/multiqc_sgr/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_scflv_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FASTQC {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fastqc=0.12.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads), path(match_dir)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path  "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Make list of old name and new name pairs to use for renaming in the bash while loop
    def old_new_pairs = reads instanceof Path || reads.size() == 1 ? [[ reads, "${prefix}.${reads}" ]] : reads.withIndex().collect { entry, index -> [ entry, "${prefix}.${entry}" ] }
    def rename_to = old_new_pairs*.join(' ').join(' ')
    def renamed_files = old_new_pairs.collect{ old_name, new_name -> new_name }.join(' ')
    """
    printf "%s %s\\n" $rename_to | while read old_name new_name; do
        [ -f "\${new_name}" ] || ln -s \$old_name \$new_name
    done

    fastqc \\
        $args \\
        --threads $task.cpus \\
        $renamed_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}

process PROTOCOL_CMD  {
    tag "$meta.id"
    label 'process_single'

    conda 'bioconda::pyfastx=2.1.0'
    container "biocontainers/pyfastx:2.1.0--py39h3d4b85c_0"

    conda 'conda-forge::pandas==1.5.2'
    container "biocontainers/pandas:1.5.2"

    conda 'conda-forge::biopython==1.82'
    container "biocontainers/biopython:1.82"

    input:
    tuple val(meta), path(reads), path(match_dir)
    path assets_dir
    val protocol

    output:
    tuple val(meta), path("${meta.id}_R*.fq*"),  emit: out_reads
    path "${meta.id}.protocol_stats.json", emit: protocol_stats_json
    path  "versions.yml" , emit: versions

    script:
    // separate forward from reverse pairs
    def (r1,r2) = reads.collate(2).transpose()
    """
    protocol_cmd.py \\
        --sample ${meta.id} \\
        --fq1 ${r1.join( "," )} \\
        --fq2 ${r2.join( "," )} \\
        --match_dir ${match_dir} \\
        --assets_dir ${assets_dir} \\
        --protocol ${protocol} 


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyfastx: \$(pyfastx --version | sed -e "s/pyfastx version //g")
        pandas: \$(pandas --version | sed -e "s/pandas version //g")
        biopython: \$(biopython --version | sed -e "s/biopython version //g")
    END_VERSIONS
    """
}

process RUN_TRUST4 {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::trust4=1.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trust4:1.1.1--h43eeafb_0':
        'biocontainers/trust4:1.1.1--h43eeafb_0' }"

    conda "conda-forge::gawk"
    
    input:
    tuple val(meta), path(reads)
    val (species)

    output:
    tuple val(meta), path("*_toassemble*")         , emit: candidate_reads
    tuple val(meta), path("*_report.tsv")          , emit: report_tsv
    tuple val(meta), path("*_filter_report.tsv")   , emit: filter_report_tsv
    tuple val(meta), path("*_raw.out")             , emit: raw_out
    tuple val(meta), path("*_final.out")           , emit: final_out
    tuple val(meta), path("*_cdr3.out")            , emit: cdr3_out
    tuple val(meta), path("*_assign.out")          , emit: assign_out
    tuple val(meta), path("*_barcode_report.tsv")  , emit: barcode_report
    tuple val(meta), path("*_barcode_airr.tsv")    , emit: barcode_airr
    tuple val(meta), path("*_b.csv")               , emit: barcode_report_b
    tuple val(meta), path("*_t.csv")               , emit: barcode_report_t
    tuple val(meta), path("*_assembled_reads.fa")  , emit: assembled_reads
    tuple val(meta), path("*_annot.fa")            , emit: annot_fa
    tuple val(meta), path("*_airr.tsv")            , emit: airr_tsv
    tuple val(meta), path("*_airr_align.tsv")      , emit: airr_alin
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def (r1, r2) = reads.collate(2).transpose()
    def readformat = task.ext.readformat ?: "bc:0:23,um:24:-1"
    def run_trust4_cmd = "-u ${r2[0]} --barcode ${r1[0]} --UMI ${r1[0]} --outputReadAssignment"
    def fasta = null
    def vdj_reference = null
    if( species == 'human' ) {
        fasta = "${projectDir}/assets/human/hg38_bcrtcr.fa"
        vdj_reference = "${projectDir}/assets/human/human_IMGT+C.fa"
    } else {
        fasta = "${projectDir}/assets/mouse/GRCm38_bcrtcr.fa.fa"
        vdj_reference = "${projectDir}/assets/mouse/mouse_IMGT+C.fa"
    }
    """
    echo $vdj_reference
    run-trust4 \\
        ${run_trust4_cmd} \\
        -t $task.cpus \\
        -f ${fasta} \\
        -o ${prefix} \\
        --ref ${vdj_reference} \\
        --readFormat ${readformat} \\
        $args

    gawk '$4!~"_" && $4!~"?"' ${prefix}_report.tsv > ${prefix}_filter_report.tsv
    
    perl trust-barcoderep-to-10X.pl ${prefix}_barcode_report.tsv ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trust4: \$(run-trust4 2>&1 | grep -o 'v[0-9.]*-r[0-9]*' | sed 's/^/TRUST4 using /' )
    END_VERSIONS
    """
}

process SUMMARIZE {
    tag "$meta.id"
    label 'process_single'

    conda 'conda-forge::pandas==1.5.2'
    container "biocontainers/pandas:1.5.2"

    conda 'bioconda::pyfastx=2.1.0'
    container "biocontainers/pyfastx:2.1.0--py39h3d4b85c_0"

    conda 'conda-forge::numpy==1.24.4'
    container "biocontainers/numpy:1.24.4"

    input:
    tuple val(meta), path(reads), path(match_dir)
    val seqtype
    val coef
    val expected_target_cell_num
    val target_cell_barcode
    val target_weight
    path fq2
    path assembled_reads
    path filter_report_tsv
    path annot_fa
    path barcode_report

    output:
    path "${meta.id}.count.txt", emit: umi_count_txt
    path "${meta.id}.cells_stats.json", emit: summary_json
    path "${meta.id}.umi_count.json", emit: umi_count_json

    script:

    """
    summarize.py \\
        --sample ${meta.id} \\
        --seqtype ${seqtype} \\
        --coef ${coef} \\
        --expected_target_cell_num ${expected_target_cell_num} \\
        --target_cell_barcode ${target_cell_barcode} \\
        --target_weight ${target_weight} \\
        --fq2 ${fq2} \\
        --assembled_reads ${assembled_reads} \\
        --filter_report_tsv ${filter_report_tsv} \\
        --annot_fa ${annot_fa} \\
        --barcode_report ${barcode_report} \\
    """
}

// process ANNOTATE {

// }

workflow scflv {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // MODULE: Run FastQC
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // protocol
    PROTOCOL_CMD (
        ch_samplesheet,
        "${projectDir}/assets/",
        params.protocol,
    )
    ch_versions = ch_versions.mix(PROTOCOL_CMD.out.versions.first())

    // trust4
    RUN_TRUST4 (
        PROTOCOL_CMD.out.out_reads,
        params.species,
    )

    // ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.results.collect{it[1]})
    ch_versions = ch_versions.mix(RUN_TRUST4.out.versions.first())
    
    // SUMMARIZE

    // ANNOTATE


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
