#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.publishdir = "$launchDir"
params.outdir = ""
params.biomix_root = "$moduleDir/bin/NextflowModules/bin/BiomiX2.5"
params.command_dir = ""
params.transcriptomics_matrix = ""
params.methylomics_matrix = ""
params.metadata = ""
params.group_1 = ""
params.group_2 = ""
params.transcriptomics_label = "RNA"
params.methylomics_label = "METHY"
params.run_transcriptomics = true
params.run_methylomics = false
params.run_mofa = false
params.compare_gold = false
params.gold_standard_dir = ""
params.gold_manifest = ""

process PREPARE_BIOMIX_WORKSPACE {
    publishDir "${params.outdir ?: params.publishdir}/nf_output", mode: 'copy', pattern: 'biomix_workspace/COMBINED_COMMANDS.json'

    conda "${moduleDir}/bin/conda_biomix_transcriptomics.yml"

    input:
    path biomix_root
    path command_dir
    path transcriptomics_matrix
    val methylomics_matrix
    path metadata
    val transcriptomics_label
    val methylomics_label
    val group_1
    val group_2

    output:
    path "biomix_workspace", emit: workspace

    script:
    """
    python ${moduleDir}/bin/biomix_prepare_workspace.py \
      --source-root "$biomix_root" \
      --dest-root biomix_workspace \
      --command-dir "$command_dir" \
      --transcriptomics-matrix "$transcriptomics_matrix" \
      --methylomics-matrix "${methylomics_matrix}" \
      --metadata "$metadata" \
      --transcriptomics-label "${transcriptomics_label}" \
      --methylomics-label "${methylomics_label}" \
      --output-dir "\$PWD/biomix_workspace"

    (
      cd biomix_workspace
      Rscript Converter_JSON.r "${group_1}" "${group_2}" "\$PWD"
    )
    """
}

process RUN_TRANSCRIPTOMICS {
    publishDir "${params.outdir ?: params.publishdir}/nf_output", mode: 'copy'

    conda "${moduleDir}/bin/conda_biomix_transcriptomics.yml"

    input:
    path biomix_workspace
    val transcriptomics_label
    val group_1
    val group_2

    output:
    path "biomix_workspace_post_transcriptomics", emit: workspace
    path "Transcriptomics", emit: transcriptomics
    path "Integration", emit: integration

    script:
    """
    cp -r "$biomix_workspace" biomix_workspace_post_transcriptomics

    Rscript ${moduleDir}/bin/biomix_run_transcriptomics.R \
      --workspace "\$PWD/biomix_workspace_post_transcriptomics" \
      --group1 "${group_1}" \
      --group2 "${group_2}" \
      --label "${transcriptomics_label}"

    cp -r biomix_workspace_post_transcriptomics/Transcriptomics ./Transcriptomics
    cp -r biomix_workspace_post_transcriptomics/Integration ./Integration
    """
}

process RUN_METHYLOMICS {
    publishDir "${params.outdir ?: params.publishdir}/nf_output", mode: 'copy'

    conda "${moduleDir}/bin/conda_biomix_methylomics.yml"

    input:
    path biomix_workspace
    val methylomics_label
    val group_1
    val group_2

    output:
    path "biomix_workspace_post_methylomics", emit: workspace
    path "Methylomics", emit: methylomics
    path "Integration", emit: integration

    script:
    """
    cp -r "$biomix_workspace" biomix_workspace_post_methylomics

    Rscript ${moduleDir}/bin/biomix_run_methylomics.R \
      --workspace "\$PWD/biomix_workspace_post_methylomics" \
      --group1 "${group_1}" \
      --group2 "${group_2}" \
      --label "${methylomics_label}"

    cp -r biomix_workspace_post_methylomics/Methylomics ./Methylomics
    cp -r biomix_workspace_post_methylomics/Integration ./Integration
    """
}

process RUN_MOFA {
    publishDir "${params.outdir ?: params.publishdir}/nf_output", mode: 'copy'

    conda "${moduleDir}/bin/conda_biomix_mofa.yml"

    input:
    path biomix_workspace
    val group_1
    val group_2

    output:
    path "biomix_workspace_post_mofa", emit: workspace
    path "Integration", emit: integration

    script:
    """
    cp -r "$biomix_workspace" biomix_workspace_post_mofa

    Rscript ${moduleDir}/bin/biomix_run_mofa.R \
      --workspace "\$PWD/biomix_workspace_post_mofa" \
      --group1 "${group_1}" \
      --group2 "${group_2}"

    cp -r biomix_workspace_post_mofa/Integration ./Integration
    """
}

process COMPARE_GOLD_STANDARD {
    conda "${moduleDir}/bin/conda_biomix_transcriptomics.yml"

    input:
    path actual_root
    path gold_root
    path manifest

    output:
    path "gold_comparison_report.json", emit: report

    script:
    """
    python ${moduleDir}/bin/compare_biomix_gold.py \
      --actual-root "$actual_root" \
      --gold-root "$gold_root" \
      --manifest "$manifest" \
      --report gold_comparison_report.json
    """
}

workflow Main {
    take:
    input_map

    main:
    prepared_workspace = PREPARE_BIOMIX_WORKSPACE(
        input_map.biomix_root,
        input_map.command_dir,
        input_map.transcriptomics_matrix,
        input_map.methylomics_matrix,
        input_map.metadata,
        input_map.transcriptomics_label,
        input_map.methylomics_label,
        input_map.group_1,
        input_map.group_2
    )

    active_workspace = prepared_workspace.workspace
    transcriptomics_outputs = null
    methylomics_outputs = null

    if (input_map.run_transcriptomics) {
        transcriptomics_outputs = RUN_TRANSCRIPTOMICS(
            active_workspace,
            input_map.transcriptomics_label,
            input_map.group_1,
            input_map.group_2
        )
        active_workspace = transcriptomics_outputs.workspace
    }

    if (input_map.run_methylomics) {
        methylomics_outputs = RUN_METHYLOMICS(
            active_workspace,
            input_map.methylomics_label,
            input_map.group_1,
            input_map.group_2
        )
        active_workspace = methylomics_outputs.workspace
    }

    if (input_map.run_mofa) {
        mofa_outputs = RUN_MOFA(
            active_workspace,
            input_map.group_1,
            input_map.group_2
        )
        active_workspace = mofa_outputs.workspace
    }

    if (input_map.compare_gold) {
        COMPARE_GOLD_STANDARD(
            active_workspace,
            input_map.gold_standard_dir,
            input_map.gold_manifest
        )
    }

    emit:
    workspace = active_workspace
    transcriptomics = transcriptomics_outputs ? transcriptomics_outputs.transcriptomics : Channel.empty()
    methylomics = methylomics_outputs ? methylomics_outputs.methylomics : Channel.empty()
}

workflow {
    [
        biomix_root: params.biomix_root,
        command_dir: params.command_dir,
        transcriptomics_matrix: params.transcriptomics_matrix,
        metadata: params.metadata,
        group_1: params.group_1,
        group_2: params.group_2
    ].each { key, value ->
        if (!value?.toString()?.trim()) {
            error "Missing required parameter: --${key}"
        }
    }

    if (params.compare_gold) {
        if (!params.gold_standard_dir?.toString()?.trim()) {
            error "Missing required parameter: --gold_standard_dir"
        }
        if (!params.gold_manifest?.toString()?.trim()) {
            error "Missing required parameter: --gold_manifest"
        }
    }

    input_map = [
        biomix_root: file(params.biomix_root, checkIfExists: true),
        command_dir: file(params.command_dir, checkIfExists: true),
        transcriptomics_matrix: file(params.transcriptomics_matrix, checkIfExists: true),
        methylomics_matrix: params.methylomics_matrix ? file(params.methylomics_matrix, checkIfExists: true) : "",
        metadata: file(params.metadata, checkIfExists: true),
        group_1: params.group_1,
        group_2: params.group_2,
        transcriptomics_label: params.transcriptomics_label,
        methylomics_label: params.methylomics_label,
        run_mofa: params.run_mofa,
        run_transcriptomics: params.run_transcriptomics,
        run_methylomics: params.run_methylomics,
        compare_gold: params.compare_gold,
        gold_standard_dir: params.compare_gold ? file(params.gold_standard_dir, checkIfExists: true) : null,
        gold_manifest: params.compare_gold ? file(params.gold_manifest, checkIfExists: true) : null
    ]

    Main(input_map)
}
