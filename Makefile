run:
	nextflow run ./nf_workflow.nf -resume -c nextflow.config \
		--biomix_root ./bin/BiomiX2.5 \
		--command_dir ./test/fixtures/egas_transcriptomics_mutated_vs_unmutated \
		--transcriptomics_matrix ./bin/BiomiX2.5/Example_dataset/EGAS00001001746/RNA_seq/EGAS00001001746_transcriptomics.tsv \
		--metadata ./bin/BiomiX2.5/Example_dataset/EGAS00001001746/Metadata/EGAS00001001746_metadata_CLL.tsv \
		--group_1 mutated \
		--group_2 unmutated

run_importer_workflow:
	nextflow run ./nf_workflow_importer.nf -resume -c nextflow.config

run_slurm:
	nextflow run ./nf_workflow.nf -resume -c nextflow_slurm.config \
		--biomix_root ./bin/NextflowModules/bin/BiomiX2.5 \
		--command_dir ./test/fixtures/egas_transcriptomics_mutated_vs_unmutated \
		--transcriptomics_matrix ./bin/BiomiX2.5/Example_dataset/EGAS00001001746/RNA_seq/EGAS00001001746_transcriptomics.tsv \
		--metadata ./bin/BiomiX2.5/Example_dataset/EGAS00001001746/Metadata/EGAS00001001746_metadata_CLL.tsv \
		--group_1 mutated \
		--group_2 unmutated

run_docker:
	nextflow run ./nf_workflow.nf -resume -with-docker <CONTAINER NAME>

init_modules:
	git submodule update --init --recursive
