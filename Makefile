.PHONY: test-fastq test-cram test-cram-multi test-cram-multisv test-e2e clean lint help pixi test-cram-sv-manta test-cram-sv-delly test-cram-sv-tiddit test-cram-sv-lumpy test-cram-sv-cnvnator test-cram-sv-manta-delly test-cram-sv-manta-lumpy test-cram-sv-all

${HOME}/.pixi/bin/pixi:
	curl -sSL https://pixi.sh/install.sh | sh

# snapshot
# nf-test snapshot tests
test-fastq-snapshot: ${HOME}/.pixi/bin/pixi
	export NXF_FILE_ROOT=${PWD}; ${HOME}/.pixi/bin/pixi run nf-test test \
		--verbose \
		--profile docker,test_fastq

# Update nf-test snapshots
test-fastq-update-snapshot: ${HOME}/.pixi/bin/pixi
	export NXF_FILE_ROOT=${PWD}; ${HOME}/.pixi/bin/pixi run nf-test test \
		tests/default.nf.test \
		--verbose \
		--update-snapshot \
		--profile docker,test_fastq

# FASTQ input test - full pipeline with alignment
test-fastq: ${HOME}/.pixi/bin/pixi
	${HOME}/.pixi/bin/pixi run nextflow run main.nf \\
		-profile docker,test_fastq \
		-resume

test-fastq-multi: ${HOME}/.pixi/bin/pixi
	${HOME}/.pixi/bin/pixi run nextflow run main.nf \
		-profile docker,test_fastq \
		-resume

# CRAM input test - single sample
test-cram: ${HOME}/.pixi/bin/pixi
	${HOME}/.pixi/bin/pixi run nextflow run main.nf \
		-profile docker,test_cram \
		--index_bwa2_reference false \
		-resume

# CRAM input test - multiple samples
test-cram-multi: ${HOME}/.pixi/bin/pixi
	${HOME}/.pixi/bin/pixi run nextflow run main.nf \
		-profile docker,test_cram_multi \
		--index_bwa2_reference false \
		-resume
# Lint
lint: ${HOME}/.pixi/bin/pixi
	${HOME}/.pixi/bin/pixi run nextflow lint $(shell find . -name "*.nf") -format

# Clean
clean:
	rm -rf work results .nextflow* *.log