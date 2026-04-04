.PHONY: test-fastq test-cram test-cram-multi test-cram-multisv test-e2e clean lint help pixi test-cram-sv-manta test-cram-sv-delly test-cram-sv-tiddit test-cram-sv-lumpy test-cram-sv-cnvnator test-cram-sv-manta-delly test-cram-sv-manta-lumpy test-cram-sv-all

${HOME}/.pixi/bin/pixi:
	curl -sSL https://pixi.sh/install.sh | sh


# FASTQ input test - full pipeline with alignment
test-fastq: ${HOME}/.pixi/bin/pixi
	${HOME}/.pixi/bin/pixi run nextflow run main.nf \
		--input assets/samplesheet_fastq.csv \
		-profile docker \
		-resume

# Lint
lint: ${HOME}/.pixi/bin/pixi
	${HOME}/.pixi/bin/pixi run nextflow lint $(shell find . -name "*.nf") -format

# Clean
clean:
	rm -rf work results .nextflow* *.log