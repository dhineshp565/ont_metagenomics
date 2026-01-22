#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include local modules
include { make_csv } from './modules/local/make_csv.nf'
include { merge_fastq } from './modules/local/merge_fastq.nf'
include { nanoplot } from './modules/local/nanoplot.nf'
include { porechop } from './modules/local/porechop.nf'
include { minimap2 } from './modules/local/minimap2.nf'
include { kraken2 } from './modules/local/kraken2.nf'
include { krona_kraken } from './modules/local/krona_kraken.nf'
include { make_report } from './modules/local/make_report.nf'
include { bracken } from './modules/local/bracken.nf'

// Include dehost subworkflow
include { DEHOST } from './subworkflows/dehost.nf'
include { METAGENOMICS } from './subworkflows/metagenomics.nf'


workflow {
	
	QCREADS(params.input, params.qscore, params.trim_barcodes)
		
		// Dehost after porechop (if enabled)
	if (params.dehost) {
			host_ref = file(params.host_db)
			DEHOST(QCREADS.out.reads, host_ref)
			reads_for_classification= DEHOST.out.dehosted_reads
			stats=DEHOST.out.stats
		} else {
			reads_for_classification = QCREADS.out.reads
	}
		
	
	METAGENOMICS (reads_for_classification,params.kraken_db,params.blastdb_path,params.blastdb_name)
	//generate report
	rmd_file=file("${baseDir}/ont_metagenomics.Rmd")
	meta_bracken=METAGENOMICS.out.bracken_output.map{ sample, file -> file }.collect()
	meta_blast= METAGENOMICS.out.blast_best.map {sample, file -> file}.collect()
	make_report(QCREADS.out.csv,METAGENOMICS.out.krona_output,rmd_file,meta_bracken,meta_blast)
	
	
	
}
