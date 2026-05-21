#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include local modules
include { make_csv } from './modules/local/make_csv.nf'
include { merge_filter_fastq } from './modules/local/merge_filter_fastq.nf'
include { porechop } from './modules/local/porechop.nf'
include { minimap2 } from './modules/local/minimap2.nf'
include { kraken2 } from './modules/local/kraken2.nf'
include { krona_kraken } from './modules/local/krona_kraken.nf'
include { make_report } from './modules/local/make_report.nf'
include { make_metagenomics_report } from './modules/local/make_metagenomics_report.nf'
include { bracken } from './modules/local/bracken.nf'
include { multiqc } from './modules/local/multiqc.nf'

// Include dehost subworkflow
include { DEHOST } from './subworkflows/dehost.nf'
include { METAGENOMICS } from './subworkflows/metagenomics.nf'
include { QCREADS } from './subworkflows/qcreads.nf'


workflow {
	
	QCREADS(params.input, params.qscore, params.trim_barcodes)
		
		// Dehost after porechop (if enabled)
	if (params.dehost) {
			host_ref = file(params.host_db)
			DEHOST(QCREADS.out.reads, host_ref)
			reads_for_classification= DEHOST.out.dehosted_reads
			dehost_stats = DEHOST.out.stats
		} else {
			reads_for_classification = QCREADS.out.reads
			dehost_stats = Channel.empty()
	}

	// Reuse NanoPlot stats from QCREADS and add optional dehost stats.
	multiqc_inputs = QCREADS.out.read_stats.mix(dehost_stats).collect()
	multiqc(multiqc_inputs)
	
	
	METAGENOMICS (reads_for_classification,params.kraken_db,params.blastdb_path,params.blastdb_name,params.kraken_confidence)
	
	// Generate metagenomics report
	rmd_meta_file=file("${baseDir}/ont_metagenomics.Rmd")
	meta_summary=METAGENOMICS.out.merged_report.map{ _sample, _blast, summary -> summary }.collect()
	blast_hits=METAGENOMICS.out.merged_report.map{ _sample, blast, _summary -> blast }.collect()
	make_metagenomics_report(QCREADS.out.csv,METAGENOMICS.out.krona_output,rmd_meta_file,meta_summary,blast_hits)
	
	
	
}
