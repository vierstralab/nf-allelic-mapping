#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include  { waspRealigning; add_snp_files_to_meta; get_container } from './main'

/* 
   Variation of the pipeline that performs pre-merge of aggregations by indiv_id and cell_type.
   Combination of indiv_id and cell_type is used in place of the aggregation id during waspRealigning workflow.
*/

// TODO: ADD TEST CASE
process merge_by_indiv_and_cell_type {
	tag "${indiv_id}:${cell_type}"
	publishDir params.outdir + '/merged'
	container "${params.container}"
	containerOptions "${get_container(params.genome_fasta_file)}"
	cpus 2

	input:
		tuple val(indiv_id), val(cell_type), val(bam_files)

	output:
		tuple val(ag_placeholder), val(indiv_id), path(name), path("${name}.crai")

	script:
	ag_placeholder="${indiv_id}_${cell_type}"
	name = "${ag_placeholder}.cram"
	"""
	samtools merge -f -@${task.cpus} \
		-O CRAM --write-index \
		--reference ${params.genome_fasta_file} \
		${name} ${bam_files}
	"""


}

workflow {
	bams_grouped_by_indiv_and_cell_type = Channel.fromPath(params.samples_file)
		| splitCsv(header:true, sep:'\t')
		| map(row -> tuple(row.indiv_id, row.cell_type, file(row.bam_file)))
		| groupTuple(by: [0,1])
		| map(it -> tuple(it[0], it[1], it[2].flatten().join(" ")))
		| merge_by_indiv_and_cell_type
		| waspRealigning

	add_snp_files_to_meta() 
}
