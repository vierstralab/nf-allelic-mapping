#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/* 
   Variation of the pipeline that performs pre-merge of aggregations by indiv_id and cell_type.
   Combination of indiv_id and cell_type is used in place of the aggregation id during waspRealigning workflow.
*/

include  { calcInitialReadCounts; waspRealigning; filter_grouped_channel; set_key_for_group_tuple; add_snp_files_to_meta } from './main'
process merge_by_indiv_and_cell_type {
        tag "${indiv_id}:${cell_type}"

        publishDir params.outdir + '/merged', mode: 'symlink' 

        module "samtools/1.14"

        cpus 4

        input:
        tuple val(indiv_id), val(cell_type), val(bam_files)

        output:
        tuple val(indiv_id), val(ag_placeholder), file('*.cram'), file('*.cram.crai')

        script:
	ag_placeholder="${indiv_id}_${cell_type}"
        """
        samtools merge -f -@${task.cpus} -O CRAM --write-index --reference ${params.genome_fasta_file} ${ag_placeholder}.cram ${bam_files}
        """


}

workflow {
        bams_grouped_by_indiv_and_cell_type = Channel
            .fromPath(params.samples_file)
            .splitCsv(header:true, sep:'\t')
            .map{ row -> tuple(row.indiv_id, row.cell_type, file(row.bam_file)) }
            .groupTuple(by: [0,1])
            .map{ it -> tuple(it[0], it[1], it[2].flatten().join(" ")) }
	
	samples_aggregations=merge_by_indiv_and_cell_type(bams_grouped_by_indiv_and_cell_type)

	waspRealigning(set_key_for_group_tuple(samples_aggregations))
	add_snp_files_to_meta() 
}
