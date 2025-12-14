#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ADD more aligners in the process below if needed. See example in the script body
process align_reads {
	container "${params.container}"
	// All external files (passed as params) should be wrapped 
	// in "get_container" function and added to containerOptions, see below
    // TODO: change to conventional tuple val(meta), val(other params).
	containerOptions "${get_container(params.genome_fasta_file)} ${get_container(params.nuclear_chroms)} ${params.aligner == 'bowtie-chip' ? get_container(params.bowtie_idx) : ''}"
	scratch true
	tag "${sample_id}:${r_tag}"
	cpus 3
	// fastq1 is the same as fastq2 if r_tag == "se"
	input:
		tuple val(sample_id), val(r_tag), path(fastq1), path(fastq2)
		// "sample_id" - ID of initial bam file that was splitted
		// into pe and se reads.

		// "r_tag" - indicates which reads were extracted from the bam file.
		// Can be either "se" or "pe"

		// fastq1, fastq2 reads from the bam file. 
		// Use any of them for se reads aligning

	output:
		tuple val(sample_id), val(r_tag), path(name)

	script:
	// Name of the output bam file
	name = "${sample_id}.${r_tag}.realigned.bam"
	switch(params.aligner) {
		case 'bwa-altius-dnase':
			// PE reads alignment
			template (r_tag == 'pe' ? 'bwa_aln_pe.sh' : 'bwa_aln_se.sh')
			break;
		case "bowtie-chip":
			template (r_tag == 'pe' ? 'bowtie2_pe.sh' : 'bowtie2_se.sh')
			break;
		default: 
			error "Aligning with ${params.aligner} is not implemented. You can add it to 'align_reads' process"
			break;
	}
}


wasp_path = '/opt/WASP' // Path in the container

def get_container(file_name) {
  parent = file(file_name).parent
  container = "--bind ${parent}"
  if (file(file_name).exists()) {
	old_parent = file(file_name).toRealPath().parent
	if (old_parent != parent) {
		container += ",${old_parent}"
	}
  } 
  return container
}

process filter_variants {
	tag "${indiv_id}"
	container "${params.container}"
	containerOptions "${get_container(params.genotype_file)} ${get_container(params.nuclear_chroms)}" 
	publishDir "${params.outdir}/target_variants"

	input:
		val indiv_id

	output:
		tuple val(indiv_id), path(outname), path("${outname}.tbi")

	script:
	outname = "${indiv_id}.bed.gz"
	"""
    echo -e "#chr\tstart\tend\tID\tref\talt\tAAF\tRAF\tGT\tGQ\tDP\tref_count\talt_count\tindiv_id" > file.tmp.bed
	bcftools query \
		-s ${indiv_id} \
		-i'GT="alt"' \
		-f"%CHROM\t%POS0\t%POS\t%ID\t%REF\t%ALT\t%INFO/AAF\t%INFO/RAF\t[%GT\t%GQ\t%DP\t%AD{0}\t%AD{1}]\t%SAMPLE\n" \
		${params.genotype_file} \
        | awk -v OFS="\\t" \
            -v min_GQ=${params.min_GQ} \
            -v min_AD=${params.min_AD} \
            -v min_DP=${params.min_DP} \
            '\$10<min_GQ { next; } \
            \$11<min_DP { next; } \
            (\$9=="0/1" || \$9=="1/0" || \$9=="0|1" || \$9=="1|0") && (\$12<min_AD || \$13<min_AD) { next; } \
            { print; }' \
        | awk 'NR==FNR { ok[\$1]=1; next } (\$1 in ok)' \
            ${params.nuclear_chroms} \
            - \
        | sort-bed - >> file.tmp.bed
        
        bgzip -c file.tmp.bed > ${outname}
        tabix -f -p bed "${outname}"
	"""
}

process collect_snv_files {
	scratch true
	container "${params.container}"
	publishDir "${params.outdir}"

	input:
		path snv_files

	output:
		tuple path(name), path("${name}.tbi")

	script:
	name = "all_variants_stats.bed.gz"
	"""
    printf '%s\n' ${snv_files} > filelist.txt
    first=\$(head -n 1 filelist.txt)
    zcat "\$first" | grep -m1 '^#' > tmp.bed

    while read file; do
        zcat "\$file" \
            | grep -v '^#' \
    done < filelist.txt \
        | sort-bed - >> "tmp.bed"

	bgzip -c tmp.bed > ${name}
	tabix ${name}
	"""
}

process generate_h5_tables {
	scratch true
	container "${params.container}"
	containerOptions "${get_container(params.genotype_file)} ${get_container(params.chrom_sizes)}"

	output:
		path '*.h5'

	script:
	"""
	for chrom in `tabix -l ${params.genotype_file}`; do
		bcftools view -r \${chrom} -Oz ${params.genotype_file} > \${chrom}.vcf.gz
		bcftools index \${chrom}.vcf.gz
	done

	gzip -c ${params.chrom_sizes} > chrom_sizes.txt.gz

	${wasp_path}/snp2h5/snp2h5 \
		--chrom chrom_sizes.txt.gz \
		--format vcf \
		--haplotype haplotypes.h5 \
		--snp_index snp_index.h5 \
		--snp_tab snp_tab.h5 \
		chr*.vcf.gz
	"""
}

process split_reads {
	tag "${sample_id}:${r_tag}"
	container "${params.container}"
	containerOptions "${get_container(params.genome_fasta_file)}"

	input:
		tuple val(sample_id), val(indiv_id), path(bam_file), path(bam_index_file), val(r_tag)

	output:
		tuple val(sample_id), val(indiv_id), val(r_tag), path(name), path("${name}.bai"), env(n_counts)

	script:
	name = "${sample_id}.${r_tag}.bam"
	pars = r_tag == 'pe' ? "-f 1" : "-F 1"
	"""
	samtools view -O bam -h ${pars} \
		--reference ${params.genome_fasta_file} ${bam_file} > ${name}
	samtools index ${name}
	n_counts=\$(samtools view -c ${name})
	"""
}

process extract_to_remap_reads {
	tag "${sample_id}:${r_tag}"
	container "${params.container}"
	cpus 2
	scratch true

	input:
		tuple val(sample_id), val(indiv_id), val(r_tag), path(bam_file), path(bam_file_index), env(n_counts)
		path h5_tables

	output:
		tuple val(sample_id), val(r_tag), path(out_bam_file), path("${out_bam_file}.bai"), emit: bamfile
		tuple val(sample_id), val(r_tag), path(fasta1), path(fasta2), emit: fastq
	
	script:
	name = "${sample_id}.${r_tag}.rmdup"
	out_bam_file = "${name}.bam"
	fasta1 = "${name}.remap.fq1.gz"
	fasta2 = "${name}.remap.fq2.gz"
	if (r_tag == 'pe') {
		"""
		python3 ${wasp_path}/mapping/rmdup_pe.py \
			${bam_file} pe.reads.rmdup.bam

		samtools sort \
			-@${task.cpus} \
			-o ${out_bam_file} \
			-O bam \
			pe.reads.rmdup.bam

		samtools index ${out_bam_file}

		python3 ${wasp_path}/mapping/find_intersecting_snps.py \
			--is_paired_end \
			--is_sorted \
			--output_dir ./ \
			--snp_tab snp_tab.h5 \
			--snp_index snp_index.h5  \
			--haplotype haplotypes.h5 \
			--samples ${indiv_id} \
			${out_bam_file}
		"""
	} else {
		"""
		# an ugly hack to deal with repeated read names on legacy SOLEXA GA1 data
		python3 $moduleDir/bin/hash_se_reads.py ${bam_file} se.hashed.bam

		python3 ${wasp_path}/mapping/rmdup.py \
			se.hashed.bam  se.reads.rmdup.bam
		
		samtools sort \
			-@${task.cpus} \
			-o ${out_bam_file} \
			-O bam \
			se.reads.rmdup.bam
		
		samtools index ${out_bam_file}

		### Creates 3 following files:
		### se.reads.rmdup.sorted.to.remap.bam (reads to remap)
		### se.reads.rmdup.sorted.keep.bam (reads to keep)
		### se.reads.rmdup.sorted.remap.fq.gz (fastq file containing the reads with flipped alleles to remap)
		python3 ${wasp_path}/mapping/find_intersecting_snps.py \
			--is_sorted \
			--output_dir ./ \
			--snp_tab snp_tab.h5 \
			--snp_index snp_index.h5  \
			--haplotype haplotypes.h5 \
			--samples ${indiv_id} \
			${out_bam_file}
		
		mv ${name}.remap.fq.gz ${fasta1}
		ln -s ${fasta1} ${fasta2}
		"""
	}
}


process wasp_filter_reads {
	container "${params.container}"
	scratch true
	tag "${sample_id}:${r_tag}"

	input:
		tuple val(sample_id), val(r_tag), path(bam_file), path(initial_bam_file), path(initial_bam_file_index)
	
	output:
		tuple val(sample_id), path(name)

	script:
	name = "${sample_id}.${r_tag}.result.bam"
	"""
	python3 ${wasp_path}/mapping/filter_remapped_reads.py \
		${initial_bam_file} \
		${bam_file} \
		${name}
	"""
}

process merge_bam_files {
	container "${params.container}"
	scratch true
	tag "${prefix}:${sample_id}"
	cpus 2
	publishDir "${params.outdir}/filtered_bam", pattern: "${name}*"

	input:
		tuple val(prefix), val(sample_id), path(bam_files)

	output:
		tuple val(prefix), val(sample_id), path(name), path("${name}.bai")

	script:
	name = "${sample_id}.${prefix}.merged.bam"
	// There is probably a cleaner way, but this will do for now
	non_empty_bam_files = bam_files.stream().filter(
			s -> s.name != 'empty.bam'
		).toArray()
	if (non_empty_bam_files.size() >= 2)
		"""
		samtools merge -f reads.rmdup.original.bam \
			${non_empty_bam_files.join(' ')}

		samtools sort \
			-@${task.cpus} \
			-o ${name} \
			reads.rmdup.original.bam
		samtools index ${name}
		"""
	else
		"""
		ln -s ${non_empty_bam_files.join(' ')} reads.passing.bam
		samtools sort \
			-@${task.cpus} \
			-o ${name}  \
			reads.passing.bam
		samtools index ${name}
		"""
}

process count_initial_reads {
	container "${params.container}"
	tag "${sample_id}"

	input:
		tuple val(sample_id), path(bam_file), path(bam_file_index), path(filtered_sites_file), path(filtered_sites_file_index)

	output:
		tuple val(sample_id), path(name)

	script:
	name = "${sample_id}.coverage.tsv.gz"
	"""
	python3 $moduleDir/bin/count_tags_pileup.py \
		${filtered_sites_file} \
		${bam_file} \
        ${name} \
		--only_coverage
	"""
}

process count_remapped_reads {
	tag "${sample_id}"
	container "${params.container}"
	publishDir "${params.outdir}/remapped_files"

	input:
		tuple val(sample_id), path(rmdup_counts), path(bam_passing_file), path(bam_passing_file_index), path(filtered_sites_file), path(filtered_sites_file_index)

	output:
		tuple val(sample_id), path(name), path("${name}.tbi")

	script:
	name = "${sample_id}.bed.gz"
	"""
	python3 $moduleDir/bin/count_tags_pileup.py \
		${filtered_sites_file} \
        ${bam_passing_file} \
        ${name} \
		--original_dedup_cover ${rmdup_counts}
	"""
}

process reformat2babachi {
	publishDir "${params.outdir}/babachi_files"
	tag "${sample_id}"
	container "${params.container}"

	input:
		tuple val(sample_id), path(bed_file), path(bed_file_index)

	output:
		tuple val(sample_id), path(name)

	script:
	name = "${sample_id}.snps.bed"
	"""
	python3 $moduleDir/bin/tags_to_babachi_format.py ${bed_file} ${name}
	"""
}

process add_snp_files_to_meta {
	publishDir "${params.outdir}"
	container "${params.container}"

    input:
        path samples_file

	output:
		path name

	script:
	name = "meta+sample_ids.tsv"
	"""
	python3 $moduleDir/bin/add_meta.py \
		${samples_file} \
		${name} \
		${params.outdir}/babachi_files
	"""
}

workflow waspRealigning {
	take:
		samples_aggregations // sample_id, indiv_id, bam, bam_index
	main:
		r_tags = Channel.of('pe', 'se') // flag that indicates paired-end/single-end reads

		h5_tables = generate_h5_tables().collect(sort: true) // h5 files
		snps_sites = samples_aggregations
			| map(it -> it[1]) // indiv_id
			| unique { it } 
			| filter_variants // indiv_id, variants, variants_index

		snps_sites
			| map(it -> it[1]) 
			| collect(sort: true) // varaints files
			| collect_snv_files
		
		snp_sites_by_sample_id = samples_aggregations
			| map(it -> tuple(it[1], it[0])) // indiv_id, sample_id
			| combine(snps_sites, by: 0) // indiv_id, sample_id, variants, variants_index
			| map(it -> tuple(*it[1..(it.size()-1)])) // sample_id, variants, variants_index
		
		split_rs = samples_aggregations
			| combine(r_tags) // sample_id, indiv_id, bam, bam_index, r_tag
			| split_reads // sample_id, indiv_id, r_tag, r_tag_bam, bam_index, read_count
			| branch {
				files: it[5].toInteger() > 0
        		nodata: true
			}
		
		to_remap_reads_and_initial_bam = extract_to_remap_reads(
            split_rs.files,
            h5_tables
        ) 
				// bamfile: sample_id, r_tag, r_tag_bam_dedup, bam_inex
				// fastq: sample_id, r_tag, fastq1, fastq2

		dedup_bam = to_remap_reads_and_initial_bam.bamfile // sample_id, r_tag, r_tag_bam_dedup, bam_inex

		filtered_bam = to_remap_reads_and_initial_bam.fastq // sample_id, r_tag, fastq1, fastq2
			| align_reads // sample_id, r_tag, realigned_bam
			| join(dedup_bam, by: [0, 1]) // sample_id, r_tag, realigned_bam, r_tag_bam_dedup, bam_index
			| wasp_filter_reads // sample_id, filtered_bam
			| map(it -> tuple('remapped', *it)) // type, sample_id, filtered_bam

		nodata = split_rs.nodata
			| map(it -> tuple(it[0], file('empty.bam'))) // sample_id, 'empty.bam'
			| multiMap { it ->
				remapped: tuple('remapped', *it)
				initial: tuple('initial', *it)
			} // type, sample_id, filtered_bam

		merged_out_bam = dedup_bam
			| map(it -> tuple('initial', it[0], it[2]))
			| mix(filtered_bam, nodata.initial, nodata.remapped) // type, sample_id, bam
			| groupTuple(size: 2, by: [0, 1])
			| merge_bam_files // type, sample_id, bam, bam_index
			| branch {
				initial: it[0] == 'initial'
					return tuple(*it[1..(it.size()-1)])
				remapped: it[0] == 'remapped'
					return tuple(*it[1..(it.size()-1)])
			} // sample_id, bam, bam_index

		out = merged_out_bam.initial // sample_id, initial_bam, bam_index
			| join(snp_sites_by_sample_id) // sample_id, initial_bam, bam_index, variants, variants_index
			| count_initial_reads // sample_id, initial_counts, counts_index
			| join(merged_out_bam.remapped) // sample_id, initial_counts, counts_index, remapped_bam, bam_index, variants, variants_index
			| join(snp_sites_by_sample_id) // sample_id, initial_counts, counts_index, remapped_bam, bam_index, variants, variants_index
			| count_remapped_reads // sample_id, summary_file, summary_file_index
			| reformat2babachi // sample_id, babachi_formatted_summary_file

	emit:
		out
}


workflow {
	samples_aggregations = Channel.fromPath(params.samples_file)
		| splitCsv(header:true, sep:'\t')
		| map(row -> tuple(row.sample_id, row.indiv_id,
                file(row.cram_file),
                file(row?.cram_index ?: "${row.cram_file}.crai")
            )
        )
		| filter { !it[1].isEmpty() }


	samples_aggregations
		| map(it -> it[1])
		| unique()
		| count()
		| view {
			it -> """There are ${it} unique INDIV_IDs in the ${params.samples_file}. Please, check that they correspond to IDs in ${params.genotype_file}"""
		}

	add_snp_files_to_meta(params.samples_file) 
	samples_aggregations | waspRealigning
}


workflow reformat {
	add_snp_files_to_meta(params.samples_file)
	params.remapped_files_dir = "${params.outdir}/remapped_files/"
	Channel.fromPath("${params.remapped_files_dir}/*.bed.gz")
		| map(
			it -> tuple(
				it.name.replaceAll('.bed.gz', ''), 
				it, 
				file("${it.parent}/${it.name}.tbi")
			))
		| reformat2babachi
}


process merge_by_indiv_and_cell_type {
	tag "${indiv_id}:${cell_type}"
	publishDir "${params.outdir}/merged"

	container "${params.container}"
	containerOptions "${get_container(params.genome_fasta_file)}"
	cpus 2

	input:
		tuple val(indiv_id), val(cell_type), val(bam_files)

	output:
		tuple val(sample_id), val(indiv_id), path(name), path("${name}.crai")

	script:
	sample_id = "${indiv_id}_${cell_type}"
	name = "${sample_id}.cram"
	"""
	samtools merge -f -@${task.cpus} \
		-O CRAM --write-index \
		--reference ${params.genome_fasta_file} \
		${name} ${bam_files}
	"""


}

workflow mergeByIndivAndCellType{
	bams_grouped_by_indiv_and_cell_type = Channel.fromPath(params.samples_file)
		| splitCsv(header:true, sep:'\t')
		| map(row -> tuple(row.indiv_id, row.cell_type, file(row.bam_file)))
		| groupTuple(by: [0,1])
		| map(it -> tuple(it[0], it[1], it[2].flatten().join(" ")))
		| merge_by_indiv_and_cell_type
		| waspRealigning

	add_snp_files_to_meta() 
}