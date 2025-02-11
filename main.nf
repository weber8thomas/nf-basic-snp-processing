nextflow.enable.dsl = 2

def valid_chroms = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                    'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                    'chr20', 'chr21', 'chr22', 
                    'chrX'
                    ]

// Validate required parameters
if (!params.samplesheet) {
    error "Please provide a samplesheet with --samplesheet"
}

// Create input channel
ch_samples = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true, sep:'\t')
    .map { row -> [ 
        row.sample, 
        file(row.vcf),
        file("${row.vcf}.tbi")
    ]}

// Create 1000G VCFs channel from the map
ch_onekg_vcfs = Channel
    .from(params.onekg_vcfs.collect { chrom, vcf_path -> 
        [
            chrom,                     // Chromosome (e.g., "chr1")
            file(vcf_path),           // VCF file
            file("${vcf_path}.tbi")   // Index file
        ]
    })

process BCF_FILTER_AND_NORM {
    tag { sample }
    conda "bioconda::bcftools bioconda::tabix"
    publishDir "${params.outdir}/normalized", mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(tbi)

    output:
    tuple val(sample), path("${sample}.normalized.vcf.gz"), path("${sample}.normalized.vcf.gz.tbi"), emit: filtered_normalized
    
    script:
    // def is_ont = vcf.toString().contains('ONT-clair3')
    // def filter_cmd = is_ont ? 
    //     "bcftools view -v snps -i 'QUAL >= 20 && GT==\"0/1\" && DP>=5'" :
    //     ""
    
    """
    # Filter SNPs with appropriate criteria
    # bcftools view -v snps -i 'QUAL >= 20 && GT==\"0/1\" && FORMAT/DP>=5' $vcf > filtered.vcf
    # bcftools view -v snps -i 'QUAL >= 20 && FORMAT/DP>=5' $vcf > filtered.vcf
    bcftools view -v snps -i 'QUAL >= 20 && FORMAT/DP>=5' ${vcf} > filtered.vcf  # Changed from vcf.tbi to vcf

    # Split multi-allelic sites
    bcftools norm -m-any filtered.vcf > split.vcf
    
    # Ensure single-base REF/ALT
    awk '(/^#/ || (length(\$4)==1 && length(\$5)==1))' split.vcf > clean.vcf
    
    # Compress and index
    bgzip -c clean.vcf > ${sample}.normalized.vcf.gz
    tabix -p vcf ${sample}.normalized.vcf.gz
    """
}

// First, create a helper process to get chromosomes
process GET_CHROMS {
    tag { sample }
    conda "bioconda::bcftools"
    
    input:
    tuple val(sample), path(vcf), path(tbi)
    
    output:
    tuple val(sample), path(vcf), path(tbi), stdout, emit: chroms
    
    script:
    """
    bcftools query -f '%CHROM\n' $vcf | sort | uniq
    """
}

// Then split process handles one chromosome at a time
process SPLIT_VCF_BY_CHROMOSOME {
    tag { "${sample}_${chrom}" }
    conda "bioconda::bcftools bioconda::tabix"
    publishDir "${params.outdir}/split_by_chrom", mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(tbi), val(chrom)

    output:
    tuple val(sample), val(chrom), path("${sample}.${chrom}.vcf.gz"), path("${sample}.${chrom}.vcf.gz.tbi"), emit: split_vcfs

    script:
    """
    bcftools view -I $vcf \
        --regions ${chrom} \
        | bcftools annotate -x INFO \
        | awk '{if(\$8==".") \$8="AC=1000;AF=0"; print}' OFS='\t' \
        | bcftools reheader -h <(bcftools view -h $vcf) \
        | bcftools view -O z \
        > ${sample}.${chrom}.vcf.gz
    tabix -p vcf ${sample}.${chrom}.vcf.gz
    """
}


process ANNOTATE_WITH_VCFANNO {
    tag { "${sample}_${chrom}" }
    conda "bioconda::vcfanno bioconda::tabix"
    publishDir "${params.outdir}/annotated_by_chrom", mode: 'copy'

    input:
    tuple val(sample), val(chrom), path(split_vcf), path(split_tbi), path(onekg_vcf), path(onekg_tbi)


    output:
    tuple val(sample), val(chrom), path("${sample}.${chrom}.annotated.vcf.gz"), path("${sample}.${chrom}.annotated.vcf.gz.tbi"), emit: annotated_vcfs

    script:
    """
    # Create temporary vcfanno configuration
    cat <<EOF > vcfanno.conf
    [[annotation]]
    file = "${onekg_vcf}"
    fields = ["AF", "AC"]
    ops = ["self", "self"]
    EOF

    # Annotate using vcfanno
    vcfanno -p ${task.cpus} vcfanno.conf ${split_vcf} | bgzip > ${sample}.${chrom}.annotated.vcf.gz
    tabix -p vcf ${sample}.${chrom}.annotated.vcf.gz
    """
}

// task.cpus

process MERGE_ANNOTATED_VCFS {
    tag { sample }
    conda "bioconda::bcftools bioconda::tabix"
    publishDir "${params.outdir}/merged_annotated", mode: 'copy'

    input:
    tuple val(sample), path(vcfs)

    output:
    tuple val(sample), path("${sample}.merged.annotated.vcf.gz"), path("${sample}.merged.annotated.vcf.gz.tbi"), emit: merged_vcfs

    script:
    """
    # Merge multiple VCFs
    bcftools concat \
        --output ${sample}.merged.annotated.vcf.gz \
        --output-type z \
        ${vcfs.join(' ')}

    # Index the merged VCF
    tabix -p vcf ${sample}.merged.annotated.vcf.gz
    """
}

    // #bcftools view -i 'GT!="0/0"' ${sample}.merged.vcf.gz | \
    // #bcftools view -i 'INFO/AF < 0.05' | \
    // # bgzip > ${sample}.merged.annotated.vcf.gz
    // bcftools view -i ${sample}.merged.vcf.gz | bgzip > ${sample}.merged.annotated.vcf.gz


process EXPORT_TABLE {
    tag { sample }
    conda "bioconda::bcftools bioconda::tabix"
    publishDir "${params.outdir}/tables", mode: 'copy'
    publishDir "/g/korbel2/weber/MosaiCatcher_files/DEMULTIPLEXING_POOLS/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE", mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(tbi)

    output:
    path("${sample}.txt.gz"), emit: table

    script:
    """
    # Write header
    echo -e "ID\\tAC\\tAF\\tGT\\tSAMPLE" > "${sample}.txt"

    # Extract fields including genotype and append sample name.
    # The %CHROM:%POS:%REF:%ALT creates a unique ID,
    # %INFO/AC and %INFO/AF extract the AC and AF values,
    # and [%GT] extracts the genotype from the sample column.
    bcftools query -f '%CHROM:%POS:%REF:%ALT\\t%INFO/AC\\t%INFO/AF\\t[%GT]\\n' "${vcf}" | \
        awk -v sample="${sample}" 'BEGIN {OFS="\\t"} {sub(/^chr/, "", \$1); print \$0, sample}' >> "${sample}.txt"

    # Compress table
    bgzip "${sample}.txt"
    """
}

workflow {
    // Step 1: Filter and normalize VCF files
    BCF_FILTER_AND_NORM(ch_samples)
    
    // Step 2: Get and process chromosomes
    GET_CHROMS(BCF_FILTER_AND_NORM.out.filtered_normalized)
        .map { sample, vcf, tbi, chroms -> 
            def chromList = chroms.toString().trim().split('\n')
            chromList.findAll { it in valid_chroms }.collect { chrom ->
                [ sample, vcf, tbi, chrom ]
            }
        }
        .flatten()
        .buffer(size: 4)
        .set { ch_split_inputs }
        
    // ch_split_inputs.view { "Filtered input: $it" }

    // Step 3: Continue with split and annotation
    split_vcfs = SPLIT_VCF_BY_CHROMOSOME(ch_split_inputs)
    
    split_vcfs.split_vcfs
        .map { [it[1], it] }  // Key by chromosome
        // .view { "Before combine: $it" }
        .combine(ch_onekg_vcfs.map { [it[0], it] }, by: 0)  // Join on chromosome key
        .map { key, left, right -> 
            [left[0], left[1], left[2], left[3], right[1], right[2]]
        }
        // .view { "After combine: $it" }
        .set { ch_to_annotate }
        
    ANNOTATE_WITH_VCFANNO(ch_to_annotate)

    
    // ANNOTATE_WITH_VCFANNO.out.annotated_vcfs.view {
    //     "STEP 5 - Annotation Output:\nsample: ${it[0]}\nchrom: ${it[1]}\nannotated_vcf: ${it[2]}\nannotated_tbi: ${it[3]}\n"
    // }

    // Step 6: Group annotated VCFs by sample for merging
    ANNOTATE_WITH_VCFANNO.out.annotated_vcfs
        .groupTuple(by: 0)
        .map { sample, chroms, vcfs, tbis -> 
            tuple(sample, vcfs, tbis)
        }
        // .view {
        //     "STEP 6 - Grouped for Merging:\nsample: ${it[0]}\nvcfs: ${it[1]}\ntbis: ${it[2]}\n"
        // }
        .set { ch_grouped_vcfs }

    // Step 7: Merge annotated VCFs
    MERGE_ANNOTATED_VCFS(ch_grouped_vcfs)
    
    // MERGE_ANNOTATED_VCFS.out.merged_vcfs.view {
    //     "STEP 7 - Merged Output:\nsample: ${it[0]}\nmerged_vcf: ${it[1]}\nmerged_tbi: ${it[2]}\n"
    // }

    // Step 8: Export to table format
    EXPORT_TABLE(MERGE_ANNOTATED_VCFS.out.merged_vcfs)
    
    // EXPORT_TABLE.out.table.view {
    //     "STEP 8 - Table Output:\ntable: ${it}\n"
    // }
}