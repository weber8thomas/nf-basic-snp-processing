nextflow.enable.dsl = 2

// Validate required parameters
if (!params.samplesheet) {
    error "Please provide a samplesheet with --samplesheet"
}
if (!params.ref_fasta) {
    error "Please provide a reference fasta with --ref_fasta"
}
if (!params.ref_fai) {
    error "Please provide a reference index with --ref_fai"
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

process BCF_FILTER_AND_NORM {
    tag { sample }
    conda "bioconda::bcftools bioconda::tabix"
    publishDir "${params.outdir}/normalized", mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(tbi)
    path(ref_fasta)
    path(ref_fai)

    output:
    tuple val(sample), path("${sample}.normalized.vcf.gz"), path("${sample}.normalized.vcf.gz.tbi"), emit: filtered_normalized
    
    script:
    """
    # First filter SNPs
    bcftools view -v snps \
        -i 'QUAL >= 20 && TYPE="snp" && GT!="0/0" && GT!="./." && INFO/AC>0' \
        $vcf > filtered.vcf
        
    # Split multi-allelic sites
    bcftools norm -m-any filtered.vcf > split.vcf
    
    # Ensure single-base REF/ALT
    awk '(/^#/ || (length(\$4)==1 && length(\$5)==1))' split.vcf > clean.vcf
    
    # Compress and index
    bgzip -c clean.vcf > ${sample}.normalized.vcf.gz
    tabix -p vcf ${sample}.normalized.vcf.gz
    """
}

workflow {
    // Get reference channels
    ref_fasta = Channel.fromPath(params.ref_fasta)
    ref_fai = Channel.fromPath(params.ref_fai)
    
    // Execute pipeline
    filtered_normalized = BCF_FILTER_AND_NORM(ch_samples, ref_fasta, ref_fai)
    // normalized = NORMALIZE(filtered.filtered, ref_fasta, ref_fai)
    
    // Show results
    filtered_normalized.filtered_normalized.view()
}