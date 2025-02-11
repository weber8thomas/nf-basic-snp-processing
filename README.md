# nf-basic-snp-processing

## Command

```
nextflow run main.nf -profile slurm --samplesheet samplesheet.tsv --outdir /g/korbel2/weber/MosaiCatcher_files/DEMULTIPLEXING_POOLS/DEMULTIPLEXING_POOLS/1000G_SNV_with_GT/NON_REFERENCES_SAMPLES_FINAL_VCF -process.executor slurm -process.clusterOptions '-A datasci' -executor.queueSize 1000 -resume
```