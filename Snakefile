# count number of reads
rule countreads:
  output: "{indir}.{file}.fq.count"
  input: "{indir}/{file}.fq"
  shell:
    "echo $(( $(wc -l <{input}) / 4 )) > {output}"

# filter low quality reads less than 20 and shorter than 100 > edited to quality to 22
rule trimreads:
  output:"trimmed/{file}.fq"
  input:"reads/{file}.fq"
  shell:
    "fastq_quality_trimmer -t 22 -l 100 -o {output} < {input}"
  
# create genome index for kallisto
rule kallisto_index:
  output:
    idx = "transcriptome/{strain}.kallisto_index",
    log   = "transcriptome/{strain}.kallisto_log"
  input: "transcriptome/{strain}.cdna.all.fa.gz"
  shell:
    "kallisto index -i {output.idx} {input} >& {output.log}"  
  
# align reads
rule kallisto_quant:
  output:
    h5   = "kallisto.{sample}/abundance.h5",
    tsv  = "kallisto.{sample}/abundance.tsv",
    json = "kallisto.{sample}/run_info.json"
  input:
    idx  = "transcriptome/Saccharomyces_cerevisiae.R64-1-1.kallisto_index",
    fq1  = "trimmed/{sample}_1.fq",
    fq2  = "trimmed/{sample}_2.fq"
  shell:
    "kallisto quant -i {input.idx} -o kallisto.{wildcards.sample} {input.fq1} {input.fq2}"
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    