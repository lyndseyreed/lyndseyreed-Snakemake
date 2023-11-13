configfile: "config.yaml"
#config.yaml must be in pwd

# Inputs to process
CONDITIONS = config["conditions"]
REPLICATES = config["replicates"]
READ       = ["1", "2"]    

# multiqc to gather information in yeast directory that involves quant 
rule multiqc:
  output:
    directory("multiqc_out")
  input:
    fastqc   = expand("reads.{cond}_{rep}_{read}_fastqc.zip", cond=CONDITIONS, rep=REPLICATES, read=READ),
    kallisto = expand("kallisto.{cond}_{rep}",               cond=CONDITIONS, rep=REPLICATES),
    salmon   = expand("salmon.{cond}_{rep}",                cond=CONDITIONS, rep=REPLICATES)
  shell:
    "multiqc . -o {output}"

rule all_counts:
  input: 
    untrimmed = expand("reads.{cond}_{rep}_{read}.fq.count", cond=CONDITIONS, rep=REPLICATES, read=READ),
    trimmed   = expand("reads.{cond}_{rep}_{read}.fq.count", cond=CONDITIONS, rep=REPLICATES, read=READ)
  output:
    untrimmed = "all_untrimmed_counts.txt",
    trimmed   = "all_trimmed_counts.txt"
  shell:
    "cat {input.untrimmed} > {output.untrimmed} ; cat {input.trimmed} > {output.trimmed}"

# count number of reads
rule countreads:
  output: "{indir}.{file}.fq.count"
  input: "{indir}/{file}.fq"
  shell:
    "echo $(( $(wc -l <{input}) / 4 )) > {output}"

# filter low quality reads less than 20 and shorter than 100 > edited to quality 22
rule trimreads:
  output:"trimmed/{file}.fq"
  input:"reads/{file}.fq"
  params:
    qual_thresh = config["trimreads_qual_thresh"],
    min_len     = config.get("trimreads_min_len", "100")
  shell:
    "fastq_quality_trimmer -t {params.qual_thresh} -l {params.min_len} -o {output} < {input}"
  
# create genome index for kallisto
rule kallisto_index:
  output:
    idx = "transcriptome/{strain}.kallisto_index",
    log   = "transcriptome/{strain}.kallisto_log"
  input: "transcriptome/{strain}.cdna.all.fa.gz"
  shell:
    "kallisto index -i {output.idx} {input} >& {output.log}"  
  
# kallisto align reads
rule kallisto_quant:
  output:
    directory("kallisto.{sample}")
  input:
    idx  = "transcriptome/Saccharomyces_cerevisiae.R64-1-1.kallisto_index",
    fq1  = "trimmed/{sample}_1.fq",
    fq2  = "trimmed/{sample}_2.fq"
  shell:
    """
      mkdir {output}
      kallisto quant -i {input.idx} \
      -o {output} {input.fq1} {input.fq2} \
      >& {output}/kallisto_quant.log
    """
# challenge 1 fast qc
rule fastqc:
  output: 
    html = "{indir}.{file}_fastqc.html",
    zip =  "{indir}.{file}_fastqc.zip"
  input: "{indir}/{file}.fq"
  shell:
    """
      fastqc -o . {input} 
      mv {wildcards.file}_fastqc.html {output.html} 
      mv {wildcards.file}_fastqc.zip {output.zip}
    """
    
# for shell you can do semi colon to connect or three """ and proper indent for each line as 1 command
# for shell can clean up long commands using """ and then a \ at the end of the line before making the next one same indentation 
    
    
    
# create genome index for salmon
rule salmon_index:
  output: directory("transcriptome/{strain}.salmon_index")
  input: "transcriptome/{strain}.cdna.all.fa.gz"
  params:
    kmer = config.get("salmon_kmer", "29")
  shell:
    "salmon index -t {input} -i {output} -k {params.kmer}"
   
# salmon alignment
rule salmon_quant:
  output: directory("salmon.{sample}")
  input:
    idx = "transcriptome/Saccharomyces_cerevisiae.R64-1-1.salmon_index",
    fq1 = "trimmed/{sample}_1.fq",
    fq2 = "trimmed/{sample}_2.fq"
  shell:
    """
      salmon quant -i {input.idx} -l A \
      -1 {input.fq1} -2 {input.fq2} \
      --validateMappings -o {output}
    """
    

    
    
    
    
    