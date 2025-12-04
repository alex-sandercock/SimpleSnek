import os

# Define the samples and their respective input files

#FASTQ directory
fastq_dir = "/workdir/ams866/Alfalfa/validation_plates/FASTQ"

#samples = ["B"]  # Add your sample names here
#Use this instead if you want to use all of the FASTQ in a directory as the sample names
# Get the list of sample names from the filenames in the directory
samples = [os.path.splitext(file)[0] for file in os.listdir(fastq_dir) if file.endswith(".FASTQ.gz")]
samples = [sample.replace(".FASTQ", "") for sample in samples]

# Define a default rule to execute the entire pipeline
rule all:
    input:
        "final_combined_filtered_snps.vcf.gz",
        "final_combined_filtered_indels.vcf.gz"

# Additional target rules for partial pipeline execution
rule map_all:
    input:
        expand("../../mapped/{sample}.bam", sample=samples)

rule sort_all:
    input:
        expand("../../sorted/{sample}.sorted.bam", sample=samples)

rule mark_all:
    input:
        expand("../../marked/{sample}.marked.bam", sample=samples)

rule call_all:
    input:
        expand("../../vcf/{sample}.g.vcf", sample=samples)

# Rule to index the genome
#rule index_genome:
#    input:
#        "../../reference_genome/final_genome.fa"
#    output:
#        "final_genome.fa.bwt",
#        "final_genome.fa.pac",
#        "final_genome.fa.ann",
#        "final_genome.fa.sa"
#    shell:
#        """
#        bwa index {input}
#        java -jar ../../bin/picard.jar CreateSequenceDictionary R={input} O=final_genome.dict
#        samtools faidx {input}
#        """

#Define a rule to trim raw fastq files
rule trim_fastq:
    input:
        R1 = fastq_dir + "/{sample}.FASTQ.gz",
    output:
        "../../trim/{sample}.trim.fastq.gz"
    shell:
        "java -jar ../../bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 {input.R1} {output} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36"

# Define a rule for mapping with BWA
rule bwa_mapping:
    input:
        R = "../../trim/{sample}.trim.fastq.gz"
    params:
        reference = "../../reference_genome/final_genome.fa",
        cores = 10,
        rg = "'@RG\\tID:Sample_{sample}\\tSM:{sample}\\tLB:{sample}\\tPL:ILLUMINA'"
    output:
        "../../mapped/{sample}.bam"
    shell:
        "bwa mem -R {params.rg} -t {params.cores} {params.reference} {input.R} | samtools view -@ 4 -Sb -F 4 - > {output}"

# Define a rule for sorting and indexing with Samtools
rule samtools_sort_index:
    input:
        "../../mapped/{sample}.bam"
    output:
        "../../sorted/{sample}.sorted.bam"
    shell:
        "samtools sort {input} -o {output} && samtools index {output}"

# Define a rule for removing PCR duplicates
rule mark_duplicates:
    input:
        "../../sorted/{sample}.sorted.bam"
    output:
        "../../marked/{sample}.marked.bam"
    shell:
        """
        java -Xmx20g -jar ../../bin/picard.jar MarkDuplicates \
        I={input} \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true \
        O={output} \
        M={wildcards.sample}_marked_dup_metrics.txt
        """

# Define a rule for SNP calling with GATK
rule gatk_snp_calling:
    input:
        bam = "../../marked/{sample}.marked.bam"
    params:
        java_option = "-Xmx20g",
        reference = "../../reference_genome/final_genome.fa",
        ploidy = 4
    output:
        "../../vcf/{sample}.g.vcf"
    shell:
        "../../bin/gatk-4.3.0.0/gatk --java-options {params.java_option} HaplotypeCaller -ploidy {params.ploidy} --native-pair-hmm-threads 2 -R {params.reference} -I {input.bam} -ERC GVCF -O {output}"

#Define a rule for making a list of GVCF file names
rule gvcf_list:
    input:
        gvcfs = expand("../../vcf/{sample}.g.vcf", sample=samples)
    output:
        "../../vcf/sample_files.list"
    shell:
        """
        echo {input.gvcfs} | tr ' ' '\\n' > {output}
        """

#Define a rule for combining GVCF files for all samples
rule gatk_merge:
    input:
        vcf_files = expand("../../vcf/{sample}.g.vcf", sample=samples),
        reference = "../../reference_genome/final_genome.fa",
        list = "../../vcf/sample_files.list"
    output:
        "combined_gvcf.vcf"
    shell:
        "../../bin/gatk-4.3.0.0/gatk CombineGVCFs -R {input.reference} -V {input.list} -O {output}"

# Define a rule for joint genotyping with GATK
rule gatk_joint_genotyping:
    input:
        combined_vcf = "combined_gvcf.vcf",
        reference = "../../reference_genome/final_genome.fa"
    output:
        "joint_genotyped.vcf"
    shell:
        "../../bin/gatk-4.3.0.0/gatk GenotypeGVCFs -R {input.reference} -V {input.combined_vcf} -O {output}"

#Define a rule to filter SNPs
rule gatk_select_snps:
    input:
        joint_vcf = "joint_genotyped.vcf",
        reference = "../../reference_genome/final_genome.fa"
    output:
        "combined_snp_raw.vcf.gz"
    shell:
        "../../bin/gatk-4.3.0.0/gatk SelectVariants -R {input.reference} -V {input.joint_vcf} --select-type-to-include SNP -O {output}"

#Define a rule to filter indels
rule gatk_select_indels:
    input:
        joint_vcf = "joint_genotyped.vcf",
        reference = "../../reference_genome/final_genome.fa"
    output:
        "combined_indel_raw.vcf.gz"
    shell:
        "../../bin/gatk-4.3.0.0/gatk SelectVariants -R {input.reference} -V {input.joint_vcf} --select-type-to-include INDEL -O {output}"


# Define a rule for filtering with BCFtools
rule bcftools_filter:
    input:
        raw_snps = "combined_snp_raw.vcf.gz",
        raw_indels = "combined_indel_raw.vcf.gz"
    params:
        filter= "'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'"
    output:
        snps = "final_combined_filtered_snps.vcf.gz",
        indels = "final_combined_filtered_indels.vcf.gz"
    shell:
        """
        bcftools filter -e {params.filter} {input.raw_snps} -Oz -o {output.snps}
        bcftools filter -e {params.filter} {input.raw_indels} -Oz -o {output.indels}
        """
#Defile a rule for filtering missing data with VCFtools?


# Define a rule to clean up intermediate files (optional)
#rule cleanup:
#    input:
#        expand("sorted/{sample}.sorted.bam", sample=samples),
#        "joint_genotyped.vcf"
#    shell:
#        "rm -f {input}"
