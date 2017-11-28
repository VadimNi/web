IDS, = glob_wildcards("./input/{id}.bam")

rule all:
    input:
        expand("./output/reports/{id}", id=IDS)

rule reports:
    input:
        "./output/{id}_for_0_1.v",
        "./output/{id}_for_0_2.v",
        "./output/{id}_for_1_2.v",
        "./output/haplotypecaller/{id}.v",
        "./output/freebayes/{id}.v",
        "./output/samtools/{id}.v"
    output:
        "./output/reports/{id}"
    run:
        for j in range(len(output)):
            a = []
            print(input)
            for i in range(6):
                with open(input[i], 'r') as inp:
                    for line in inp:
                        line = line.rstrip()
                        a.append(line)
            with open(output[j], 'w') as out:
                out.write("{0}\t{1}\t{2}\n".format(a[4],a[0],a[1]))
                out.write("{0}\t{1}\t{2}\n".format(a[3],None,a[2]))
                out.write("{0}\t{1}\t{2}".format(a[5],None,None))

rule count_intersect_variants:
    input:
        "./output/{id}_for_0_1.vcf.gz",
        "./output/{id}_for_0_2.vcf.gz",
        "./output/{id}_for_1_2.vcf.gz"
    output:
        "./output/{id}_for_0_1.v",
        "./output/{id}_for_0_2.v",
        "./output/{id}_for_1_2.v"
    shell:
        "vcftools --gzvcf {input[0]} 2>&1 | cut -d' ' -f9 | xargs > {output[0]} ; "
        "vcftools --gzvcf {input[1]} 2>&1 | cut -d' ' -f9 | xargs > {output[1]} ; "
        "vcftools --gzvcf {input[2]} 2>&1 | cut -d' ' -f9 | xargs > {output[2]} "

rule intersect_vcf:
    input:
        "./output/freebayes/{id}.vcf.gz",
        "./output/haplotypecaller/{id}.vcf.gz",
        "./output/samtools/{id}.vcf.gz"
    output:
        "./output/{id}_for_0_1.vcf.gz",
        "./output/{id}_for_0_2.vcf.gz",
        "./output/{id}_for_1_2.vcf.gz"
    shell:
        "vcf-isec -f -n +2 {input[0]} {input[1]} | bgzip -c > {output[0]} ; "
        "vcf-isec -f -n +2 {input[0]} {input[2]} | bgzip -c > {output[1]} ; "
        "vcf-isec -f -n +2 {input[1]} {input[2]} | bgzip -c > {output[2]} "

rule count_individual_variants:
    input:
        "./output/haplotypecaller/{id}.vcf.gz",
        "./output/freebayes/{id}.vcf.gz",
        "./output/samtools/{id}.vcf.gz"
    output:
        "./output/haplotypecaller/{id}.v",
        "./output/freebayes/{id}.v",
        "./output/samtools/{id}.v"
    shell:
        "vcftools --gzvcf {input[0]} 2>&1 | cut -d' ' -f9 | xargs > {output[0]} ; "
        "vcftools --gzvcf {input[1]} 2>&1 | cut -d' ' -f9 | xargs > {output[1]} ; "
        "vcftools --gzvcf {input[2]} 2>&1 | cut -d' ' -f9 | xargs > {output[2]} "

rule index_vcf:
    input:
        "./output/haplotypecaller/{id}.vcf",
        "./output/freebayes/{id}.vcf",
        "./output/samtools/{id}.vcf"
    output:
        "./output/haplotypecaller/{id}.vcf.gz",
        "./output/freebayes/{id}.vcf.gz",
        "./output/samtools/{id}.vcf.gz"
    shell:
        "bgzip {input[0]} ; "
        "tabix -p vcf {output[0]} ; "
        "bgzip {input[1]} ; "
        "tabix -p vcf {output[1]} ; "
        "bgzip {input[2]} ; "
        "tabix -p vcf {output[2]}"

rule vc_samtools:
    input:
        "./data/22.fa",
        "./input/{id}.bam",
        "./data/22.fa.fai"
    output:
        "./output/samtools/{id}.vcf"
    shell:
        "samtools mpileup -uf {input[0]} {input[1]} | bcftools view -vcg - > {output}"

rule vc_freebayes:
    input:
        "./data/22.fa",
        "./input/{id}.bam",
        "./data/22.fa.fai"
    output:
        "./output/freebayes/{id}.vcf"
    shell:
        "freebayes -f {input[0]} {input[1]} > {output}"

rule vc_haplotypecaller:
    input:
        "./data/22.fa",
        "./input/{id}.bam",
        "./data/22.fa.fai",
        "./data/22.dict",
        "./input/{id}.bam.bai"
    output:
        "./output/haplotypecaller/{id}.vcf"
    shell:
        "java -jar $GATK -R {input[0]} -T HaplotypeCaller -I {input[1]} -o {output}"


rule index_bam:
    input:
        "./input/{id}.bam"
    output:
        "./input/{id}.bam.bai"
    shell:
        "samtools index {input}"

rule reference_dict:
    input:
        "./data/22.fa"
    output:
        "./data/22.dict"
    shell:
        "java -jar $PICARD CreateSequenceDictionary R={input} O={output}"

rule reference_index:
    input: 
        "./data/22.fa"
    output:
        "./data/22.fa.fai"
    shell:
        "samtools faidx {input}"

rule folders:
    output:
        "./output/haplotypecaller",
        "./output/freebayes",
        "./output/samtools",
        "./output/reports"
    shell: "mkdir {output[0]} {output[1]} {output[2]} {output[2]}"
