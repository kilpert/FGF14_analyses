if config["downsample_n"]:

    print("Downsample FASTQ:", config["downsample_n"])

    rule input_downsample:
        input:
            lambda wildcards: config["sample"][wildcards.sample]
        output:
            "{results}/{ref}/input/{sample}.fastq.gz"
        params:
            n=config["downsample_n"],
            seed=config["seed"]
        conda:
            "../envs/seqtk.yaml"
        threads:
            2
        shell:
            "seqtk sample "
            "-s {params.seed} "
            "{input} "
            "{params.n} "
            "| pigz -9 -p {threads} "
            ">{output} "

else:

    rule input:
        input:
            lambda wildcards: config["sample"][wildcards.sample]
        output:
            "{results}/{ref}/input/{sample}.fastq.gz"
        params:
            n=config["downsample_n"],
            seed=config["seed"]
        conda:
            "../envs/seqtk.yaml"
        shell:
            "ln -s "
            "$(realpath {input}) "
            "{output} "


rule input_fastqc:
    input:
        "{results}/{ref}/input/{sample}.fastq.gz"
    output:
        "{results}/{ref}/input/fastqc/{sample}_fastqc.html"
    params:
        outdir="{results}/{ref}/input/fastqc",
        adapters="-a "+config["fastqc"]["adapters"] if config["fastqc"]["adapters"] else ""
    log:
        "{results}/{ref}/log/input_fastqc.{sample}.log"
    conda:
        "../envs/fastqc.yaml"
    threads:
        4
    shell:
        "export _JAVA_OPTIONS='-Xmx8g'; "
        "fastqc "
        "--threads {threads} "
        "{params.adapters} "
        "-o {params.outdir} "
        "{input} "
        ">{log} 2>&1 "


rule input_fastq_read_phred:
    input:
        "{results}/{ref}/input/{sample}.fastq.gz"
    output:
        "{results}/{ref}/input/fastq_read_phred/{sample}.fastq_read_phred.html"
    params:
        n=10
    log:
        "{results}/{ref}/log/input_fastq_read_phred.{sample}.log"
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/fastq_read_phred.py "
        "{input} "
        "{params.n} "
        "{output} "
        ">{log} 2>&1 "


rule input_read_counts:
    input:
        expand("{{results}}/{{ref}}/input/{sample}.fastq.gz", 
            sample=samples,
        )
    output:
        temp("{results}/{ref}/read_counts/input.read_counts.tsv")
    params:
        data_colname="input"
    shell:
        "( "
        """echo -e "sample\tprocess\tn"; """
        "for f in $(printf '%s\n' {input} | sort); do "
        "bname=$(basename $f | cut -d'.' -f1); "
        "n=$(zcat $f | sed -n '1~4p' | wc -l); "
        """echo -e "$bname\t{params.data_colname}\t$n"; """
        "done; "
        ") "
        ">{output} "

