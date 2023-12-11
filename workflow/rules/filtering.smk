## bbduk filtering

rule bbduk_filtering:
    input:
        "{results}/{ref}/input/{sample}.fastq.gz"
    output:
        "{results}/{ref}/filtering/bbduk/{sample}.fastq.gz"
    params:
        extra=config["bbduk"] if config["bbduk"] else ""
    log:
        "{results}/{ref}/log/filtering_bbduk.{sample}.log"
    conda:
        "../envs/bbmap.yaml"
    shell:
        "bbduk.sh "
        "{params.extra} "
        "in={input} "
        "out={output} "
        ">{log} 2>&1 "


rule bbduk_fastqc:
    input:
        "{results}/{ref}/filtering/bbduk/{sample}.fastq.gz"
    output:
        "{results}/{ref}/filtering/bbduk/fastqc/{sample}_fastqc.html"
    params:
        outdir="{results}/{ref}/filtering/bbduk/fastqc",
        adapters="-a "+config["fastqc"]["adapters"] if config["fastqc"]["adapters"] else ""
    log:
        "{results}/{ref}/log/bbduk_fastqc.{sample}.log"
    conda:
        "../envs/fastqc.yaml"
    threads:
        4
    shell:
        "fastqc "
        "--threads {threads} "
        "{params.adapters} "
        "-o {params.outdir} "
        "{input} "
        ">{log} 2>&1 "


rule filtering_bbduk_read_counts:
    input:
        expand("{{results}}/{{ref}}/filtering/bbduk/{sample}.fastq.gz", 
            sample=samples,
        )
    output:
        temp("{results}/{ref}/read_counts/filtering.read_counts.tsv")
    params:
        data_colname="filtering"
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

