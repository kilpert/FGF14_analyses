## by flanking regions


def read_orientation_upstream(wildcards):
    seq = ""
    if wildcards.orientation[0] == "F":
        seq = Seq(config["flanking_region"]["upstream"])
    if wildcards.orientation[0] == "R":
        seq = Seq(config["flanking_region"]["upstream"]).reverse_complement()
    return seq


def read_orientation_downstream(wildcards):
    seq = ""
    if wildcards.orientation[1] == "F":
        seq = Seq(config["flanking_region"]["downstream"])
    if wildcards.orientation[1] == "R":
        seq = Seq(config["flanking_region"]["downstream"]).reverse_complement()
    return seq


rule read_orientation_by_flanking_regions:
    input:
        "{results}/{ref}/filtering/bbduk/{sample}.fastq.gz"
    output:
        temp("{results}/{ref}/read_orientation/by_flanking_regions/{sample}/{sample}.{orientation}.tsv") # temp
    params:
        flanking_region_upstream=read_orientation_upstream,
        flanking_region_downstream=read_orientation_downstream,
        mismatches=config["repeat_region_spanning"]["agrep_mismatches"]
    benchmark:
        "{results}/{ref}/.benchmark/read_orientation_by_flanking_regions.{sample}.{orientation}.benchmark.tsv"
    conda:
        "../envs/python.yaml"
    threads:
        4
    shell:
        "( "
        """echo -e "Sample\tFlanking Region Orientation\tN"; """
        "N=$(zcat {input} "
        "| sed -n '2~4p' "
        "| python workflow/scripts/approximate_grep.py -n {params.mismatches} {params.flanking_region_upstream} "
        "| python workflow/scripts/approximate_grep.py -n {params.mismatches} {params.flanking_region_downstream} "
        "| wc -l); "
        """echo -e "{wildcards.sample}\t{wildcards.orientation}\t$N"; """
        ") "
        ">{output} "


rule read_orientation_by_flanking_regions_tsv:
    input: 
        sorted(list(expand("{{results}}/{{ref}}/read_orientation/by_flanking_regions/{sample}/{sample}.{orientation}.tsv",
            sample = samples,
            orientation=["FF", "FR", "RF", "RR"],
        )))
    output:
        "{results}/{ref}/read_orientation/by_flanking_regions/flanking_regions_orientation.tsv"
    log:
        "{results}/{ref}/read_orientation/by_flanking_regions/read_orientation_by_flanking_regions.log"
    conda:
        "../envs/R.yaml"
    shell:
        "Rscript workflow/scripts/concatenate_read_counts.R "
        "{input} "
        "{output} "
        ">{log} 2>&1 "


## by motif


rule read_orientation_by_motif:
    input:
        "{results}/{ref}/filtering/bbduk/{sample}.fastq.gz",
    output:
        "{results}/{ref}/read_orientation/by_motif/{sample}/motif_orientation.tsv" # temp
    params:
        motif_forward=config["motifs"][0],
        motif_reverse_complement=Seq(config["motifs"][0]).reverse_complement()
    shell:
        "( "
        """echo -e "Sample\tMotif\tStrand\tCount (N)"; """
        ## +strand (forward)
        "N=$(zcat {input} "
        "| grep -o '{params.motif_forward}' "
        "| wc -l); "
        """echo -e "{wildcards.sample}\t{params.motif_forward}\t+\t$N"; """
        ## -strand (reverse complement)
        "N=$(zcat {input} "
        "| grep -o '{params.motif_reverse_complement}' "
        "| wc -l); "
        """echo -e "{wildcards.sample}\t{params.motif_reverse_complement}\t-\t$N" """
        ") "
        ">{output} "


rule read_orientation_by_motif_tsv:
    input:
        expand("{{results}}/{{ref}}/read_orientation/by_motif/{sample}/motif_orientation.tsv",
            sample=samples,
        )
    output:
        "{results}/{ref}/read_orientation/by_motif/motif_orientation.tsv"
    log:
        "{results}/{ref}/log/read_orientation_by_motif.log"
    conda:
        "../envs/R.yaml"
    shell:
        "Rscript workflow/scripts/concatenate_read_counts.R "
        "{input} "
        "{output} "
        ">{log} 2>&1 "

