rule repeat_region_spanning_flanking_regions_info_txt:
    output:
        "{results}/{ref}/repeat_region_spanning/flanking_regions.info.txt"
    params:
        upstream_F = config["flanking_region"]["upstream"],
        upstream_R = Seq(config["flanking_region"]["upstream"]).reverse_complement(),
        downstream_F = config["flanking_region"]["downstream"],
        downstream_R = Seq(config["flanking_region"]["downstream"]).reverse_complement(),
        mismatches = config["repeat_region_spanning"]["edit_distance"]
    shell:
        "( "
        "echo -e 'upstream: {params.upstream_F}'; "
        "echo -e 'upstream_rc: {params.upstream_R}'; "
        "echo -e 'downstream: {params.downstream_F}'; "
        "echo -e 'downstream_rc: {params.downstream_R}'; "
        "echo -e 'allowed mismatches: {params.mismatches}'; "
        ") "
        ">{output} "


def orientation_upstream(wildcards):
    seq = ""
    if wildcards.orientation[0] == "F":
        seq = Seq(config["flanking_region"]["upstream"])
    if wildcards.orientation[0] == "R":
        seq = Seq(config["flanking_region"]["upstream"]).reverse_complement()
    return seq


def orientation_downstream(wildcards):
    seq = ""
    if wildcards.orientation[1] == "F":
        seq = Seq(config["flanking_region"]["downstream"])
    if wildcards.orientation[1] == "R":
        seq = Seq(config["flanking_region"]["downstream"]).reverse_complement()
    return seq


## all repeat spanning reads
rule filter_spanning_by_orientation:
    input:
        "{results}/{ref}/filtering/bbduk/{sample}.fastq.gz"
    output:
        "{results}/{ref}/repeat_region_spanning/orientation/{orientation}/{sample}.{orientation}.repeat_region_spanning.fastq.gz"
    params:
        ##flanking_region_upstream=config["flanking_region"]["upstream"],
        ##flanking_region_downstream=config["flanking_region"]["downstream"],
        orientation_upstream=orientation_upstream,
        orientation_downstream=orientation_downstream,
        mismatches=config["repeat_region_spanning"]["edit_distance"],
        extra="mm=f rcomp=f",
        ##k=max(len(config["flanking_region"]["upstream"]),len(config["flanking_region"]["downstream"]))
        k_upstream=len(config["flanking_region"]["upstream"]),
        k_downstream=len(config["flanking_region"]["downstream"]),
    log:
        "{results}/{ref}/log/filter_spanning_by_orientation.{sample}.{orientation}.log"
    benchmark:
        "{results}/{ref}/.benchmark/filter_spanning_by_orientation.{sample}.{orientation}.benchmark.tsv"
    conda:
        "../envs/bbmap.yaml"
    threads:
        4
    shell:
        ## upstream
        "bbduk.sh "
        "ordered "
        "{params.extra} "
        "k={params.k_upstream} "
        "edist={params.mismatches} "
        "literal={params.orientation_upstream} "
        "in={input} "
        "outm=stdout.fq "
        "2>{log} "
        ## downstream
        "| bbduk.sh "
        "ordered "
        "{params.extra} "
        "k={params.k_downstream} "
        "edist={params.mismatches} "
        "literal={params.orientation_downstream} "
        "in=stdin.fq "
        "int=f " # no interleaving for unpaired
        "outm={output} "
        "2>>{log} "


rule repeat_region_spanning_by_orientation_read_counts:
    input:
        expand("{{results}}/{{ref}}/repeat_region_spanning/orientation/{orientation}/{sample}.{orientation}.repeat_region_spanning.fastq.gz", 
            orientation=orientations,
            sample=samples,
        )
    output:
        temp("{results}/{ref}/read_counts/repeat_region_spanning_by_orientation.read_counts.tsv") # temp
    shell:
        """echo -e "sample\tprocess\tn" """
        ">{output}; "
        "( "
        "for f in $(printf '%s\n' {input} | sort); do "
        "bname=$(basename $f | cut -d'.' -f1); "
        "orientation=$(basename $f | cut -d'.' -f2); "
        "n=$(zcat $f | sed -n '1~4p' | wc -l); "
        """echo -e "$bname\trepeat_region_spanning_$orientation\t$n"; """
        "done; "
        ") "
        "| sort -k1,1 -k2,2 "
        ">>{output} "


rule orientation_blacklist_txt:
    input:
        FF="{results}/{ref}/repeat_region_spanning/orientation/FF/{sample}.FF.repeat_region_spanning.fastq.gz", 
        RR="{results}/{ref}/repeat_region_spanning/orientation/RR/{sample}.RR.repeat_region_spanning.fastq.gz", 
    output:
        "{results}/{ref}/repeat_region_spanning/orientation/blacklist/{sample}.blacklist.txt"
    shell:
        "comm -12 "
        "<(zcat {input.FF} | sed -n '1~4p' | sort | cut -d' ' -f1 | sed 's/^@//') "
        "<(zcat {input.RR} | sed -n '1~4p' | sort | cut -d' ' -f1 | sed 's/^@//') "
        "| sort "
        ">{output} "


rule bbduk_filtering_and_cutadapt_trimming:
    input:
        fastq="{results}/{ref}/repeat_region_spanning/orientation/{orientation}/{sample}.{orientation}.repeat_region_spanning.fastq.gz",
        blacklist="{results}/{ref}/repeat_region_spanning/orientation/blacklist/{sample}.blacklist.txt"
    output:
        "{results}/{ref}/repeat_region_spanning/orientation/{orientation}/bbduk/{sample}.bbduk.fastq.gz" # temp
    conda:
        "../envs/cutadapt.yaml"
    params:
        flanking_region_upstream=config["flanking_region"]["upstream"],
        flanking_region_downstream=config["flanking_region"]["downstream"],
        to_forward_strand=lambda wildcards: "rcomp=t" if wildcards.orientation=="RR" else "",
    log:
        "{results}/{ref}/log/bbduk_filtering_and_cutadapt_trimming.{sample}.{orientation}.log"
    conda:
        "../envs/cutadapt.yaml"
    threads:
        4
    shell:
        ## 1) reverse complement (to forward strand) if RR
        "reformat.sh "
        "{params.to_forward_strand} "
        "in={input.fastq} "
        "out=stdout.fq "
        "2>{log} "
        ## 2+3) trimming with cutadapt
        "| cutadapt "
        "--discard-untrimmed "
        "-g {params.flanking_region_upstream}...{params.flanking_region_downstream} "
        "- "
        "2>>{log} "
        ## 4) remove ambiguous reads (from blacklist) that have valid flanking regions in FF and RR orientation
        "| filterbyname.sh "
        "include=f " # exclude
        "names={input.blacklist} "
        "in=stdin.fq "
        "int=f " # no interleaving for unpaired
        "out=stdout.fq "
        "2>>{log} "
        ## 5) to fastq (no changes)
        "| reformat.sh "
        "t={threads} "
        "in=stdin.fq "
        "int=f " # no interleaving for unpaired
        "out={output} "
        "2>>{log} "



## In a third step: fastq to fastq (and/or reverse complement!)
## reverse complement of region of interst if on -strand!
rule bbduk_concat_and_fix_strand:
    input:
        sorted(expand("{{results}}/{{ref}}/repeat_region_spanning/orientation/{orientation}/bbduk/{{sample}}.bbduk.fastq.gz",
            orientation=config["orientations"],
        )), # temp
    output:
        "{results}/{ref}/repeat_region_spanning/no_flanking/{sample}.no_flanking.fastq.gz"
    log:
        "{results}/{ref}/log/bbduk_concat_and_fix_strand.{sample}.log"
    params:
        extra="", #"-Xmx1G"
        reverse_complement="rcomp=t" if config["reverse_complement"] else "",
    conda:
        "../envs/bbmap.yaml"
    shell:
        "zcat {input} "
        "| reformat.sh "
        "{params.extra} "
        "{params.reverse_complement} "
        "in=stdin.fq "
        "int=f " # no interleaving for unpaired
        "out={output} "
        ">{log} 2>&1 "


## Create an indicator file in folder if reverse complement
if config["reverse_complement"]:
    rule reverse_complement_indicator_file:
        output:
            "{results}/{ref}/repeat_region_spanning/no_flanking/REVERSE_COMPLEMENT.txt"
        shell:
            "echo 'Sequences in this folder are *reverse complements* of the input sequences!' "
            ">{output} "


## additionally filter for repeats
rule repeat_region_spanning_filter_repeats2:
    input:
        ##"{results}/{ref}/repeat_region_spanning/no_flanking/{sample}.no_flanking.fastq.gz"
        rules.bbduk_concat_and_fix_strand.output
    output:
        "{results}/{ref}/repeat_region_spanning/no_flanking/{filter_repeats}/{sample}.fastq.gz"
    params:
        extra="mm=f rcomp=f",
        k=lambda wildcards: len(wildcards.filter_repeats.replace("_absent", "")),
        seq=lambda wildcards: wildcards.filter_repeats if not wildcards.filter_repeats.endswith("_absent") else wildcards.filter_repeats.replace("_absent", ""),
        out=lambda wildcards: "out" if wildcards.filter_repeats.endswith("_absent") else "outm",
    log:
        "{results}/{ref}/log/repeat_region_spanning_filter_repeats2.{sample}.{filter_repeats}.log"
    conda:
        "../envs/bbmap.yaml"
    threads:
        4
    shell:
        "bbduk.sh "
        "ordered "
        "t={threads} "
        "{params.extra} "
        "k={params.k} "
        "edist=0 "
        "literal={params.seq} "
        "in={input} "
        "{params.out}={output} "
        ">{log} 2>&1 "


rule repeat_region_spanning_repeat_stats:
    input:
        ##"{results}/{ref}/repeat_region_spanning/no_flanking/{filter_repeats}/{sample}.fastq.gz"
        rules.repeat_region_spanning_filter_repeats2.output
    output:
        "{results}/{ref}/repeat_region_spanning/no_flanking/{filter_repeats}/repeat_stats/{sample}.repeat_stats.tsv"
    params:
        motifs=" ".join(config["motifs"])
    shell:
        "zcat {input} "
        "| python workflow/scripts/repeat_region_stats2.py "
        "{params.motifs} "
        ">{output} "


rule repeat_region_spanning_repeat_stats_plots:
    input:
        "{results}/{ref}/repeat_region_spanning/no_flanking/{filter_repeats}/repeat_stats/{sample}.repeat_stats.tsv"
    output:
        "{results}/{ref}/repeat_region_spanning/no_flanking/{filter_repeats}/repeat_stats/html/{sample}.repeat_stats.html"
    params:
        motif_colors=json.dumps({"colors":config["motif_colors"]})
    log:
        "{results}/{ref}/log/repeat_region_spanning_repeat_stats_plots.{filter_repeats}.{sample}.log"
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/repeat_region_stats_plots2.py "
        "'{params.motif_colors}' "
        "{input} "
        "{output} "
        ">{log} 2>&1 "


rule repeat_region_spanning_repeat_stats_density_plot:
    input:
        "{results}/{ref}/repeat_region_spanning/no_flanking/{filter_repeats}/repeat_stats/{sample}.repeat_stats.tsv"
    output:
        "{results}/{ref}/repeat_region_spanning/no_flanking/{filter_repeats}/repeat_stats/density_plot/{sample}.density_plot.png"
    log:
        "{results}/{ref}/log/repeat_region_spanning_repeat_stats_density_plot.{filter_repeats}.{sample}.log"
    conda:
        "../envs/R.yaml"
    shell:
        "Rscript workflow/scripts/density_plot.R "
        "{wildcards.sample} "
        "{input} "
        "{output} "
        ">{log} 2>&1 "
        "|| true "
        "&& touch {output} "


rule repeat_region_spanning_motif_positions_tsv:
    input:
        ##"{results}/{ref}/repeat_region_spanning/no_flanking/{filter_repeats}/{sample}.fastq.gz"
        rules.repeat_region_spanning_filter_repeats2.output
    output:
        "{results}/{ref}/repeat_region_spanning/no_flanking/{filter_repeats}/motif_positions/{sample}.motif_positions.tsv"
    log:
        "{results}/{ref}/log/repeat_region_spanning_motif_positions_tsv.{filter_repeats}.{sample}.log"
    params:
        motifs=config["motifs"]
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/waterfall_plot_tsv.py "
        "{input} "
        "{params.motifs} "
        "{output} "
        ">{log} 2>&1 "


rule repeat_region_spanning_motif_positions_barplot_horizontal:
    input:
        "{results}/{ref}/repeat_region_spanning/no_flanking/{filter_repeats}/motif_positions/{sample}.motif_positions.tsv"
    output:
        "{results}/{ref}/repeat_region_spanning/no_flanking/{filter_repeats}/motif_positions/png_horizontal/{sample}.motif_positions_barplot.png"
    params:
        subsample_n=500,
        motif_colors=json.dumps({"colors":config["motif_colors"]})
    log:
        "{results}/{ref}/log/repeat_region_spanning_motif_positions_barplot_horizontal.{filter_repeats}.{sample}.log"
    benchmark:
        "{results}/{ref}/.benchmark/repeat_region_spanning_motif_positions_barplot_horizontal.{filter_repeats}.{sample}.benchmark.tsv"
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/motif_positions_barplot_horizontal2.py "
        "'{params.motif_colors}' "
        "{input} "
        "{params.subsample_n} "
        "{output} "
        ">{log} 2>&1 "
        "&& touch {output} "


rule repeat_region_spanning_motif_positions_mean_motif_length:
    input:
        expand("{{results}}/{{ref}}/repeat_region_spanning/no_flanking/{{filter_repeats}}/motif_positions/{sample}.motif_positions.tsv",
            sample=samples
        )
    output:
        "{results}/{ref}/repeat_region_spanning/no_flanking/{filter_repeats}/motif_positions/mean_motif_length_from_position/mean_motif_length_from_position.html"
    params:
        motif_colors=json.dumps({"colors":config["motif_colors"]})
    log:
        "{results}/{ref}/log/mean_motif_length_from_position.{filter_repeats}.log"
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/mean_motif_length_from_positions2.py "
        "'{params.motif_colors}' "
        "{input} "
        "{output} "
        ">{log} 2>&1 "

