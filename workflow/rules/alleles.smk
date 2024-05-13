## alleles


rule alleles:
    input:
        ##"{results}/{ref}/repeat_region_spanning/no_flanking/{sample}.no_flanking.fastq.gz"
        rules.bbduk_concat_and_fix_strand.output
    output:
        directory("{results}/{ref}/repeat_region_spanning/no_flanking/alleles/{sample}")
    params:
        alleles=lambda wildcards: json.dumps(config["alleles"][wildcards.sample]),
        log="{results}/{ref}/log/alleles.{sample}.log"
    conda:
        "../envs/python.yaml"
    threads:
        2
    shell:
        "mkdir -p {output}; "
        "python workflow/scripts/filter_by_alleles.py "
        "{wildcards.sample} "
        "'{params.alleles}' "
        "{input} "
        "{output} "
        ">{params.log} 2>&1 "


rule alleles_n_alleles_tsv:
    input:
        rules.alleles.output
    output:
        "{results}/{ref}/repeat_region_spanning/no_flanking/alleles/n_alleles/{sample}.n_alleles.tsv"
    shell:
        "( "
        "for f in {input}/*.fastq.gz; do "
        "bname=$(basename $f | cut -d '.' -f2); "
        "n=$(zcat $f | sed -n '1~4p' | wc -l); "
        """echo -e "{wildcards.sample}\t$bname\t$n"; """
        "done "
        ") "
        ##"| sort -k3,3nr "
        ">{output} "
        "&& touch {output} "


rule alleles_motif_positions:
    input:
        rules.alleles.output
    output:
        directory("{results}/{ref}/repeat_region_spanning/no_flanking/alleles/alleles_motif_positions/{sample}")
    log:
        "{results}/{ref}/log/alleles_motif_positions.{sample}.log"
    params:
        motifs=config["motifs"]
    conda:
        "../envs/python.yaml"
    shell:
        "mkdir -p {output}; "
        "( "
        """echo "Sample: {wildcards.sample}"; """
        "for f in {input}/*.fastq.gz; do "
        "bname=$(basename ${{f%.fastq.gz}}); "
        "echo $bname; "
        "python workflow/scripts/waterfall_plot_tsv.py "
        "$f "
        "{params.motifs} "
        "{output}/${{bname}}.motif_position.tsv; "
        "done; "
        "echo; "
        ") "
        ">{log} 2>&1 "


rule alleles_motif_positions_plots:
    input:
        rules.alleles_motif_positions.output
    output:
        "{results}/{ref}/repeat_region_spanning/no_flanking/alleles/alleles_motif_positions_plot/{sample}.motif_positions.png"
    params:
        json_str=lambda wildcards: json.dumps({
            "colors": config["motif_colors"], 
            "subsample": 300,
            "width": 600,
            "x_range": [0, 1800],
            "x_tickvals": list(range(0, 2000, 500)),
            "sample":{
                wildcards.sample:{
                    "alleles": config["alleles"][wildcards.sample]
                }
            }
        }),
    log:
        "{results}/{ref}/log/alleles_motif_positions_plots2.{sample}.log"
    benchmark:
        "{results}/{ref}/.benchmark/alleles_motif_positions_plots2.{sample}.benchmark.tsv"
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/motif_positions_alleles_plot2.py "
        "'{params.json_str}' "
        "{input}/*.motif_position.tsv "
        "{output} "
        ">{log} 2>&1 "

