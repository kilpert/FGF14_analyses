import os
import math
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pandas as pd
import random
import sys
import json
from tabulate import tabulate
import re


infile = "data.tsv"
outfile = "data.html"


print("{:#^60}".format(f" Args "))
for arg in sys.argv[1:]:
    print(arg)


try:
    config = json.loads(sys.argv[1])
    infiles = sys.argv[2:-1]
    outfile = sys.argv[-1]
except:
    print("Error! Invalid input!")
    exit(1)

print("{:#^60}".format(f" Config "))
print(json.dumps(config, indent=4))


## config to json (for paper)
outfile_json = f"{os.path.splitext(outfile)[0]}.json"
with open(outfile_json, "w") as f:
    f.write(json.dumps(config))


try:
    user_colors_dict = config["colors"]
except:
    user_colors_dict = None

try:
    sample_name = list(config["sample"].keys())[0]
except:
    sample_name = ""

try:
    d_alleles = config["sample"][sample_name]["alleles"]
except:
    d_alleles = {}

try:
    width = config["width"]
except:
    width = None

try:
    x_range = config["x_range"]
except:
    x_range = None

try:
    x_tickvals = config["x_tickvals"]
except:
    x_tickvals = None

try:
    y_tickvals = config["y_tickvals"]
except:
    y_tickvals = None


## Path from json_str
def json_str_to_allele_str(json_str):
    allele_str = json_str
    allele_str = re.sub('\s', '', allele_str)
    allele_str = re.sub(":{", "", allele_str)
    allele_str = re.sub("[{}]", "", allele_str)
    allele_str = re.sub('"{2,}', ' ', allele_str)
    allele_str = re.sub('"', '', allele_str)
    allele_str = re.sub('include ', '', allele_str)
    ## print(f"allele path: '{allele_str}'")
    return allele_str


## Collect info for fig
d_figinfo = {
    "alleles": {},
    "total": {}
}


## allele info str
for allele in d_alleles.keys():
    info = json_str_to_allele_str(json.dumps(d_alleles[allele]))
    ##print("info:", info)
    d_figinfo["alleles"][allele] = {"info": info}
d_figinfo["alleles"]["other"] = {"info": "other"}


print("{:#^60}".format(f" Variables "))
print("sample_name:", sample_name)
print("infiles:", infiles, len(infiles))
print("outfile:", outfile)
print("d_alleles:", json.dumps(d_alleles, indent=4))
print("d_figinfo:", json.dumps(d_figinfo, indent=4))


## Get motif position for each allele

def get_dict_of_df(infiles):
    print("{:#^60}".format(f" get_dict_of_df "))

    d_df = {}

    for infile in infiles:
        print("infile:", infile)

        bname = os.path.basename(infile)
        ## print("{:=^60}".format(f" {bname} "))

        ## get sample and allele from tsv filename
        m = re.search("(\w+)\.(\w+)\.motif_position\.tsv", bname)
        sample = m.group(1)
        allele = m.group(2)
        ## print(sample, allele)

        ##print("{:-^60}".format(" df "))
        df = pd.read_csv(infile, sep="\t")

        ## NOT include empty data frames (files)
        if len(df.index) > 0:
            d_df[allele] = df
        ## print(tabulate(df, headers='keys'))
        ## print(df)

        n_reads = len(df["read"].unique())
        print("reads:", n_reads)

    return d_df


d_df = get_dict_of_df(infiles)
print("{:#^60}".format(f" d_df "))
print(d_df, d_df.keys())


## n_reads (total and for each allele)
n_total = 0
for allele, df in d_df.items():
    n_reads = len(df["read"].unique())
    n_total += n_reads
    d_figinfo["alleles"][allele]["n_reads"] = n_reads
d_figinfo["total"]["n_reads"] = n_total


## allele_perc
for allele, df in d_df.items():
    n_reads = d_figinfo["alleles"][allele]["n_reads"]
    allele_perc = round(n_reads / n_total * 100, 2)
    d_figinfo["alleles"][allele]["allele_perc"] = allele_perc


## df_all
## df_all = pd.concat([x[1] for x in d_df.items()])


## Subsample

def get_dict_of_df_subsampled(d_df, n_subsample=100, seed=123):
    print("{:#^60}".format(f" get_dict_of_df_subsampled "))

    print(f"n_subsample: {n_subsample}")

    d_df_subsampled = {}

    for allele, df in d_df.items():
        print("{:=^60}".format(f" {allele} "))
        print(df)
        ## read names (sorted by max read length)
        ##read_names = df.groupby(["read"])["end"].max().to_frame().sort_values(by=['end'], ascending=False) #.index.values.tolist()
        read_names = df["read"].unique().tolist()
        ##print(read_names)
        print("n_reads:", len(read_names))

        ## n_subsample = 200
        ## if "subsample_total" in config:
        ##     print("subsample_total")
        ##     n_subsample = round(n_reads * config["subsample_total"] / n_total)
        ## elif "subsample" in config:
        ##     print("subsample")
        ##     n_subsample = config["subsample"]
        ## print("n_subsample:", n_subsample, type(n_subsample))

        ## subsample df (by read name) ########################################################################################################
        random.seed(seed)
        subsampled_ordered_read_names = random.sample(read_names, min(n_subsample, len(read_names)))
        subsampled_ordered_read_names = [x for x in read_names if (x in subsampled_ordered_read_names)] # reorder according to size
        n_subsampled_reads = len(subsampled_ordered_read_names)
        ##print("n subsampled read names:", n_subsampled_reads)

        n_reads = d_figinfo["alleles"][allele]["n_reads"]
        d_figinfo["alleles"][allele]["n_subsampled"] = n_subsampled_reads
        d_figinfo["alleles"][allele]["allele_subsampled_perc"] = round(n_subsampled_reads * 100 / n_reads, 2)

        print("{:-^60}".format(" df (subsampled) "))
        df = df[df['read'].isin(subsampled_ordered_read_names)]
        ##print(tabulate(df, headers='keys'))
        print(df)
        print("n_reads:", len(df["read"].unique().tolist()))
        ##print(df.shape)

        ## outfile_tsv (for paper)
        outfile_tsv = f"{os.path.join(os.path.dirname(outfile), os.path.basename(outfile).split('.')[0])}.{allele}.motif_positions.tsv"
        print(outfile_tsv)
        df.to_csv(outfile_tsv, sep="\t", index=False)

        d_df_subsampled[allele] = df

    return d_df_subsampled


d_df_subsampled = get_dict_of_df_subsampled(d_df, n_subsample=config["subsample"])
print("{:#^60}".format(f" d_df_subsampled "))
print(d_df_subsampled)


## get read-size-ordered read_names (longest to shortest)

def get_read_size_ordered_ordered_read_names(d_df):
    print("{:#^60}".format(f" get_read_size_ordered_ordered_read_names "))
    
    d_ordered_read_names = {}

    for allele, df in d_df.items():
        print("{:=^60}".format(f" {allele} "))

        ## size ordered read names
        ##read_names = df.groupby(["read"])["end"].max().to_frame().sort_values(by=['end'], ascending=False).index.values.tolist()
        df_size = df.groupby(["read"])["end"].max().to_frame().sort_values(by=['end'], ascending=False)
        print(df_size.head())
        read_names = df_size.index.values.tolist()
        d_ordered_read_names[allele] = read_names

    return d_ordered_read_names
    

d_ordered_read_names = get_read_size_ordered_ordered_read_names(d_df_subsampled)
print("{:#^60}".format(f" d_ordered_read_names "))
##print(d_ordered_read_names)


## get max_motifs

def get_max_motifs(d_df):
    print("{:#^60}".format(f" get_max_motifs(d_df) "))

    d_max_motifs = {}

    for allele, df in d_df.items():
        print("{:=^60}".format(f" {allele} "))

        print("motifs:")
        ## for motif in df.groupby(["read"])["motif"]:
        ##     print(motif)
        ## print(df.groupby(["read"])["motif"].count())

        max_motifs = df.groupby(["read"])["motif"].count().max()
        ## print("max_motifs:")
        ## print(max_motifs)
        ## print(type(max_motifs))
        print("isnan:", math.isnan(max_motifs))
        d_max_motifs[allele] = max_motifs

    return d_max_motifs


d_max_motifs = get_max_motifs(d_df_subsampled)
print("{:#^60}".format(f" d_max_motifs "))
##print(d_max_motifs)


## get d_df_dict

def get_d_df_dict(d_df):
    print("{:#^60}".format(f" get_d_df_dict "))

    d_df_dict = {}

    for allele, df in d_df.items():
        ## print("{:=^60}".format(f" {allele} "))

        ##print("{:-^60}".format(" df_dict "))
        df_dict = {}

        for read in d_ordered_read_names[allele]:
            print("{:-^50}".format(f" {read} "))
            dfx = df.groupby(["read"]).get_group(read).reset_index(drop=True)
            dfx["bp"] = dfx["end"].diff() # relativ value
            dfx.loc[0, "bp"] = dfx.loc[0, "end"]
            dfx["bp"] = dfx["bp"].astype(int)
            ##print(dfx)
            df_dict[read] = dfx

        d_df_dict[allele] = df_dict
            
    return d_df_dict


d_df_dict = get_d_df_dict(d_df_subsampled)
print("{:#^60}".format(f" d_df_dict "))

for allele, df in d_df_dict.items():
    print("{:=^60}".format(" df "))
    ##print(df)


## remove alleles of empty data frames from d_figinfo
d_figinfo["alleles"] = {allele: d_figinfo["alleles"][allele] for allele in list(d_df_subsampled.keys())}


print("{:#^60}".format(f" d_figinfo "))
print(json.dumps(d_figinfo, indent=4))


## generate plots for each allele ########################################################################################

print("{:#^60}".format(f" Fig "))


## colors

def get_colors(varnames, default_colors_list=px.colors.qualitative.G10, user_colors_dict=None, return_list=True):
    color_dict = {}

    ## default colors
    for i in range(0, len(varnames)):
        color_dict[varnames[i]] =  default_colors_list[i]

    ## overwrite colors
    if user_colors_dict:
        for varname in color_dict.keys():
            if varname in user_colors_dict.keys():
                color_dict[varname] = user_colors_dict[varname]

    if return_list:
        return(list(color_dict.values()))
    else:
        return(color_dict)


## print( [a for a in list(d_figinfo["alleles"].keys())] )

## subplot titles
subplot_titles = []
for allele in d_figinfo["alleles"].keys():
    info = d_figinfo['alleles'][allele]['info']
    n_subsampled = d_figinfo['alleles'][allele]['n_subsampled']
    allele_subsampled_perc = d_figinfo['alleles'][allele]['allele_subsampled_perc']
    n_reads = d_figinfo['alleles'][allele]['n_reads']
    allele_perc = d_figinfo['alleles'][allele]['allele_perc']
    subplot_title = f"{info}<br><sup>{n_subsampled} ({allele_subsampled_perc}%) of {n_reads} ({allele_perc}%) reads</sup>"
    print(subplot_title)
    subplot_titles.append(subplot_title)


## Figure

fig = make_subplots(
    rows=1,
    cols=len(d_figinfo["alleles"].keys()),
    ##subplot_titles = list(d_figinfo["alleles"].keys()),
    subplot_titles = subplot_titles,
    shared_xaxes = "rows",
    horizontal_spacing = 0.05
)


print("{:-^60}".format(f" my_colors "))
df_subsampled_all = pd.concat([x[1] for x in d_df_subsampled.items()])
## print(df_subsampled_all)
my_colors = get_colors(sorted(df_subsampled_all["motif"].unique()), user_colors_dict=user_colors_dict, return_list=False)
## print("my_colors:", my_colors)
## for motif, color in my_colors.items():
##     print(f"{motif}: {color}")


current_col = 0


for allele, df in d_df_subsampled.items():
    print("{:=^60}".format(f" {allele} "))
    ##print(df)

    current_col += 1
    print(f"current_col: {current_col}")

    ordered_read_names = list(d_df_dict[allele].keys())
    ## print("ordered_read_names:", ordered_read_names, len(ordered_read_names))

    max_motifs = d_max_motifs[allele]
    print(f"max_motifs: {max_motifs}")

    if not math.isnan(max_motifs):
    
        for i in range(0, max_motifs):
            print("{:-^50}".format(f" {i+1} {round((i+1)/max_motifs*100, 1)}% ")+"\r")

            bp = []
            motif = []
            start = []
            end = []

            df_dict = d_df_dict[allele]
            ## print(df_dict)
            ## print(type(df_dict))

            
            for read in ordered_read_names:
                ## print("{:.^50}".format(f" {read} "))

                try:
                    bp.append(df_dict[read].iloc[i]["bp"])
                    motif.append(df_dict[read].iloc[i]["motif"])
                    start.append(df_dict[read].iloc[i]["start"])
                    end.append(df_dict[read].iloc[i]["end"])
                except:
                    bp.append(None)
                    motif.append(None)
                    start.append(None)
                    end.append(None)

            color = [my_colors[x] if x in my_colors else "white" for x in motif]

            ## print("ordered_read_names:", ordered_read_names, len(ordered_read_names), type(ordered_read_names))
            ## print("motif:", motif)
            ## print("bp:", bp)
            ## print("color:", color)
            ## print("start:", start)
            ## print()

            fig.add_trace(go.Bar(
                y=ordered_read_names,
                x=bp,
                name=f"#{i+1}",
                orientation='h',
                marker_color=color,
                marker_line_width=0,
                ##text=f"#{i+1}",
                ##hoverinfo="skip",
                ## hovertemplate=
                ## "<b>%{x}</b><br>" +
                ## "%{hovertext}",
                ##hovertemplate=None,
                ##hovertext=[f"{m} ({b} bp)<br>{s}-{e}" for m,b,s,e in zip(motif, bp, start, end)],
                ##width=1,
            ), row=1, col=current_col)

            fig.update_xaxes(title="bp", row=1, col=current_col)
            fig.update_yaxes(categoryorder='array', categoryarray=list(reversed(ordered_read_names)), row=1, col=current_col)


fig.update_layout(
    barmode = "stack",
    showlegend = False,
    yaxis=dict(title='Reads'),
    ##autosize=False,
    bargap=0,
    title = f"{sample_name}<br><sup>{d_figinfo['total']['n_reads']} reads (total)</sup>",
    title_x=0.5,
    width=width*len(d_figinfo["alleles"].keys()),
    height=200+len(df_dict)*3,
    ##height=1000,
)

if "x_range" in config:
    fig.update_xaxes(range=[config["x_range"][0], config["x_range"][1]])

if "x_tickvals" in config:
    print("x_tickvals:", x_tickvals)
    fig.update_xaxes(
        tickmode = 'array',
        tickvals = x_tickvals
    )
    
fig.update_traces(width=1)
fig.update_yaxes(showticklabels=False)
fig.update_layout(margin=dict(t=130))
fig.show()

## png
fig.write_image(outfile)
print("{:#^60}".format(f" Outfile "))
print(f"{outfile}\n")

