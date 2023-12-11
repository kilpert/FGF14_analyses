#!/usr/bin/env python
# coding: utf-8


import os
import pandas as pd
import plotly.express as px
import sys
import json
from tabulate import tabulate


os.getcwd()


print("{:#^60}".format(f" Args "))
for arg in sys.argv[1:]:
    print(arg)


try:
    config = json.loads(sys.argv[1])
    infile = sys.argv[2]
    outfile = sys.argv[3]
except:
    config = None
    infile = sys.argv[1]
    outfile = sys.argv[2]
## print(config)
## print(infile)
## print(outfile)


try:
    user_colors_dict = config["colors"]
except:
    user_colors_dict = None


print("{:#^60}".format(f" config "))
print(json.dumps(config, indent=4))


print("{:#^60}".format(f" Variables "))
sample_name = os.path.basename(infile).split(".")[0]
print("sample_name:", sample_name)


print("{:#^60}".format(f" df "))
df = pd.read_csv(infile, sep="\t")
df = df.sort_values(by=["repeat_bp"], ascending=False)
##print(tabulate(df, headers='keys'))


## Calculate motif lenghts in bp
motif_bp_list = []

## add columns with length of motif in bp
df_colnames = df.columns.values.tolist()
for motif in [motif for motif in df_colnames if motif.endswith("_n")]:
    motif_seq = motif.replace("_n", "")
    ## print(motif, motif_seq, len(motif_seq))

    motif_bp = f"{motif_seq}_bp"
    df[motif_bp] = df[motif] * len(motif_seq)
    motif_bp_list.append(motif_bp)

## other_bp
df["other_bp"] = df["repeat_bp"]
for motif_bp in motif_bp_list:
    df["other_bp"] -= df[motif_bp]
print(tabulate(df, headers='keys'))



## Figures
figures = [] # colect figures


## Histogram of repeat regions lengths
print("{:#^60}".format(f" Histogram "))

fig = px.histogram(df,
    x="repeat_bp",
    nbins=100,
    color="starts_with_first_motif",
    color_discrete_map={
       True: 'cadetblue',
       False: 'darksalmon'
    },
    marginal="rug",
    title=f"{sample_name} - Histogram of repeat lenghts (bp)<br><sup>n_reads={len(df)}, median_bp={df['repeat_bp'].median()}</sup>"
)
##fig.show()
figures.append(fig)


## Repeat region compositon bar plot
print("{:#^60}".format(f" Bar plot "))

df_bp = df[['read_name'] + motif_bp_list + ['other_bp']]

df_bp_long = df_bp.melt(id_vars=["read_name"], var_name="Motif", value_name="bp")
df_bp_long["Motif"] = df_bp_long["Motif"].str.replace("_bp", "") # remove "_bp" from variable
print(tabulate(df_bp_long, headers='keys'))



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


my_colors = get_colors(df_bp_long["Motif"].unique(), user_colors_dict=user_colors_dict)
print("my_colors:", my_colors)

fig = px.bar(df_bp_long,
             x="read_name",
             y="bp",
             color="Motif",
             color_discrete_sequence=my_colors,
             ## color_discrete_map={
             ##    'AAG': 'green',
             ##    'TACC': 'blue',
             ##    'other': 'black'
             ## },
             title=f"{sample_name} - Repeat region composition<br><sup>n_reads={len(df)}, median_bp={df['repeat_bp'].median()}</sup>",
)
## fig.update_layout(
##     xaxis_title="Read name"
## )
##fig.show()
figures.append(fig)


## Write all figures to html
with open(str(outfile), 'w') as f:
    for fig in figures:
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', config={'displaylogo':False}))


print("{:#^60}".format(f" Outfile "))
print(f"{outfile}\n")

