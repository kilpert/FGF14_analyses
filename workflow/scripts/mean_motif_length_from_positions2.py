#!/usr/bin/env python
# coding: utf-8


import json
import os
import pandas as pd
import plotly.express as px
import sys
from tabulate import tabulate


import os
os.getcwd()


print("{:#^60}".format(f" Args "))
for arg in sys.argv[1:]:
    print(arg)


try:
    config = json.loads(sys.argv[1])
    infiles = sys.argv[2:-1]
    outfile = sys.argv[-1]
except:
    config = None
    infiles = sys.argv[1:-1]
    outfile = sys.argv[-1]

infiles = sorted(infiles) # sort input files

try:
    user_colors_dict = config["colors"]
except:
    user_colors_dict = None


print("{:#^60}".format(f" config "))
print(json.dumps(config, indent=4))


print("{:#^60}".format(f" Variables "))
print("user_colors_dict:", user_colors_dict)
print("infiles:", infiles, len(infiles))
print("outfile:", outfile)


def get_motif_stats(df, sample_name):
    # print("{:#^50}".format(f" {sample_name} "))

    df["bp"] = df["end"] - df["start"]
    # print(df)

    df = df.groupby(["read", "motif"])["bp"].sum().to_frame().reset_index()
    ##print(df)

    # print("{:#^40}".format(f" fix missing (NaN) "))
    ## fix missing motifs by converting to wide ant back to long
    ##print("{:#^40}".format(f" pivot (wide) "))
    df = df.pivot(index="read", columns="motif", values="bp")
    df["read"] = df.index

    ##print("{:#^40}".format(f" melt (long) "))
    df = pd.melt(df, id_vars="read", value_name="bp")

    df = df.fillna(0) # replace NaN by 0!!!
    # print(df)

    # print("{:#^40}".format(f" calculate stats (e.g. mean) "))
    df_mean = df.groupby(["motif"])["bp"].mean().to_frame().reset_index().rename(columns={"bp": "mean_bp"})
    df_median = df.groupby(["motif"])["bp"].median().to_frame().reset_index().rename(columns={"bp": "median_bp"})

    ## merge
    df = pd.merge(df_mean, df_median, on="motif")

    df.mean_bp = df.mean_bp.round(1)
    df["sample"] = sample_name
    df = df.sort_values(by=['motif'], ignore_index=True)  # sort
    df = df.reindex(columns=["sample", "motif", "mean_bp", "median_bp"])
    print(df, "\n")

    return df



df_dict = {}
for infile in infiles:
    sample_name = os.path.basename(infile).split(".")[0]
    print("{:#^60}".format(f" {sample_name} "))

    df = pd.read_csv(infile, sep="\t")
    print(df)
    df_dict[sample_name] = df

print("{:#^60}".format(f" df_dict "))
##print(tabulate(df, headers='keys'))
print(df)


## get all motifs
print("{:#^60}".format(f" motif_list "))
motif_list = []

for sample_name, df in df_dict.items():
    for motif in df["motif"].unique():
        if motif not in motif_list:
            motif_list.append(motif)

print(motif_list, len(motif_list))


## calculate mean
print("{:#^60}".format(f" calculate mean "))

df_motif_length_list = []

for sample_name, df in df_dict.items():
    try:
        df_motif_length_list.append(get_motif_stats(df, sample_name))
    except:
        pass


print("{:#^60}".format(f" Motif length df "))
df = pd.concat(df_motif_length_list, ignore_index=True)
df = df.sort_index(ascending=True)
## print(df)
df["motif_length"] = [len(motif) if motif!="other" else None for motif in df['motif']] # motif length
df["mean_x_motif"] = df["mean_bp"] / df["motif_length"]
df["median_x_motif"] = df["median_bp"] / df["motif_length"]

df = df.round(1)

print(df)
print(df.shape)


## Plots ###############################################################################################################

print("{:#^60}".format(f" Plot "))

figures = []


## colors ##############################################################################################################

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


my_colors = get_colors(sorted(df["motif"].unique()), user_colors_dict=user_colors_dict, return_list=False)
print("my_colors:", my_colors)



## median ################################################################################################################

## grouped

print("{:#^60}".format(f" Median (grouped) "))

fig = px.bar(
    df.fillna(0),
    x="sample",
    y="median_bp",
    color="motif",
    barmode="group",
    color_discrete_map=my_colors,
    hover_name="sample",
    hover_data=["median_x_motif"],
    labels={
        "motif" : "Motif",
        "sample" : "Samples",
        "median_bp" : "Median motif length in reads (bp)",
        "median_x_motif" : "Motifs (x)",
    },
    title=f"Median motif lenghts (grouped)",
    category_orders={"motif": list(my_colors.keys())}
)
fig.update_xaxes(categoryorder='category ascending') # sort x axis
figures.append(fig)



## stacked

print("{:#^60}".format(f" Median (stacked) "))

fig = px.bar(
    df.fillna(0),
    x="sample",
    y="median_bp",
    color="motif",
    barmode="stack",
    color_discrete_map=my_colors,
    hover_name="sample",
    hover_data=["median_x_motif"],
    labels={
        "motif" : "Motif",
        "sample" : "Samples",
        "median_bp" : "Median motif length in reads (bp)",
        "median_x_motif" : "Motifs (x)",
    },
    title=f"Median motif lenghts (stacked)",
    category_orders={"motif": list(my_colors.keys())}
)
fig.update_xaxes(categoryorder='category ascending') # sort x axis
figures.append(fig)


## mean ################################################################################################################

## grouped

print("{:#^60}".format(f" Mean (grouped) "))
print(df)

fig = px.bar(
    df.fillna(0),
    x="sample",
    y="mean_bp",
    color="motif",
    barmode="group",
    color_discrete_map=my_colors,
    hover_name="sample",
    hover_data=["mean_x_motif"],
    labels={
        "motif" : "Motif",
        "sample" : "Samples",
        "mean_bp" : "Mean motif length in reads (bp)",
        "mean_x_motif" : "Motifs (x)",
    },
    title=f"Mean motif lenghts (grouped)",
    category_orders={"motif": list(my_colors.keys())}
)
fig.update_xaxes(categoryorder='category ascending') # sort x axis
figures.append(fig)



## stacked

print("{:#^60}".format(f" Mean (stacked) "))

fig = px.bar(
    df.fillna(0),
    x="sample",
    y="mean_bp",
    color="motif",
    barmode="stack",
    color_discrete_map=my_colors,
    hover_name="sample",
    hover_data=["mean_x_motif"],
    labels={
        "motif" : "Motif",
        "sample" : "Samples",
        "mean_bp" : "Mean motif length in reads (bp)",
        "mean_x_motif" : "Motifs (x)",
    },
    title=f"Mean motif lenghts (stacked)",
    category_orders={"motif": list(my_colors.keys())}
)
fig.update_xaxes(categoryorder='category ascending') # sort x axis
figures.append(fig)


## save plot to html ###################################################################################################

with open(str(outfile), 'w') as f:
    for fig in figures:
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', config={'displaylogo':False}))


os.getcwd()

print(outfile)


## save df to tsv

print("{:#^60}".format(f" Saving df to tsv "))


## long
print("{:#^60}".format(f" df (long) "))
df.to_csv(f"{os.path.splitext(outfile)[0]}.tsv", sep="\t", index=False)
print(df)

# ## wide
# print("{:#^60}".format(f" df (wide) "))
# ##df = df.pivot(index="sample", columns="motif", values="bp")
# df = df.pivot(index="sample", columns="motif")
# df.to_csv(f"{os.path.splitext(outfile)[0]}.wide.tsv", sep="\t", index=True)
# print(df)

