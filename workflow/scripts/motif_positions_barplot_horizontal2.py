import os
import math
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import random
import sys
import json
from tabulate import tabulate


infile = "data.tsv"
outfile = "data.html"


os.getcwd()


print("{:#^60}".format(f" Args "))
for arg in sys.argv[1:]:
    print(arg)


try:
    config = json.loads(sys.argv[1])
    infile = sys.argv[2]
    subsample_n = int(sys.argv[3])
    outfile = sys.argv[4]
except:
    config = None
    infile = sys.argv[1]
    subsample_n = int(sys.argv[2])
    outfile = sys.argv[3]
## print(config)
## print(infile)
## print(subsample_n)
## print(outfile)


try:
    user_colors_dict = config["colors"]
except:
    user_colors_dict = None


print("{:#^60}".format(f" config "))
print(json.dumps(config, indent=4))


print("{:#^60}".format(f" Variables "))
# bname = os.path.basename(infile).split(".")[0]
sample_name = os.path.basename(infile).split(".")[0]
print("sample_name:", sample_name)
print("outfile:", outfile)


print("{:#^60}".format(" df "))
df = pd.read_csv(infile, sep="\t")
## print(tabulate(df, headers='keys'))
print(df)

if df.shape[0] == 0:
    print("Warning: Empty input file!")
    exit(0)

# print("{:#^60}".format(" n_reads "))
# n_reads = len(df["read"].unique())
# print("reads (unique):", n_reads)

print("{:#^60}".format(" read names (sorted by max length)"))
read_names = df.groupby(["read"])["end"].max().to_frame().sort_values(by=['end'], ascending=False).index.values.tolist()
print("total:", len(read_names))


## subsample df ########################################################################################################
random.seed(123)
subsampled_read_names = random.sample(read_names, min(subsample_n, len(read_names)))
subsampled_read_names = [x for x in read_names if (x in subsampled_read_names)] # reorder according to size
print("subsampled:", len(subsampled_read_names))

print("{:#^60}".format(" df (subsampled) "))
df = df[df['read'].isin(subsampled_read_names)]
##print(tabulate(df, headers='keys'))
print(df)
print(df.shape)

print("{:#^60}".format(" n_reads "))
n_reads = len(df["read"].unique())
print("reads (unique):", n_reads)

print("{:#^60}".format(" max_motifs "))
max_motifs = df.groupby(["read"])["motif"].count().max()
print(max_motifs)
##print(type(max_motifs))
print("isnan:", math.isnan(max_motifs))

print("{:#^60}".format(" df_dict "))
df_dict = {}
for read in subsampled_read_names:
    ##print("{:-^40}".format(" "+read+" "))
    dfx = df.groupby(["read"]).get_group(read).reset_index(drop=True)
    dfx["bp"] = dfx["end"].diff() # relativ value
    dfx.loc[0, "bp"] = dfx.loc[0, "end"]
    dfx["bp"] = dfx["bp"].astype(int)
    df_dict[read] = dfx
    ##print(dfx)

print("{:#^60}".format(" df_dict (grouped) "))
for read,dfx in df_dict.items():
   print("{:-^30}".format(" "+read+" "))
   print(dfx)


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


## Fig

print("{:#^60}".format(" fig "))

figures = []

fig = go.Figure()


print("user_colors_dict:", user_colors_dict)

if not math.isnan(max_motifs):

    for i in range(0, max_motifs):
        # print("{:=^50}".format(" "+str(i)+" "))

        bp = []
        motif = []
        start = []
        end = []

        for read in subsampled_read_names:
            # print("{:-^40}".format(" "+read+" "))

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
            ##print()

        color = [my_colors[x] if x in my_colors else "white" for x in motif]

        ## print("subsampled_read_names:", subsampled_read_names)
        ## print("motif:", motif)
        ## print("bp:", bp)
        ## print("color:", color)
        ## print("start:", start)
        ## print()

        fig.add_trace(go.Bar(
            y=subsampled_read_names,
            x=bp,
            orientation='h',
            name=f"#{i+1}",
            marker_color=color,
            marker_line_width=0,
            ##text=motif,
            ##hoverinfo="skip",
            hovertemplate=
            "<b>%{x}</b><br>" +
            "%{hovertext}",
            ##hovertemplate=None,
            hovertext=[f"{m} ({b} bp)<br>{s}-{e}" for m,b,s,e in zip(motif, bp, start, end)],
            ##width=1,
        ))

fig.update_layout(
    barmode = "stack",
    showlegend = False,
    title = f"{sample_name} - Motif positions<br><sup>{len(subsampled_read_names)} (of {len(read_names)}) reads</sup>",
    yaxis=dict(title='Reads', autorange='reversed', visible=True, showticklabels=False),
    xaxis=dict(title='bp'),
    ##autosize=False,
    width=600,
    height=250+len(df_dict)*2,
    ##height=1000,
    bargap=0,
)

##fig.update_traces(width=1)

##fig.show()

## png
fig.write_image(os.path.splitext(outfile)[0]+".png")
# fig.write_image(outfile)

# ## html
# figures.append(fig)
#
# print(len(figures))
#
# print("{:#^60}".format(" outfile "))
#
# with open(str(outfile), 'w') as f:
#     for fig in figures:
#         f.write(fig.to_html(full_html=True, include_plotlyjs='cdn', config={'displaylogo':False}))

##fig.write_image(f"{bname}.png")
# df.to_csv(outfile, sep="\t", index=False)


print("{:#^60}".format(f" Outfile "))
print(f"{outfile}\n")

