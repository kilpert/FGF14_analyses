import gzip
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import statistics
import sys

print("{:#^60}".format(" Argv "))
print(sys.argv)


print("{:#^60}".format(" Variables "))
infile = sys.argv[1]

reads = sys.argv[2:-1]

outfile = sys.argv[-1]

n_reads = len(reads)
# try:
if len(reads)==1:
    try:
        n_reads = int(reads[0])
        reads = []
    except:
        pass

if not reads and n_reads==0:
    n_reads = 3
    reads = []

print(f"infile: {infile}")
print(f"n_reads: {n_reads}")
print(f"reads: {reads}")
print(f"outfile: {outfile}")

sample = os.path.basename(infile).split(".")[0]
print(f"sample: {sample}")


seq_dict = {}

read_names = []

print("{:#^60}".format(" Reads "))
with gzip.open(infile, 'rt') as f:
    i = 0
    while True:
        identifier = f.readline().strip()
        if not identifier:
            break
        if not identifier.startswith("@"):
            break

        # if identifier == "@55657bf8-6170-4039-8526-2fb731a253f6":

        id = identifier.split()[0].replace("@", "")
        ##print("{:#^60}".format(" "+id+" "))

        seq = f.readline().strip()
        ##print(seq)

        desc = f.readline().strip()
        qual = f.readline().strip()

        if not reads and len(read_names) < n_reads:
            read_names.append(id)

            if id in read_names:
                print(id)
                # print(seq)
                # print(desc)
                # print(qual)

                if id not in seq_dict:
                    seq_dict[id] = {}
                seq_dict[id]["seq"] = seq
                seq_dict[id]["desc"] = desc
                seq_dict[id]["qual"] = qual
        else:
            if id in reads:
                print(id)
                # print(seq)
                # print(desc)
                # print(qual)

                if id not in seq_dict:
                    seq_dict[id] = {}
                seq_dict[id]["seq"] = seq
                seq_dict[id]["desc"] = desc
                seq_dict[id]["qual"] = qual


print("{:#^60}".format(f" seq_dict "))
print(len(seq_dict))


print("{:#^60}".format(f" plots "))
figures = []

for k, v in seq_dict.items():
    print("{:#^60}".format(f" {k} "))

    read_name = k
    seq = seq_dict[k]["seq"]
    qual = seq_dict[k]["qual"]

    ## plot data frame
    d = {
        "seq": [x for x in seq],
        "qual": [x for x in qual],
    }
    df = pd.DataFrame(d)

    df["phred"] = [ord(x) - 33 for x in df["qual"]]

    df = df.drop(columns=["qual"])
    df['bp'] = range(1, 1 + len(df))

    df["floor"] = -2

    print( df )

    phred_mean = statistics.mean(df["phred"])
    phred_max = max(df["phred"])
    print( "phred mean:", phred_mean )
    print( "phred max:", phred_max )

    fig1 = px.line(df, x="bp", y="phred")
    fig1.update_traces(line_color='grey')

    fig2 = px.scatter(df, x="bp", y="phred", color="seq", text="seq",
                      color_discrete_map={
                          'A': 'green',
                          'T': 'red',
                          'C': 'blue',
                          'G': 'orange'
                      }
                      )
    fig2.update_traces(
        showlegend=False,
        textfont_color="white",
        textfont_size=12,
        textfont_family="Arial",
        marker=dict(
            size=13,
            ##opacity=0,
            line=dict(width=0,
                      color='white')),
        selector=dict(mode='markers+text'))

    fig3 = px.scatter(df, x="bp", y="floor", color="seq", text="seq",
                      color_discrete_map={
                          'A': 'green',
                          'T': 'red',
                          'C': 'blue',
                          'G': 'orange'
                      }
                      )
    fig3.for_each_trace(lambda t: t.update(textfont_color=t.marker.color))
    fig3.update_traces(
        showlegend=False,
        ##textfont_color="black",
        textfont_size=12,
        textfont_family="Arial",
        marker=dict(
            size=13,
            opacity=0,
            line=dict(width=0,
                      color='white')),
        selector=dict(mode='markers+text'))

    fig = go.Figure(data=fig1.data + fig2.data + fig3.data)
    ##fig.add_hline(y=phred_mean, line_dash="dash", annotation_text=f"Mean: {phred_mean:.1f}", annotation_position="top left")
    fig.add_hline(y=phred_mean, line_dash="dash", opacity=0.2)
    fig.update_layout(
        xaxis_range=[0, 100],
        yaxis_range=[-5, max(50, phred_max+phred_max*0.02)],
        title=f"{sample} - {read_name}<br><br><sup>Read length: {len(df)}, Phred mean: {phred_mean:.1f}</sup>",
        xaxis_title="Read sequence",
        yaxis_title="Phred score",
        yaxis_fixedrange=True,  # disallow vertical pan
        dragmode="pan"  # default drag mode
    )

    fig.show()

    figures.append(fig)

##print(figures)
print(len(figures))

print(f"outfile: {outfile}")

with open(str(outfile), 'w') as f:
    for fig in figures:
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', config={'displaylogo':False}))

