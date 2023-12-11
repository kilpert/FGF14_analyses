import gzip
import pandas as pd
import re
import sys


## functions ###########################################################################################################

def get_motif_df(seq, motifs):
    motif = []
    start = []
    end = []

    for query in query_motifs:
        for m in re.finditer(query, seq):
            motif.append(query)
            start.append(m.start())
            end.append(m.end())

    df = pd.DataFrame({"motif": motif, "start": start, "end": end})
    return df.sort_values(by=["start"]).reset_index(drop=True)


def show_motifs(seq, df):
    try:
        print("{:6} {}".format("seq", seq))
        for index, row in df.iterrows():
            upstream = "."*row["start"]
            downstream = "."*(len(seq)-row["end"])
            print("{:6} {}".format(row["motif"], upstream + "*"*(row["end"]-row["start"]) + downstream))
    except:
        pass


def dist_bp(a, b):
    return min(b) - max(a)


def adjacent(tuple_list):
    bools = []
    last = None

    for e in tuple_list:
        if last == None:
            last = e
            bools.append(0)
            continue

        bools.append(dist_bp(last, e))
        last = e

    return bools


def dist_to_prev(tuple_list):
    bools = []
    last = None

    for e in tuple_list:
        if last == None:
            last = e
            bools.append(0)
            continue

        bools.append((min(e) - max(last)) * -1)
        last = e

    return bools


def dist_to_next(tuple_list):
    bools = []
    last = None

    for e in tuple_list:
        if last == None:
            last = e
            # bools.append(0)
            continue

        bools.append(dist_bp(last, e))
        last = e

    bools.append(None)

    return bools


def join_motifs_df(df):
    df["motif_groups"] = df['motif'].ne(df['motif'].shift()).cumsum()
    df["dist_to_prev_bp"] = dist_to_prev(list(zip(df.start, df.end))) # dist to prev mofif
    df["dist_groups"] = (df['dist_to_prev_bp']!=0).cumsum()
    ##print(df)
    return df.groupby(["motif_groups", "dist_groups"]).agg({'motif':'first', 'start':'min', 'end':'max'}).reset_index(drop=True)


def fix_overlaps(df):
    df["dist_to_next_bp"] = dist_to_next(list(zip(df.start, df.end)))
    df.loc[df["dist_to_next_bp"]<0, ["end"]] = df[df["dist_to_next_bp"]<0].apply(lambda x: int(x.end + x.dist_to_next_bp), axis=1)
    df = df.drop(columns=['dist_to_next_bp'])
    ##print(df)
    return df


def fill_gaps(df, length):
    last = None
    motif = []
    start = []
    end = []

    tuple_list = list(zip(df.start, df.end))

    ## upstream first motif
    if tuple_list[0][0] > 0:
        motif.append("other")
        start.append(0)
        end.append(tuple_list[0][0])

    ## within motifs
    for x in tuple_list:
        if last == None:
            last = x
            continue

        ##print(last, x)
        if min(x) - max(last) > 0:
            ##print(max(last), min(x))
            motif.append("other")
            start.append(max(last))
            end.append(min(x))
        last = x

    ## downstream first motif
    tuple_list = list(zip(df.start, df.end))
    motif.append("other")
    start.append(max(x))
    end.append(length)

    df2 = pd.DataFrame({"motif": motif, "start": start, "end": end})

    df = pd.concat([df, df2])
    df = df.sort_values(by=['start', 'end'])
    ##print(df)
    return df


## main ################################################################################################################

if __name__ == "__main__":

    print("{:#^60}".format(" Argv "))
    print(sys.argv)

    infile = sys.argv[1]
    print("infile:", infile)

    query_motifs = sys.argv[2:-1]
    print("query_motifs:", query_motifs)

    outfile = sys.argv[-1]
    print("outfile:", outfile)

    print("{:#^60}".format(" Seq "))

    df_list = []

    with gzip.open(infile, 'rb') as f:
        i = 0
        while True:
            identifier = f.readline().decode('utf-8').strip()
            if not identifier:
                break
            if not identifier.startswith("@"):
                break

            # if identifier == "@55657bf8-6170-4039-8526-2fb731a253f6":

            id = identifier.split()[0].replace("@", "")
            ##print("{:#^60}".format(" "+id+" "))

            seq = f.readline().decode('utf-8').strip()
            ##print(seq)

            desc = f.readline()
            qual = f.readline()

            ## df

            # if id != "55657bf8-6170-4039-8526-2fb731a253f6":
            #     continue

            ## print("{:#^60}".format(" All motifs "))
            df = get_motif_df(seq, query_motifs)

            if df.shape[0] == 0: # if there are no motifs in read!!!
                continue
            ##show_motifs(seq, df)

            ## print("{:#^60}".format(" Join motifs "))
            df = join_motifs_df(df)
            ##print(df)
            ##show_motifs(seq, df)

            ## print("{:#^60}".format(" Fix overlaps "))
            df = fix_overlaps(df)
            ##print(df)
            ##show_motifs(seq, df)

            ## print("{:#^60}".format(" Fill gaps "))
            df = fill_gaps(df, len(seq))
            ##print(df)
            ##show_motifs(seq, df)

            df.insert(0, 'read', id)
            ##print(df)
            df_list.append(df.reset_index(drop=True))

            # if id == "55657bf8-6170-4039-8526-2fb731a253f6":
            #     break



    ## for i in range(0, len(df_list)):
    ##     print("{:#^60}".format(f" {i} "))
    ##     print(df_list[i])

    if df_list:
        df = pd.concat(df_list).reset_index(drop=True)
    else:
        df = pd.DataFrame({"read": [], "motif": [], "start": [], "end": []})

    print("{:#^60}".format(" df "))
    print(df)
    print("{:#^60}".format(" df info "))
    print(df.info())

    ## problematic reads
    print("{:#^60}".format(" start NA "))
    print( "start any NA:", df["start"].isna().any() )
    print( df[df["start"].isna()] )
    print("{:#^60}".format(" end NA "))
    print( "end any NA:", df["end"].isna().any() )
    print( df[df["end"].isna()] )

    ## some info
    ## print("{:#^60}".format(" head "))
    ## print(df.head(10))
    ## print("{:#^60}".format(" tail "))
    ## print(df.tail(10))
    ## print(df.shape)
    ## print("n reads:", len(df["read"].unique()))

    df.to_csv(outfile, sep="\t", index=False)

