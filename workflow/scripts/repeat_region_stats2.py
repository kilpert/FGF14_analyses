import json
import re
import sys


repeats = sys.argv[1:]
## print("repeats:", repeats, len(repeats))

i = 0

def vprint(a, b):
    print( f"{a:20} {b}" )


reads_dict = {}

while True:
    i += 1

    id = sys.stdin.readline().strip()
    seq = sys.stdin.readline().strip()
    desc = sys.stdin.readline().strip()
    qual = sys.stdin.readline().strip()

    if not id or not seq or not desc or not qual:
        break

    ## print("{:#^60}".format(f" Read {i} "))

    ## vprint("1 - id:", id)
    ## vprint("2 - seq:", seq)
    ## vprint("3 - desc:", desc)
    ## vprint("4 - qual:", qual)

    read_name = id.split()[0][1:]
    ## vprint("read_name:", read_name)
    if not read_name in reads_dict:
        reads_dict[read_name] = {}

    repeat_len = len(seq)
    ## vprint("length:", repeat_len)
    reads_dict[read_name]["length"] = repeat_len

    starts_with_X = seq.startswith(repeats[0])
    ## vprint("starts_with_X:", starts_with_X)
    reads_dict[read_name]["starts_with_X"] = starts_with_X

    reads_dict[read_name]["repeats"] = {}
    for repeat in repeats:
        m = re.findall(repeat, seq)
        reads_dict[read_name]["repeats"][repeat] = len(m)
        ## vprint(f"{repeat}_n:", len(m))


## print("{:#^60}".format(f" reads_dict "))
## print(json.dumps(reads_dict, indent=4))
## print(len(reads_dict))

## print("{:#^60}".format(f" Output "))

# print( reads_dict[list(reads_dict.keys())[0]]["repeats"].keys() )
print( "read_name\tstarts_with_first_motif\trepeat_bp\t{}".format(
        "\t".join([str(x)+"_n" for x in repeats])
    )
)

for read in reads_dict:
    print( "{read}\t{starts_with_X}\t{repeat_bp}\t{reps}".format(
            read=read,
            starts_with_X=reads_dict[read]["starts_with_X"],
            repeat_bp=reads_dict[read]["length"],
            reps="\t".join([str(x) for x in reads_dict[read]["repeats"].values()])
        )
    )
