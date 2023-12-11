#!/usr/bin/env python


## Usage: cat seq.fasta | python approximate_grep.py [n] [queries]


import argparse
import regex # not re!!!
import sys

parser = argparse.ArgumentParser(
                    prog='approximate_grep.py',
                    description='Grep from stdin using regex approximate "fuzzy" matching.',
                    epilog='')

parser.add_argument("queries", nargs='+', help="search queries")
parser.add_argument("-n", type=int, default=0, help="N allowed mismatches. Default: %(default)s")
parser.add_argument("-t", "--type", type=str, default="e", help="Type of mismatch. 'e': all errors, 'i': insertions only, 'd': deletions only, 's': substitutions only. Default: '%(default)s'")
parser.add_argument('-v', '--verbose', action='store_true')
args = parser.parse_args()

pattern = fr"({f'|'.join(args.queries)}){{{args.type}<={args.n}}}"

if args.verbose:
    print("{:#^60}".format(" args "))
    print("queries:", args.queries)
    print(f"n: {args.n}")
    print(f"type: '{args.type}'")
    print(f"regex: '{pattern}'")

if args.verbose:
    print("{:#^60}".format(" grep "))

for line in sys.stdin:
    line = line.rstrip()
    m = regex.findall(pattern, line)
    if m:
        print(line)

        if args.verbose:
            print(f"{args.type}:", m)
            print("-"*60)
