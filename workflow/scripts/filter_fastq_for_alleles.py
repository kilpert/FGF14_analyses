#!/usr/bin/env python

import argparse
import sys
import yaml
import json
import re
import gzip


## Usage:
## python filter_fastq_for_alleles.py --yaml-file alleles.yaml > >(gzip >include.fastq.gz) 2> >(gzip >exclude.fastq.gz)


parser = argparse.ArgumentParser(
    prog='filter_fastq_for.py',
    description='Filter FASTQ file for exact sequence and min/max lenghts',
    epilog='')
parser.add_argument("--json-str", dest="json_str", type=str, default=None, help="txt string (json)). Default: %(default)s")
parser.add_argument("--json-file", dest="json_file", type=str, default="", help="allele file (json). Default: %(default)s")
parser.add_argument("--yaml-file", dest="yaml_file", type=str, default="", help="allele file (yaml). Default: %(default)s")
## parser.add_argument("--out", dest="include_fastq", type=str, default="", help="output of included reads (fastq.gz). Default: '%(default)s'")
## parser.add_argument("--out2", dest="exclude_fastq", type=str, default="", help="output of excluded reaas (fastq.gz). Default: '%(default)s'")
parser.add_argument('-v', '--verbose', action='store_true')
args = parser.parse_args()


if args.verbose:
    print("{:#^60}".format(" Args "))
    print(args)


if args.json_str:
    json_str = args.json_str
    d = json.loads(json_str)

elif args.json_file:
    with open(args.json_file, 'r') as f:
        d = json.load(f)

elif args.yaml_file:
    with open(args.yaml_file, 'r') as f:
        d = yaml.safe_load(f)

else:
    print("No allele(s) specified. Exit.")
    exit(0)


if args.verbose:
    print("{:#^60}".format(" Allels (dict) "))
    print(json.dumps(d, indent=4))


## Variables
try:
    minlen = d["read"]["minlen"]
except:
    minlen = None

try:
    maxlen = d["read"]["maxlen"]
except:
    maxlen = None

try:
    n_include = len(d["include"].keys())
except:
    n_include = 0

try:
    n_exclude = len(d["exclude"].keys())
except:
    n_exclude = 0

if args.verbose:
    print("{:#^60}".format(f" Variables "))

    if minlen:
        print(f"minlen: {minlen}")

    if maxlen:
        print(f"maxlen: {maxlen}")

    if n_include:
        print(f"n_include: {n_include}")

    if n_exclude:
        print(f"n_exclude: {n_exclude}")


## Path from json_str
def json_str_to_path(json_str):
    allele_path = json_str
    allele_path = re.sub("[{}]", "", allele_path)
    allele_path = re.sub('"', '', allele_path)
    allele_path = re.sub('\s', '', allele_path)
    allele_path = re.sub('\W', '_', allele_path)
    ## print(f"allele path: '{allele_path}'")
    return allele_path


if args.verbose:
    print("{:#^60}".format(f" Path "))
    print(f"path: '{json_str_to_path(json.dumps(d))}'")


while True:
    name = sys.stdin.readline().strip()
    seq = sys.stdin.readline().strip()
    opt = sys.stdin.readline().strip()
    phred = sys.stdin.readline().strip()


    ## normal break of loop
    if not name:
        break


    ## check read name
    if not name.startswith("@"):
        print("Error! Unexpected read name:", name)
        exit(1)
    

    include_read = True


    if args.verbose:
        print("{:#^60}".format(f" Read "))
        ## print(f"seq: '{seq}'")
    

    ## filter for minlen
    if minlen:
        if args.verbose:
            print("{:-^60}".format(f" minlen "))
        
        if minlen>len(seq):
            include_read = False

            if args.verbose:
                print("=> exclude")
        else:
            if args.verbose:
                 print("=> include")


    ## filter for maxlen
    if maxlen:
        if args.verbose:
            print("{:-^60}".format(f" maxlen "))

        if maxlen<len(seq):
            include_read = False

            if args.verbose:
                print("=> exclude")
        else:
            if args.verbose:
                print("=> include")


    ## include
    if n_include > 0:

        if args.verbose:
            print("{:-^60}".format(f" include "))

        include_cond_met = []

        for allele in d["include"].keys():
            m  = re.findall(allele, seq)
            ## print("thres:", d["include"][allele]["n"])
            ## print(allele, ":", m, len(m))
            if len(m) >= d["include"][allele]["n"]:
                include_cond_met.append(True)
            else:
                include_cond_met.append(False)
        ## print(f"include_cond_met: {include_cond_met}")
        ## print(all(include_cond_met))

        if not all(include_cond_met):
            ## print("=> exclude")
            include_read = False
        ## else:
            ## print("=> include")


    ## exclude
    if n_exclude > 0:

        if args.verbose:
            print("{:-^60}".format(f" exclude "))
        exclude_cond_met = []

        for allele in d["exclude"].keys():
            m  = re.findall(allele, seq)
            ## print("thres:", d["exclude"][allele]["n"])
            ## print(allele, ":", m, len(m))
            if len(m) >= d["exclude"][allele]["n"]:
                exclude_cond_met.append(True)
            else:
                exclude_cond_met.append(False)
        ## print(f"exclude_cond_met: {exclude_cond_met}")
        ## print(all(exclude_cond_met))

        if all(exclude_cond_met):
            ## print("=> exclude")
            include_read = False
        ## else:
            ## print("=> include")


    ## print("include_read:", include_read)

    ## output
    if args.verbose:
        print("{:-^60}".format(f" fastq "))
        print("1:", name)
        print("2:", seq)
        print("3:", opt)
        print("4:", phred)
        if include_read:
            print("==> INCLUDE (stdout)")
        else:
            print("==> EXCLUDE (stderr)")
    else:
        if include_read:
            ##print("==> INCLUDE")
            print(name)
            print(seq)
            print(opt)
            print(phred)
        else:
            ##print("==> EXCLUDE")
            print(name, file=sys.stderr)
            print(seq, file=sys.stderr)
            print(opt, file=sys.stderr)
            print(phred, file=sys.stderr)


    ## ## save includes to file
    ## if args.include_fastq:
    ##     print(f"Saving INCLUDES to file: {args.include_fastq}")
    ## if args.exclude_fastq:
    ##     print(f"Saving EXCLUDES to file: {args.exclude_fastq}")

