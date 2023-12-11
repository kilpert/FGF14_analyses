#!/usr/bin/env python

import os
import subprocess
import sys
import re
import json


print("{:#^60}".format(" Args "))

sample = sys.argv[1]
alleles_json_str = sys.argv[2]
input_fastq = sys.argv[3]
output_dir = sys.argv[4]


print("sample:", sample)
print("alleles_json_str:", alleles_json_str)
print("input_fastq:", input_fastq)
print("output_dir", output_dir)


print("{:#^60}".format(" Alleles "))
d = json.loads(alleles_json_str)
print(json.dumps(d, indent=4))


print("{:#^60}".format(f" Filter by alleles "))

remove_list = []

is_first = True

for allele in list(d.keys()) + ["other"]:
    print("{:#^60}".format(f" {allele} "))

    print("is_first:", is_first)

    ## filter by allele
    if is_first:
        infile = input_fastq
    else:
        infile = outfile_other

    outfile = os.path.join(output_dir, f"{sample}.{allele}.fastq.gz")

    if allele != "other":
        outfile_other = os.path.join(output_dir, f"{sample}.{allele}.other.fastq.gz")
    else:
        outfile_other = None

    print("infile:", infile)
    print("outfile:", outfile)
    if outfile_other:
        print("outfile_other:", outfile_other)

    if allele != "other":
        cmd = f"zcat {infile} | python workflow/scripts/filter_fastq_for_alleles.py --json-str '{json.dumps(d[allele])}' > >(gzip >{outfile}) 2> >(gzip >{outfile_other})"

        remove_list.append(outfile_other)
    else:
        cmd = f"cp {infile} {outfile}"
    
    ##cmd = " ".join(cmd)
    print(f"cmd: {cmd}")
    subprocess.run(cmd, check=True, shell=True, executable='bash')

    is_first = False


## remove files
print("{:#^60}".format(f" Remove temp files "))
print(remove_list, len(remove_list))
for f in remove_list:
    os.remove(f)

