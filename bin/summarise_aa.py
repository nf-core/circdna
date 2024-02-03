#!/usr/bin/env python
import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser()

parser.add_argument("--summary_file", "-s", type=str, required=True)
parser.add_argument("--class_file", "-c", type=str, required=True)
parser.add_argument("--id", type=str, required=True)
parser.add_argument("--output", type=str, required=True)
args = parser.parse_args()

file = args.summary_file
lines = open(file, "r").readlines()[2:]

d = {}
i = 0
for line in lines:
    line = line.strip().strip("-")
    res = re.sub(r"\[amplicon[0-9]*\] ", "", line).strip()
    if res == "":
        continue

    l = res.split(" = ")
    header = l[0]
    value = l[1]

    if header in d.keys():
        d[header].append(value)
    else:
        d[header] = [value]

df_summary = pd.DataFrame.from_dict(d)
df_summary.insert(loc=0, column="id", value=args.id)

df_class = pd.read_table(args.class_file, sep="\t")
df_class["amplicon_number"] = df_class["amplicon_number"].str.replace("amplicon", "", regex=False)
df_class = df_class.rename(columns={"sample_name": "id", "amplicon_number": "AmpliconID"})

df_full = pd.merge(df_class, df_summary)
df_full.to_csv(args.output, sep="\t", index=False)
