#!/usr/bin/env python3
import csv
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--macrel", required=True)
parser.add_argument("--anticp", required=True)
parser.add_argument("--fasta", required=True)
args = parser.parse_args()

# Coleta dos IDs
macrel_ids = set()
with open(args.macrel) as f:
    for line in f:
        if line.startswith("Access"):
            header = line.strip().split("\t")
            break
    reader = csv.DictReader(f, fieldnames=header, delimiter="\t")
    for row in reader:
        if row["AMP_family"] != "None":
            macrel_ids.add(row["Access"])

anticp_ids = set()
with open(args.anticp) as f:
    for line in f:
        if line.startswith(">"):
            parts = line[1:].strip().split(",")
            seq_id = parts[0].split()[0]
            anticp_ids.add(seq_id)

# Junta tudo
all_ids = macrel_ids.union(anticp_ids)

# Escreve arquivo
with open("good_ids.txt", "w") as f:
    for seq_id in all_ids:
        f.write(seq_id + "\n")