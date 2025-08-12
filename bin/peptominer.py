#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys
from pathlib import Path

def get_script_dir():
    script_path = Path(os.path.realpath(__file__))
    return script_path.parent

def detect_mode(input_dir):
    fastq_files = [
        f for f in os.listdir(input_dir)
        if f.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))
    ]
    if len(fastq_files) == 1:
        return "single"
    elif len(fastq_files) == 2:
        return "paired"
    elif len(fastq_files) == 3:
        return "mixed"
    else:
        print(f"Erro: Esperado 1 (single-end), 2 (paired-end) arquivos, ou 3 (mixed), mas foram encontrados {len(fastq_files)}.", file=sys.stderr)
        sys.exit(1)

def run_pipeline(input_dir, output_dir=None, k_min=21, k_max=25, threads=12, use_kraken=False, use_map=False):
    script_dir = get_script_dir()
    pipeline_path = script_dir.parent / "main.nf"
    
    if not pipeline_path.exists():
        print(f"Error: Pipeline file not found at {pipeline_path}", file=sys.stderr)
        sys.exit(1)

    mode = detect_mode(input_dir)

    command = [
        "nextflow", "run", str(pipeline_path),
        "--input", input_dir,
        "--k_min", str(k_min),
        "--k_max", str(k_max),
        "--threads", str(threads),
        "--mode", mode,
        "-resume"
    ]

    if output_dir is not None:
        command.extend(["--output", output_dir])
    
    if use_kraken:
        command.extend(["--use_kraken", "true"])

    if use_map:
        command.extend(["--map", "true"])

    subprocess.run(command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PeptoMiner: Pipeline for therapeutic peptide discovery")
    parser.add_argument("-i", "--input", required=True, help="Directory containing FASTQ files")
    parser.add_argument("-o", "--output", help="Output directory")
    parser.add_argument("--k_min", type=int, default=21, help="Minimum k-mer size (default: 21)")
    parser.add_argument("--k_max", type=int, default=25, help="Maximum k-mer size (default: 25)")
    parser.add_argument("-t", "--threads", type=int, default=12, help="Number of threads (default: 12)")
    parser.add_argument("--kraken", action="store_true", help="Enable Kraken2 taxonomic classification")
    parser.add_argument("--map", action="store_true", help="Enable genome mapping step")
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.input):
        print(f"Error: Input directory not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    run_pipeline(args.input, args.output, args.k_min, args.k_max, args.threads, args.kraken, args.map)