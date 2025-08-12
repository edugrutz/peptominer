# PeptoMiner
A Computational Tool for Mining Therapeutic Peptides from Metagenomic Data

⚠️ **Project Status: In Development**

This project is under active development and not yet ready for production use. Features and structure may change at any time.


## 📘 Project Overview

PeptoMiner is a computational tool designed to mine and identify potential therapeutic peptides from metagenomic datasets. It streamlines the discovery process by integrating preprocessing, translation, peptide prediction, and annotation into a single automated pipeline, implemented using **Nextflow**.


## 📦 Input
Raw metagenomic sequences in FASTQ or FASTA format

## 🧬 Example Usage
```
peptomine --input metagenome.fasta --output peptides.fasta
```

Optional parameters:
```
--threads 8
--k_min 8
--k_max 50
```

## 📁 Project Structure
```
PeptoMine/
├── bin
│   ├── peptomine.py
├── envs/
│   ├── environment.yml
│   ├── environment_anticp.yml
├── README.md
├── nextflow.config
└── main_pipeline.nf
```
## 👨‍💻 Author
Eduardo Grutzmann Furtado — Master's student in Biotechnology, specialized in bioinformatics, software development and front-end.

grutzmann9@gmail.com

www.linkedin.com/in/edugrutz

## 📃 License
MIT License

