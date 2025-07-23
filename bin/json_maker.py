import csv
import json
from Bio import SeqIO

# 1. Carregar taxonomia dos contigs (Kraken2)
contig_tax = {}
with open("kraken2_contigs_output.txt") as f:
    for line in f:
        parts = line.strip().split("\t")
        if parts[0] == "C":  # só pega contigs classificados
            contig_id = "_".join(parts[1].split("_")[:2])
            taxid = parts[2]
            contig_tax[contig_id] = taxid

# 2. Carregar eggnog
eggnog = {}
with open("eggnog_output.emapper.annotations") as f:
    # Pular linhas até o cabeçalho real
    for line in f:
        if line.startswith("#query"):
            header = line.strip().split("\t")
            break

    reader = csv.DictReader(f, fieldnames=header, delimiter="\t")
    for row in reader:
        query = row["#query"]
        eggnog[query] = row

# 3. Carregar proteínas para mapear contig
prot2contig = {}
for record in SeqIO.parse("proteins.faa", "fasta"):
    prot_id = record.id
    header = record.description
    contig_id = header.split()[0].split("_")[0]  # ajustar conforme seu header
    prot2contig[prot_id] = {
        "contig": contig_id,
        "sequence": str(record.seq)
    }

# 4. Carregar peptídeos antimicrobianos (MACREL)
amp_peptides = {}
with open("macrel.out.prediction.tsv") as f:
    # pular comentários iniciais
    for line in f:
        if line.startswith("Access"):
            header = line.strip().split("\t")
            break

    reader = csv.DictReader(f, fieldnames=header, delimiter="\t")
    for row in reader:
        amp_peptides[row["Access"]] = row

# 5. Carregar peptídeos anticancer (AntiCP)
anticancer_peptides = {}
with open("anticp_output.csv") as f:
    for line in f:
        if line.startswith(">"):
            # Remove o ">" e separa por vírgula
            parts = line[1:].strip().split(",")
            # Define os campos esperados
            row = {
                "Sequence_ID": parts[0].split()[0],  # >k51_2152_2 etc.
                "Sequence": parts[1],
                "Prediction": parts[2],
                "Score": float(parts[3])
            }
            anticancer_peptides[row["Sequence_ID"]] = row

# 6. Montar peptídeos finais
peptideos_final = []

# Unir os dados
all_ids = set(amp_peptides.keys()) | set(anticancer_peptides.keys())
for pid in all_ids:
    amp = amp_peptides.get(pid)
    acp = anticancer_peptides.get(pid)

    tipo = []
    if amp: tipo.append("antimicrobial")
    if acp: tipo.append("anticancer")

    # Proteína e contig de origem
    prot_info = prot2contig.get(pid, {})
    contig_id = prot_info.get("contig", "unknown")

    # Taxonomia
    taxonomy = contig_tax.get(contig_id, "unknown")

    # Anotação funcional
    annotation = eggnog.get(pid, {})

    peptideo_obj = {
        "id": pid,
        "sequence": amp["Sequence"] if amp else acp["Sequence"],
        "type": tipo,
        "scores": {},
        "origin": {
            "protein_id": pid,
            "contig": contig_id,
            "protein_file": "proteins.faa",
            "contig_file": "final.contigs.fa"
        },
        "taxonomy": {
            "lineage": taxonomy,
            "source": {
                "program": "Kraken2",
                "file": "kraken2_output.tsv"
            }
        },
        "annotations": {
            "COG": {
                "value": annotation.get("COG_category", ""),
                "source": {
                    "program": "eggNOG-mapper",
                    "file": "eggnog_output.tsv"
                }
            },
            "KEGG": annotation.get("KEGG_Pathway", "").split(","),
            "GO": annotation.get("GOs", "").split(","),
            "EC": annotation.get("EC", "")
        }
    }

    if amp:
        peptideo_obj["scores"]["amp"] = {
            "value": float(amp["AMP_probability"]),
            "source": {
                "program": "MACREL",
                "file": "macrel.out.prediction.tsv"
            }
        }
        peptideo_obj["scores"]["hemolytic"] = {
            "value": float(amp["Hemolytic_probability"]),
            "label": amp["Hemolytic"],
            "source": {
                "program": "MACREL",
                "file": "macrel.out.prediction.tsv"
            }
        }

    if acp:
        peptideo_obj["scores"]["anticancer"] = {
            "value": float(acp["Score"]),
            "source": {
                "program": "AntiCP",
                "file": "anticp_output.csv"
            }
        }

    peptideos_final.append(peptideo_obj)

# 7. Salvar em JSON
with open("peptideos_integrados.json", "w") as out:
    json.dump(peptideos_final, out, indent=2)
