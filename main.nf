process fastp {

    input:
    path fastq_files

    output:
    path "output.fastq"

    script:
    if (params.mode == "single")
        """
        fastp -i ${fastq_files} -o output.fastq \
            --detect_adapter_for_pe \
            --cut_front \
            --cut_tail \
            --length_required 50 \
            --trim_front1 10 \
            --thread 4
        """
    else
        """
        fastp -i ${fastq_files[0]} -I ${fastq_files[1]} -o output_R1.fastq -O output_R2.fastq \
            --detect_adapter_for_pe \
            --cut_front \
            --cut_tail \
            --length_required 50 \
            --trim_front1 10 \
            --trim_front2 10 \
            --thread 4
        """
}

process kraken2 {
    publishDir "${params.output}", mode: 'copy'

    input:
    path clean_fastq_files

    output:
    path kraken2_report: "kraken2_output.txt" // Dando um nome específico para esta saída
    path kraken2_class: "kraken2_classification.txt" // Dando um nome específico para esta saída

    script:
    if (params.mode == "single")
        """
        kraken2 --db ${params.kraken_db} \
                --threads 4 \
                --report kraken2_output.txt \
                --output kraken2_classification.txt \
                ${clean_fastq_files}
        """
    else
        """
        kraken2 --db ${params.kraken_db} \
                --threads 4 \
                --report kraken2_output.txt \
                --output kraken2_classification.txt \
                --paired ${clean_fastq_files[0]} ${clean_fastq_files[1]}
        """
}

process kraken2_contigs {
    publishDir "${params.output}", mode: 'copy'

    input:
    path megahit_output

    output:
    path "kraken2_contigs_output.txt"

    script:
    if (params.mode == "single")
        """
        kraken2 --db ${params.kraken_db} \
                --threads 4 \
                --report kraken2_contigs_output.txt \
                --output kraken2_contigs_classification.txt \
                ${megahit_output}
        """
    else
        """
        kraken2 --db ${params.kraken_db} \
                --threads 4 \
                --report kraken2_contigs_output.txt \
                --output kraken2_contigs_classification.txt \
                --paired ${megahit_output[0]} ${megahit_output[1]}
        """
}

process megahit {
    publishDir "${params.output}", mode: 'copy'

    input:
    path clean_fastq_files

    output:
    path "final.contigs.fa"

    script:
    if (params.mode == "single")
        """
        megahit -r ${clean_fastq_files} -o megahit_output \
                -t 8 \
                --mem-flag 0 \
                --use-gpu \
                --k-min 21 \
                --k-max 51 \
                --k-step 10
        cp megahit_output/final.contigs.fa .
        """
    else
        """
        megahit -1 ${clean_fastq_files[0]} -2 ${clean_fastq_files[1]} -o megahit_output \
                -t 8 \
                --mem-flag 0 \
                --use-gpu \
                --k-min 21 \
                --k-max 51 \
                --k-step 10
        cp megahit_output/final.contigs.fa .
        """
}

process pyrodigal {
    publishDir "${params.output}", mode: 'copy'

    input:
    path megahit_output

    output:
    path "proteins.faa"

    script:
    """
    prodigal -i ${megahit_output} -o pyrodigal_output.gff \
             -a proteins.faa \
             -d nucleotides.fna \
             -f gff
    """
}

process anticp {
    conda 'envs/peptominer_anticp.yml'
    publishDir "${params.output}", mode: 'copy'

    input:
    path proteins_faa

    output:
    path "anticp_output.csv"

    script:
    """
    anticp2 -i ${proteins_faa} -o anticp_output.csv
    """
}

process macrel {
    publishDir "${params.output}", mode: 'copy'

    input:
    path proteins_faa

    output:
    path "macrel.out.prediction.tsv"

    script:
    """
    macrel peptides -f ${proteins_faa} -o macrel_output
    gunzip -c macrel_output/macrel.out.prediction.gz > macrel.out.prediction.tsv
    """
}

process eggnog {
    publishDir "${params.output}", mode: 'copy'

    input:
    path proteins_faa

    output:
    path "eggnog_output.emapper.annotations"

    script:
    """
    emapper.py -i ${proteins_faa} -o eggnog_output \
                  --cpu 8 \
                  --sensmode fast
    """
}

process krona {
    publishDir "${params.output}", mode: 'copy'

    input:
    path kraken2_report_file 

    output:
    path "krona_output.html"

    script:
    """
    ktImportTaxonomy -o krona_output.html ${kraken2_report_file}
    """
}

workflow {
    def fastq_files = Channel.fromPath('data/*.fastq')
    fastp(fastq_files)
    kraken2(fastp.out.collect())
    megahit(fastp.out.collect())
    kraken2_contigs(megahit.out.collect())
    pyrodigal(megahit.out.collect())
    anticp(pyrodigal.out.collect())
    macrel(pyrodigal.out.collect())
    if (params.map)
        eggnog(pyrodigal.out.collect())
    krona(kraken2.out.collect())
}