process fastp_single {
    input:
    path fastq_file
    
    output:
    path "output_single.fastq"
    
    script:
    """
    fastp -i ${fastq_file} -o output_single.fastq \\
        --detect_adapter_for_pe \\
        --cut_front \\
        --cut_tail \\
        --length_required 50 \\
        --trim_front1 10 \\
        --thread 4
    """
}

process fastp_paired {
    input:
    tuple path(fastq_r1), path(fastq_r2)
    
    output:
    tuple path("output_paired_R1.fastq"), path("output_paired_R2.fastq")
    
    script:
    """
    fastp -i ${fastq_r1} -I ${fastq_r2} -o output_paired_R1.fastq -O output_paired_R2.fastq \\
        --detect_adapter_for_pe \\
        --cut_front \\
        --cut_tail \\
        --length_required 50 \\
        --trim_front1 10 \\
        --trim_front2 10 \\
        --thread 4
    """
}
process kraken2 {
    publishDir "${params.output}", mode: 'copy'

    input:
    path clean_fastq_files

    output:
    path "kraken2_output.txt"

    script:
    if (params.mode == "single")
        """
        kraken2 --db ${params.kraken2_db} \
                --threads 4 \
                --report kraken2_output.txt \
                --output kraken2_classification.txt \
                ${clean_fastq_files}
        """
    else if (params.mode == "paired")
        """
        kraken2 --db ${params.kraken2_db} \
                --threads 4 \
                --report kraken2_output.txt \
                --output kraken2_classification.txt \
                --paired ${clean_fastq_files[0]} ${clean_fastq_files[1]}
        """
    else if (params.mode == "mixed")
    """
    kraken2 --db ${params.kraken2_db} \
            --threads 4 \
            --report kraken2_output.txt \
            --output kraken2_classification.txt \
            --paired ${clean_fastq_files[0]} ${clean_fastq_files[1]} \
            ${clean_fastq_files[2]}
    """
}

process kraken2_mixed {
    publishDir "${params.output}", mode: 'copy'

    input:
    path clean_fastq_files

    output:
    path "kraken2_output.txt"

    script:
    """
    kraken2 --db ${params.kraken2_db} \\
            --threads 4 \\
            --report kraken2_paired_report.txt \\
            --output kraken2_paired_classification.txt \\
            --paired ${clean_fastq_files[0]} ${clean_fastq_files[1]}

    kraken2 --db ${params.kraken2_db} \\
            --threads 4 \\
            --report kraken2_orphans_report.txt \\
            --output kraken2_orphans_classification.txt \\
            ${clean_fastq_files[2]}

    cat kraken2_paired_report.txt kraken2_orphans_report.txt > kraken2_output.txt
    """
}

process kraken2_contigs {
    publishDir "${params.output}", mode: 'copy'

    input:
    path megahit_output

    output:
    path "kraken2_contigs_output.txt"

    script:
        """
        kraken2 --db ${params.kraken2_db} \
                --threads 4 \
                --report kraken2_contigs_output.txt \
                --output kraken2_contigs_classification.txt \
                ${megahit_output}
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
        megahit -r ${clean_fastq_files} -o megahit_output \\
                -t ${params.threads} --mem-flag 0 \\
                --k-min 21 --k-max 51 --k-step 10
        cp megahit_output/final.contigs.fa .
        """
    else if (params.mode == "paired")
        """
        megahit -1 ${clean_fastq_files[0]} -2 ${clean_fastq_files[1]} -o megahit_output \\
                -t ${params.threads} --mem-flag 0 \\
                --k-min 21 --k-max 51 --k-step 10
        cp megahit_output/final.contigs.fa .
        """
    else if (params.mode == "mixed")
        """
        megahit -1 ${clean_fastq_files[0]} -2 ${clean_fastq_files[1]} -r ${clean_fastq_files[2]} -o megahit_output \\
                -t ${params.threads} --mem-flag 0 \\
                --k-min 21 --k-max 51 --k-step 10
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
             -p meta \
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
    input_dir = file("${params.input}")
    
    // Assegurar que o diretório de entrada existe
    if (!input_dir.exists()) {
        error "Diretório de entrada não encontrado: ${params.input}"
    }

    if (params.mode == "single") {
        input_channel = Channel.fromPath("${input_dir}/*.fastq")
        fastp_single(input_channel)
        kraken2(fastp_single.out)
        megahit(fastp_single.out)
        kraken2_contigs(megahit.out)
    } 
    else if (params.mode == "paired") {
        input_channel = Channel.fromFilePairs("${input_dir}/*_{1,2}.fastq")

        fastp_paired(input_channel.map { it[1] })

        kraken2(fastp_paired.out.collect())
        megahit(fastp_paired.out.collect())
        kraken2_contigs(megahit.out)
    }
    else if (params.mode == "mixed") {
    paired_channel = Channel.fromFilePairs("${input_dir}/*_{1,2}.fastq")
    orphans_channel = Channel.fromPath("${input_dir}/*.fastq")
        .filter { !it.name.matches(".*_[12]\\.fastq") }

    fastp_paired(paired_channel.map { it[1] })
    fastp_single(orphans_channel)

    qc_outputs = fastp_paired.out.merge(fastp_single.out)

    kraken2_mixed(qc_outputs.collect())
    megahit(qc_outputs)

    kraken2_contigs(megahit.out)
    }

    pyrodigal(megahit.out)
    anticp(pyrodigal.out)
    macrel(pyrodigal.out)
    if (params.map) {
        eggnog(pyrodigal.out)
    }
    if (params.mode != "mixed") {
        krona(kraken2.out)
    }
    if (params.mode == "mixed") {
        krona(kraken2_mixed.out)
    }
}