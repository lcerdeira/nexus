#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parâmetros de entrada (podem ser alterados via linha de comando)
// params.short_reads_1 = "${projectDir}/../tests/illumina_R1.fastq.gz"
// params.short_reads_2 = "${projectDir}/../tests/illumina_R2.fastq.gz"
// params.long_reads    = "${projectDir}/../tests/nanopore.fastq.gz"
// params.reference     = "${projectDir}/../tests/reference.fasta" // Para o QUAST
arams.ont_pod5      = "${projectDir}/data/ont_reads/"          // Pasta com sinal elétrico
params.pacbio_hifi   = "${projectDir}/data/pacbio_hifi.fastq.gz" // Leituras de alta precisão
params.dorado_model  = "sup@v4.3.0"                              // Modelo Super Accuracy da Nanopore
params.genome_size   = "5.4m"                                    // Necessário para Canu (ex: Klebsiella)
params.outdir        = "${projectDir}/results_benchmark"
// params.outdir        = "${projectDir}/../results"

// --- DEFINIÇÃO DOS PROCESSOS (AS "MÁQUINAS" DA ESTEIRA) ---

// 1. BASECALLING COM INTELIGÊNCIA ARTIFICIAL (GPU)
process RUN_DORADO {
    publishDir "${params.outdir}/dorado", mode: 'copy'
    input: path pod5_dir
    output: path "ont_hq_reads.fastq.gz", emit: ont_fastq
    
    script:
    """
    # Exige que o binário do Dorado esteja instalado no HPC
    dorado basecaller ${params.dorado_model} ${pod5_dir} > ont_hq_reads.fastq
    gzip ont_hq_reads.fastq
    """
}

// 2. MONTADORES OXFORD NANOPORE (CPU)
process RUN_FLYE {
    conda "bioconda::flye=2.9.3"
    input: path ont_fastq
    output: path "flye_assembly.fasta", emit: assembly
    script: "flye --nano-hq ${ont_fastq} --out-dir . --threads ${task.cpus} && mv assembly.fasta flye_assembly.fasta"
}

process RUN_SHASTA {
    conda "bioconda::shasta=0.11.1"
    input: path ont_fastq
    output: path "shasta_assembly.fasta", emit: assembly
    script: "shasta --input ${ont_fastq} --config Nanopore-May2022 --threads ${task.cpus} && mv ShastaRun/Assembly.fasta shasta_assembly.fasta"
}

// 3. MONTADORES PACBIO HIFI (CPU)
process RUN_HIFIASM {
    conda "bioconda::hifiasm=0.19.8"
    input: path pacbio_fastq
    output: path "hifiasm_assembly.fasta", emit: assembly
    script: 
    """
    hifiasm -o hifiasm_asm -t ${task.cpus} ${pacbio_fastq}
    # Hifiasm gera grafos (.gfa), precisamos extrair a sequência consenso primária para FASTA:
    awk '/^S/{print ">"\$2"\\n"\$3}' hifiasm_asm.bp.p_ctg.gfa > hifiasm_assembly.fasta
    """
}

// 4. VETERANO OLC (Híbrido)
process RUN_CANU {
    conda "bioconda::canu=2.2"
    input: path ont_fastq
    output: path "canu_assembly.fasta", emit: assembly
    script: "canu -p canu_asm -d . genomeSize=${params.genome_size} -nanopore ${ont_fastq} maxThreads=${task.cpus} && mv canu_asm.contigs.fasta canu_assembly.fasta"
}

// 1. O Nosso Montador (Simulação de chamada do Rust)
process RUN_NEXUS {
    publishDir "${params.outdir}/nexus", mode: 'copy'

    input:
    // path short1
    // path short2
    // path long_reads
    path ont_fastq
    path pacbio_fastq

    output:
    path "nexus_assembly.fasta", emit: assembly

    script:

    """
    echo "Executando o algoritmo Dragon/Nexus com leituras ONT e/ou PacBio..."
    # AQUI ENTRA O COMANDO DO SEU ALGORITMO. Exemplo abstrato:
    # nexus_assemble --long-reads ${ont_fastq} --hifi ${pacbio_fastq} --threads ${task.cpus} -o nexus_assembly.fasta
    touch nexus_assembly.fasta 
    """
}

// 2. O Concorrente Híbrido Clássico
// process RUN_UNICYCLER {
//     publishDir "${params.outdir}/unicycler", mode: 'copy'

//     input:
//     path short1
//     path short2
//     path long_reads

//     output:
//     path "assembly.fasta", emit: assembly

//     script:
//     """
//     unicycler -1 $short1 -2 $short2 -l $long_reads -o . --keep 0
//     """
// }

// 6. BENCHMARKING FINAL (QUAST)
process RUN_QUAST {
    publishDir "${params.outdir}/quast_benchmark", mode: 'copy'
    conda "bioconda::quast=5.2.0"
    input: path assemblies
    output: path "quast_results/*"
    script: "quast.py ${assemblies.join(' ')} -o quast_results --threads ${task.cpus}"
}

// ORQUESTRAÇÃO
workflow {
    ch_pod5    = Channel.fromPath(params.ont_pod5)
    ch_pacbio  = file(params.pacbio_hifi)

    // A GPU entra em ação
    dorado_reads = RUN_DORADO(ch_pod5)

    // Paralelismo massivo nos 20 nós de CPU do HPC
    flye_asm    = RUN_FLYE(dorado_reads)
    shasta_asm  = RUN_SHASTA(dorado_reads)
    canu_asm    = RUN_CANU(dorado_reads)
    hifiasm_asm = RUN_HIFIASM(ch_pacbio)
    
    // Seu algoritmo recebendo ambas as leituras
    nexus_asm   = RUN_NEXUS(dorado_reads, ch_pacbio)

    // Agrupa todos os genomas gerados e manda para o juiz (QUAST)
    ch_all_assemblies = flye_asm.mix(shasta_asm, canu_asm, hifiasm_asm, nexus_asm).collect()
    RUN_QUAST(ch_all_assemblies)
}