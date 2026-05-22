use std::collections::{HashMap, HashSet};

// --- [CÓDIGO ANTERIOR DO DBG MANTIDO] ---
#[derive(Debug)]
pub struct Node {
    pub count: usize,
    pub in_edges: Vec<u8>,
    pub out_edges: Vec<u8>,
}

impl Node {
    fn new() -> Self {
        Node {
            count: 0,
            in_edges: Vec::new(),
            out_edges: Vec::new(),
        }
    }
}

pub struct DeBruijnGraph {
    pub k: usize,
    pub nodes: HashMap<Vec<u8>, Node>,
}

impl DeBruijnGraph {
    pub fn new(k_size: usize) -> Self {
        DeBruijnGraph { k: k_size, nodes: HashMap::new() }
    }

    pub fn add_sequence(&mut self, sequence: &[u8]) {
        if sequence.len() < self.k { return; }
        for kmer in sequence.windows(self.k) {
            let node = self.nodes.entry(kmer.to_vec()).or_insert(Node::new());
            node.count += 1;
        }
        for window in sequence.windows(self.k + 1) {
            let left_kmer = &window[0..self.k];
            let right_kmer = &window[1..=self.k];
            let in_base = window[0];
            let out_base = window[self.k];

            if let Some(node) = self.nodes.get_mut(left_kmer) {
                if !node.out_edges.contains(&out_base) { node.out_edges.push(out_base); }
            }
            if let Some(node) = self.nodes.get_mut(right_kmer) {
                if !node.in_edges.contains(&in_base) { node.in_edges.push(in_base); }
            }
        }
    }

    pub fn extract_unitigs(&self) -> Vec<String> {
        let mut unitigs = Vec::new();
        let mut visited = HashSet::new();

        for (kmer, node) in &self.nodes {
            if node.in_edges.len() != 1 {
                for &out_base in &node.out_edges {
                    let mut current_kmer = kmer.clone();
                    let mut unitig = String::from_utf8(current_kmer.clone()).unwrap();
                    let mut current_out_base = out_base;

                    loop {
                        let mut next_kmer = current_kmer[1..].to_vec();
                        next_kmer.push(current_out_base);

                        if let Some(next_node) = self.nodes.get(&next_kmer) {
                            if !visited.insert(next_kmer.clone()) { break; }
                            unitig.push(current_out_base as char);

                            if next_node.in_edges.len() != 1 || next_node.out_edges.len() != 1 { break; }
                            current_kmer = next_kmer;
                            current_out_base = next_node.out_edges[0];
                        } else { break; }
                    }
                    if unitig.len() > self.k { unitigs.push(unitig); }
                }
            }
        }
        unitigs
    }
}
// --- [FIM DO CÓDIGO DO DBG] ---

// --- [MANTENHA A FASE 1 - DBG INTOCADA ACIMA DESTA LINHA] ---

// --- FASE 2 e 3: MOTOR OLC OTIMIZADO E CONSENSO ---

#[derive(Debug)]
pub struct Overlap {
    pub seq_a_index: usize, // Usamos índices numéricos para poupar memória
    pub seq_b_index: usize,
    pub overlap_length: usize,
}

// 1. O Indexador de Prefixos O(N)
pub fn build_prefix_index(sequences: &Vec<String>, k_index: usize) -> HashMap<String, Vec<usize>> {
    let mut index: HashMap<String, Vec<usize>> = HashMap::new();
    
    for (id, seq) in sequences.iter().enumerate() {
        if seq.len() >= k_index {
            // Extrai as primeiras 'k_index' letras da sequência
            let prefix = seq[0..k_index].to_string();
            // Salva no dicionário qual sequência (ID) começa com esse prefixo
            index.entry(prefix).or_insert(Vec::new()).push(id);
        }
    }
    index
}

// 2. Busca de Overlaps Usando o Índice
pub fn find_fast_overlaps(sequences: &Vec<String>, min_overlap: usize) -> Vec<Overlap> {
    let mut overlaps = Vec::new();
    
    // Constrói o índice usando um tamanho de "semente" igual ao overlap mínimo
    let prefix_index = build_prefix_index(sequences, min_overlap);

    for (id_a, seq_a) in sequences.iter().enumerate() {
        if seq_a.len() < min_overlap { continue; }

        // Em vez de testar contra todas as reads, testamos de trás para frente no final da seq_a
        for len in (min_overlap..=seq_a.len()).rev() {
            let suffix_a = &seq_a[seq_a.len() - len..];
            
            // Pegamos apenas a semente (as primeiras letras) do sufixo atual
            let seed = &suffix_a[0..min_overlap];

            // Consultamos o dicionário em tempo O(1)
            if let Some(candidate_ids) = prefix_index.get(seed) {
                for &id_b in candidate_ids {
                    // Ignora auto-alinhamento
                    if id_a == id_b { continue; }

                    let seq_b = &sequences[id_b];
                    
                    // Verifica se o sufixo inteiro realmente bate com o prefixo inteiro
                    if seq_b.len() >= len && suffix_a == &seq_b[0..len] {
                        overlaps.push(Overlap {
                            seq_a_index: id_a,
                            seq_b_index: id_b,
                            overlap_length: len,
                        });
                        break; // Achou o maior overlap, vai para a próxima read
                    }
                }
            }
        }
    }
    overlaps
}

// 3. Fase 4: O Consenso (Gerando a Montagem Final)
pub fn generate_contig(sequences: &Vec<String>, overlaps: &Vec<Overlap>) -> String {
    if sequences.is_empty() { return String::new(); }
    if overlaps.is_empty() { return sequences[0].clone(); }

    // Começamos com a primeira sequência do overlap
    let mut current_contig = sequences[overlaps[0].seq_a_index].clone();
    
    for overlap in overlaps {
        let next_seq = &sequences[overlap.seq_b_index];
        // Adiciona apenas a parte da Sequência B que NÃO estava sobreposta
        let unique_part = &next_seq[overlap.overlap_length..];
        current_contig.push_str(unique_part);
    }
    
    current_contig
}

fn main() {
    println!("=== Nexus Assembler: Executando Motor Híbrido Otimizado ===");

    // Nossos cacos perfeitos vindos do DBG e uma long-read (agora representados em strings puras)
    let hybrid_pool = vec![
        String::from("ATGCGT"),       // Unitig 1
        String::from("CGTACGTAGCTA"), // Long-Read (A Ponte)
        String::from("GCTAGCTAG"),    // Unitig 2
    ];

    println!("\n[Indexando] Construindo dicionário de prefixos em tempo O(N)...");
    
    // Busca os overlaps usando o índice super rápido
    let overlaps = find_fast_overlaps(&hybrid_pool, 3);

    for overlap in &overlaps {
        println!(
            "-> Overlap entre Read {} e Read {} ({} bases)",
            overlap.seq_a_index, overlap.seq_b_index, overlap.overlap_length
        );
    }

    // Gera o contig colando os cacos sem duplicar as letras
    let final_genome = generate_contig(&hybrid_pool, &overlaps);
    
    println!("\n=== RESULTADO FINAL: MONTADOR NEXUS ===");
    println!("Contig montado: {}", final_genome);
    println!("Tamanho final:  {} pb", final_genome.len());
}