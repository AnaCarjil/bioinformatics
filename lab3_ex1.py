import sys
from collections import Counter

# Standard Genetic Code Table: DNA Codon -> Amino Acid
GENETIC_CODE = {
    # Phenylalanine & Leucine
    'TTT': 'Phenylalanine (Phe)', 'TTC': 'Phenylalanine (Phe)',
    'TTA': 'Leucine (Leu)', 'TTG': 'Leucine (Leu)',
    'CTT': 'Leucine (Leu)', 'CTC': 'Leucine (Leu)',
    'CTA': 'Leucine (Leu)', 'CTG': 'Leucine (Leu)',
    # Isoleucine & Methionine (Start)
    'ATT': 'Isoleucine (Ile)', 'ATC': 'Isoleucine (Ile)',
    'ATA': 'Isoleucine (Ile)', 'ATG': 'Methionine (Met / Start)',
    # Valine
    'GTT': 'Valine (Val)', 'GTC': 'Valine (Val)',
    'GTA': 'Valine (Val)', 'GTG': 'Valine (Val)',
    # Serine
    'TCT': 'Serine (Ser)', 'TCC': 'Serine (Ser)',
    'TCA': 'Serine (Ser)', 'TCG': 'Serine (Ser)',
    'AGT': 'Serine (Ser)', 'AGC': 'Serine (Ser)',
    # Proline
    'CCT': 'Proline (Pro)', 'CCC': 'Proline (Pro)',
    'CCA': 'Proline (Pro)', 'CCG': 'Proline (Pro)',
    # Threonine
    'ACT': 'Threonine (Thr)', 'ACC': 'Threonine (Thr)',
    'ACA': 'Threonine (Thr)', 'ACG': 'Threonine (Thr)',
    # Alanine
    'GCT': 'Alanine (Ala)', 'GCC': 'Alanine (Ala)',
    'GCA': 'Alanine (Ala)', 'GCG': 'Alanine (Ala)',
    # Tyrosine & Stop codons
    'TAT': 'Tyrosine (Tyr)', 'TAC': 'Tyrosine (Tyr)',
    'TAA': 'Stop', 'TAG': 'Stop', 'TGA': 'Stop',
    # Histidine & Glutamine
    'CAT': 'Histidine (His)', 'CAC': 'Histidine (His)',
    'CAA': 'Glutamine (Gln)', 'CAG': 'Glutamine (Gln)',
    # Asparagine & Lysine
    'AAT': 'Asparagine (Asn)', 'AAC': 'Asparagine (Asn)',
    'AAA': 'Lysine (Lys)', 'AAG': 'Lysine (Lys)',
    # Aspartic Acid & Glutamic Acid
    'GAT': 'Aspartic Acid (Asp)', 'GAC': 'Aspartic Acid (Asp)',
    'GAA': 'Glutamic Acid (Glu)', 'GAG': 'Glutamic Acid (Glu)',
    # Cysteine & Tryptophan
    'TGT': 'Cysteine (Cys)', 'TGC': 'Cysteine (Cys)',
    'TGG': 'Tryptophan (Trp)',
    # Arginine
    'CGT': 'Arginine (Arg)', 'CGC': 'Arginine (Arg)',
    'CGA': 'Arginine (Arg)', 'CGG': 'Arginine (Arg)',
    'AGA': 'Arginine (Arg)', 'AGG': 'Arginine (Arg)',
    # Glycine
    'GGT': 'Glycine (Gly)', 'GGC': 'Glycine (Gly)',
    'GGA': 'Glycine (Gly)', 'GGG': 'Glycine (Gly)',
}

def parse_fasta(filepath):
    """Reads a FASTA file and returns the concatenated sequence string."""
    sequence_parts = []
    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            sequence_parts.append(line.upper())
    return "".join(sequence_parts)

def count_codons(sequence, window_step=3):
    """
    Counts codon frequencies using a non-overlapping reading frame (step=3).
    Ignores incomplete trailing triplets and ambiguous IUPAC bases (N, R, Y, etc.).
    """
    codon_counts = Counter()
    valid_bases = {'A', 'C', 'G', 'T'}
    
    for i in range(0, len(sequence) - 2, window_step):
        codon = sequence[i:i + 3]
        if all(base in valid_bases for base in codon):
            codon_counts[codon] += 1
            
    return codon_counts

def display_top_codons(codon_counts, top_n=10):
    """Displays formatted results with DNA codon, RNA equivalent, count, frequency %, and amino acid."""
    total_codons = sum(codon_counts.values())
    top_codons = codon_counts.most_common(top_n)
    
    print("\n" + "=" * 76)
    print(f"{'Rank':<5}{'DNA':<6}{'RNA':<6}{'Count':<10}{'Frequency (%)':<15}{'Amino Acid':<30}")
    print("=" * 76)
    
    for rank, (dna_codon, count) in enumerate(top_codons, 1):
        rna_codon = dna_codon.replace('T', 'U')
        freq = (count / total_codons) * 100 if total_codons > 0 else 0
        amino_acid = GENETIC_CODE.get(dna_codon, "Unknown")
        print(f"{rank:<5}{dna_codon:<6}{rna_codon:<6}{count:<10}{freq:<15.2f}{amino_acid:<30}")
    print("=" * 76 + "\n")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        path = input("Enter path to FASTA file: ").strip()
        
    try:
        seq = parse_fasta(path)
        print(f"Loaded sequence length: {len(seq)} bp")
        counts = count_codons(seq, window_step=3)
        display_top_codons(counts, top_n=10)
    except FileNotFoundError:
        print(f"Error: File not found at '{path}'.")