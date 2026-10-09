# Implement an application that converts the coding region of a gene into an amino acid sequence. Use the genetic code table from the slide
GENETIC_CODE = {
    # U-row
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
    "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp",
    
    # C-row
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    
    # A-row
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    
    # G-row
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}


def translate_coding_sequence(sequence: str, stop_at_termination: bool = True) -> list[str]:
    """
    Translates an mRNA or DNA coding sequence into an amino acid sequence.
    
    :param sequence: Nucleotide sequence string (e.g. 'AUGGCU...')
    :param stop_at_termination: If True, halts translation when a Stop codon is encountered.
    :return: List of 3-letter amino acid codes.
    """

    seq = sequence.strip().upper().replace("T", "U")
    
    amino_acids = []
    
    for i in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[i : i + 3]
        amino_acid = GENETIC_CODE.get(codon, "Unknown")
        
        if amino_acid == "Stop":
            if stop_at_termination:
                break
            amino_acids.append("Stop")
        else:
            amino_acids.append(amino_acid)
            
    return amino_acids


if __name__ == "__main__":
    sample_mrna = "AUGGCUUAUUAA"
    result = translate_coding_sequence(sample_mrna)
    print("mRNA:", sample_mrna)
    print("Amino acid chain:", "-".join(result))

    
    sample_dna = "ATGGCTCCTTGAGGA"
    result_dna = translate_coding_sequence(sample_dna)
    print("\nDNA coding sequence:", sample_dna)
    print("Amino acid chain:", "-".join(result_dna))