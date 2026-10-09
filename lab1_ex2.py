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


def translate_strict_orf(sequence: str):
    """
    Traduce o secventa doar daca incepe de la un codon de start (AUG)
    si se finalizeaza obligatoriu la un codon Stop (UAA, UAG, UGA).
    """

    seq = sequence.strip().upper().replace("T", "U")

    start_pos = seq.find("AUG")
    if start_pos == -1:
        return "Eroare: Nu a fost gasit niciun codon de start (AUG)."

    coding_seq = seq[start_pos:]

    amino_acids = []
    stop_found = False

    for i in range(0, len(coding_seq) - len(coding_seq) % 3, 3):
        codon = coding_seq[i : i + 3]
        amino_acid = GENETIC_CODE.get(codon, "Unknown")

        if amino_acid == "Stop":
            stop_found = True
            break
        
        amino_acids.append(amino_acid)

    if not stop_found:
        return "Eroare: Secventa a inceput de la AUG, dar nu s-a terminat cu un codon STOP."

    return "-".join(amino_acids)


if __name__ == "__main__":
    test1 = "CCGAUGGCUUAUUAA"
    print("Test 1 (Valid):", test1)
    print("Rezultat:", translate_strict_orf(test1))

    test2 = "GCUUAUUAAGGC"
    print("\nTest 2 (Fara AUG):", test2)
    print("Rezultat:", translate_strict_orf(test2))

    test3 = "AUGGCUUAUCCC"
    print("\nTest 3 (Fara Stop):", test3)
    print("Rezultat:", translate_strict_orf(test3))