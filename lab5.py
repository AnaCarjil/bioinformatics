import random

# DNA complement lookup mapping
COMPLEMENT = str.maketrans("ACGT", "TGCA")

def reverse_complement(seq: str) -> str:
    """Returns the reverse complement of a DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]

def generate_random_dna(length: int) -> str:
    """Generates a random DNA sequence of specified length."""
    return "".join(random.choices(["A", "C", "G", "T"], k=length))

def create_transposon(internal_len: int = 30, ir_len: int = 8):
    """Creates a transposon with terminal inverted repeats (IRs)."""
    ir = generate_random_dna(ir_len)
    internal_seq = generate_random_dna(internal_len)
    transposon_seq = ir + internal_seq + reverse_complement(ir)
    return transposon_seq, ir

def insert_transposons(
    host_length: int = 300, 
    num_transposons: int = 3, 
    tsd_len: int = 5, 
    ir_len: int = 8, 
    internal_len: int = 20
):
    """Simulates insertion of transposons into host DNA creating Target Site Duplications (TSDs)."""
    host_dna = list(generate_random_dna(host_length))
    inserted_records = []

    for i in range(num_transposons):
        transposon, ir_seq = create_transposon(internal_len=internal_len, ir_len=ir_len)
        insert_pos = random.randint(10, max(11, len(host_dna) - (tsd_len + 10)))
        tsd_seq = "".join(host_dna[insert_pos : insert_pos + tsd_len])
        
        full_insertion = tsd_seq + transposon + tsd_seq
        host_dna[insert_pos : insert_pos + tsd_len] = list(full_insertion)
        
        inserted_records.append({
            "id": i + 1,
            "target_site_tsd": tsd_seq,
            "ir_seq": ir_seq,
        })

    return "".join(host_dna), inserted_records

def detect_transposons(dna_seq: str, min_ir_len: int = 6, max_ir_len: int = 15, tsd_len: int = 5):
    """Detects transposable elements without duplicate hits for the same insertion site."""
    detected = []
    n = len(dna_seq)

    for i in range(n - (2 * tsd_len + 2 * min_ir_len)):
        tsd_candidate = dna_seq[i : i + tsd_len]

        for j in range(i + tsd_len + 2 * min_ir_len, n - tsd_len + 1):
            if dna_seq[j : j + tsd_len] == tsd_candidate:
                transposon_candidate = dna_seq[i + tsd_len : j]
                
                # Check IR sizes from longest to shortest to get the maximum matching IR length
                longest_ir = None
                max_size = min(max_ir_len, len(transposon_candidate) // 2)
                
                for ir_size in range(max_size, min_ir_len - 1, -1):
                    left_ir = transposon_candidate[:ir_size]
                    right_ir = transposon_candidate[-ir_size:]
                    
                    if right_ir == reverse_complement(left_ir):
                        longest_ir = left_ir
                        break  # Stop at the longest valid IR
                
                if longest_ir:
                    detected.append({
                        "start": i,
                        "end": j + tsd_len,
                        "transposon_start": i + tsd_len,
                        "transposon_end": j,
                        "tsd": tsd_candidate,
                        "ir": longest_ir,
                    })

    return detected

if __name__ == "__main__":
    print("--- Simulating DNA Transposition ---")
    dna_sequence, _ = insert_transposons(host_length=300, num_transposons=3)
    print(f"Generated Sequence Length: {len(dna_sequence)} bp\n")
    
    print("--- Running Transposon Detection Algorithm ---")
    results = detect_transposons(dna_sequence)
    
    if not results:
        print("No elements detected.")
    for idx, match in enumerate(results, 1):
        print(f"Element #{idx}:")
        print(f"  Full Bound (with TSDs): [{match['start']}, {match['end']}]")
        print(f"  Transposon Bounds: [{match['transposon_start']}, {match['transposon_end']}]")
        print(f"  Target Site Duplication (TSD): {match['tsd']}")
        print(f"  Inverted Repeat (IR): {match['ir']}\n")