import re
import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. INPUT: Restriction Enzymes & DNA Sequence
# ==========================================

# Format: "Enzyme Name": ("5'->3' Recognition Sequence", cleavage_offset_from_5_prime)
RESTRICTION_ENZYMES = {
    "EcoRI":   ("GAATTC", 1),  # G|AATTC
    "BamHI":   ("GGATCC", 1),  # G|GATCC
    "HindIII": ("AAGCTT", 1),  # A|AGCTT
    "TaqI":    ("TCGA",   1),  # T|CGA
    "HaeIII":  ("GGCC",   2),  # GG|CC (blunt)
}

def load_fasta(filepath):
    """Optional helper to load a sequence directly from an NCBI .fasta file."""
    with open(filepath, "r") as f:
        lines = [line.strip() for line in f if not line.startswith(">")]
    return "".join(lines).upper()

# Example bacterial DNA sequence (Replace this with your NCBI sequence or use load_fasta("sequence.fasta"))
SAMPLE_NCBI_DNA = (
    "ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATCTACACCAACGTGACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACGCGAATTATTTTTGATGGCGTTAACTCGGCGTTTCATCTGTGGTGCAACGGGCGCTGGGTCGGTTACGGCCAGGACAGTCGTTTGCCGTCTGAATTTGACCTGAGCGCATTTTTACGCGCCGGAGAAAACCGCCTCGCGGTGATGGTGCTGCGCTGGAGTGACGGCAGTTATCTGGAAGATCAGGATATGTGGCGGATGAGCGGCATTTTCCGTGACGTCTCGTTGCTGCATAAACCGACTACACAAATCAGCGATTTCCATGTTGCCACTCGCTTTAATGATGATTTCAGCCGCGCTGTACTGGAGGCTGAAGTTCAGATGTGCGGCGAGTTGCGTGACTACCTACGGGTAACAGTTTCTTTATGGCAGGGTGAAACGCAGGTCGCCAGCGGCACCGCGCCTTTCGGCGGTGAAATTATCGATGAGCGTGGTGGTTATGCCGATCGCGTCACACTACGTCTGAACGTCGAAAACCCGAAACTGTGGAGCGCCGAAATCCCGAATCTCTATCGTGCGGTGGTTGAACTGCACACCGCCGACGGCACGCTGATTGAAGCAGAAGCCTGCGATGTCGGTTTCCGCGAGGTGCGGATTGAAAATGGTCTGCTGCTGCTGAACGGCAAGCCGTTGCTGATTCGAGGCGTTAACCGTCACGAGCATCATCCTCTGCATGGTCAGGTCATGGATGAGCAGACGATGGTGCAGGATATCCTGCTGATGAAGCAGAACAACTTTAACGCCGTGCGCTGTTCGCATTATCCGAACCATCCGCTGTGGTACACGCTGTGCGACCGCTACGGCCTGTATGTGGTGGATGAAGCCAATATTGAAACCCACGGCATGGTGCCAATGAATCGTCTGACCGATGATCCGCGCTGGCTACCGGCGATGAGCGAACGCGTAACGCGAATGGTGCAGCGCGATCGTAATCACCCGAGTGTGATCATCTGGTCGCTGGGGAATGAATCAGGCCACGGCGCTAATCACGACGCGCTGTATCGCTGGATCAAATCTGTCGATCCTTCCCGCCCGGTGCAGTATGAAGGCGGCGGAGCCGACACCACGGCCACCGATATTATTTGCCCGATGTACGCGCGCGTGGATGAAGACCAGCCCTTCCCGGCTGTGCCGAAATGGTCCATCAAAAAATGGCTTTCGCTACCTGGAGAGACGCGCCCGCTGATCCTTTGCGAATACGCCCACGCGATGGGTAACAGTCTTGGCGGTTTCGCTAAATACTGGCAGGCGTTTCGTCAGTATCCCCGTTTACAGGGCGGCTTCGTCTWD"
    .replace("WD", "GGATCC") # Ensures clean valid bases
    + "AAGCTTCGATCGATCGAATTCGATCGGCCATCGATCGATCGATCGATCGGATCCATCGATCGAATTC" * 4
)

# ==========================================
# 2. DIGESTION LOGIC
# ==========================================

def digest_dna(dna_seq, recognition_seq, cut_offset):
    """
    Finds all cleavage sites for a given enzyme on a linear DNA strand.
    Returns:
        cleavage_positions (list of 1-based indices where cuts occur)
        fragment_lengths (list of resulting fragment sizes in bp)
    """
    dna_seq = dna_seq.upper()
    recognition_seq = recognition_seq.upper()
    
    # Find all start indices of the recognition sequence
    cleavage_positions = []
    start = 0
    while True:
        idx = dna_seq.find(recognition_seq, start)
        if idx == -1:
            break
        cut_pos = idx + cut_offset
        cleavage_positions.append(cut_pos)
        start = idx + 1  # Advance search

    # Calculate fragment lengths based on cut positions
    cut_boundaries = [0] + cleavage_positions + [len(dna_seq)]
    fragment_lengths = [
        cut_boundaries[i + 1] - cut_boundaries[i]
        for i in range(len(cut_boundaries) - 1)
    ]
    
    return cleavage_positions, fragment_lengths

# ==========================================
# 3. OUTPUT: Console Report & Gel Simulation
# ==========================================

def run_analysis(dna_seq, enzymes):
    seq_len = len(dna_seq)
    print(f"Total DNA Sequence Length: {seq_len} bp\n")
    print("=" * 65)
    
    results = {}
    
    for name, (rec_seq, offset) in enzymes.items():
        positions, fragments = digest_dna(dna_seq, rec_seq, offset)
        num_cleavages = len(positions)
        results[name] = fragments
        
        print(f"Enzyme: {name} (5'-{rec_seq}-3', cut offset: {offset})")
        print(f"  1. Number of cleavages : {num_cleavages}")
        print(f"  2. Cleavage positions  : {positions if positions else 'None (uncut)'}")
        print(f"  3. Fragment lengths    : {sorted(fragments, reverse=True)} bp")
        print("-" * 65)
        
    plot_electrophoresis_gel(results, seq_len)

def plot_electrophoresis_gel(digest_results, total_length):
    """
    Simulates an agarose gel electrophoresis image with a DNA ladder
    and one lane per restriction enzyme.
    """
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Set dark gel background
    ax.set_facecolor("#121417")
    fig.patch.set_facecolor("#1e2229")
    
    # Define DNA Ladder markers (in bp)
    ladder = [5000, 3000, 2000, 1500, 1000, 750, 500, 250, 100]
    lanes = ["Ladder"] + list(digest_results.keys())
    num_lanes = len(lanes)
    
    # Migration physics: distance traveled is inversely proportional to log10(fragment_size)
    # Smaller fragments migrate further down (higher y-value if inverted, or lower y-value)
    def get_y_position(bp):
        bp = max(bp, 50)  # Prevent log(0)
        return np.log10(bp)

    # Draw wells at the top of each lane
    max_y = get_y_position(8000)
    min_y = get_y_position(50)
    
    for i, lane_name in enumerate(lanes):
        x_center = i + 1
        # Draw loading well
        ax.add_patch(plt.Rectangle((x_center - 0.25, max_y - 0.03), 0.5, 0.05, 
                                   color="#2c323d", ec="#555e6e", lw=1))
        
        if lane_name == "Ladder":
            fragments = ladder
            band_color = "#00e5ff"  # Cyan glow for ladder
        else:
            fragments = digest_results[lane_name]
            band_color = "#76ff03"  # Classic SYBR Green / EtBr fluorescence
            
        # Plot each DNA fragment band
        for frag in fragments:
            y_pos = get_y_position(frag)
            # Band intensity/thickness scales slightly with fragment size
            thickness = np.clip(2.0 + (frag / 1500.0), 2.0, 5.5)
            
            # Draw glowing band effect
            ax.hlines(y=y_pos, xmin=x_center - 0.28, xmax=x_center + 0.28,
                      color=band_color, linewidth=thickness, alpha=0.85)
            ax.hlines(y=y_pos, xmin=x_center - 0.28, xmax=x_center + 0.28,
                      color="white", linewidth=thickness * 0.3, alpha=0.6)
            
            # Label fragment sizes on the ladder and enzyme lanes
            if lane_name == "Ladder":
                ax.text(x_center - 0.35, y_pos, f"{frag} bp", color="#00e5ff",
                        va="center", ha="right", fontsize=8, fontweight="bold")
            else:
                ax.text(x_center + 0.31, y_pos, f"{frag}", color="#b2ff59",
                        va="center", ha="left", fontsize=7.5)

    # Formatting axes and labels
    ax.set_xlim(0.2, num_lanes + 0.8)
    ax.set_ylim(min_y, max_y + 0.1)
    ax.set_xticks(range(1, num_lanes + 1))
    ax.set_xticklabels(lanes, color="white", fontsize=11, fontweight="bold")
    ax.xaxis.tick_top()  # Put lane names at the top near the wells
    
    # Hide y-axis ticks since migration is logarithmic and labeled by the ladder
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor("#444c56")
        
    plt.title("Simulated Agarose Gel Electrophoresis (Restriction Digest)", 
              color="white", fontsize=13, fontweight="bold", pad=20)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # To use a downloaded NCBI FASTA file, uncomment the line below:
    # dna = load_fasta("ncbi_bacterial_seq.fasta")
    
    dna = SAMPLE_NCBI_DNA
    run_analysis(dna, RESTRICTION_ENZYMES)