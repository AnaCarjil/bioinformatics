import random
import sys

N_SAMPLES = 20000
MIN_LEN = 100
MAX_LEN = 150


def read_fasta(path):
    """Return the genome as one uppercase string (header lines are skipped)."""
    sequence = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if sequence:      # a second header means a second record: stop
                    break
                continue
            sequence.append(line.upper())
    return "".join(sequence)


def sample_genome(genome, n_samples, min_len, max_len):
    """Return a list of (start, sample) tuples taken at random positions."""
    samples = []
    for _ in range(n_samples):
        length = random.randint(min_len, max_len)            # 100..150 inclusive
        start = random.randint(0, len(genome) - length)      # sample stays inside the genome
        samples.append((start, genome[start:start + length]))
    return samples


def write_fasta(samples, path):
    with open(path, "w") as f:
        for i, (start, seq) in enumerate(samples, start=1):
            f.write(f">sample_{i} start={start} length={len(seq)}\n")
            f.write(seq + "\n")



def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "genome.fasta"
    genome = read_fasta(path)
    if len(genome) < MAX_LEN:
        sys.exit("Genome is shorter than the maximum sample length.")

    print(f"File: {path}")
    print(f"Genome length: {len(genome):,} bases")

    samples = sample_genome(genome, N_SAMPLES, MIN_LEN, MAX_LEN)
    write_fasta(samples, "samples.fasta")

    lengths = [len(s) for _, s in samples]
    total_bases = sum(lengths)
    print(f"Samples taken: {len(samples):,}")
    print(f"Shortest sample: {min(lengths)}   Longest sample: {max(lengths)}")
    print(f"Average length: {total_bases / len(samples):.1f}")
    print(f"Average coverage: {total_bases / len(genome):.1f}x")
    print("Samples saved to samples.fasta")


if __name__ == "__main__":
    main()