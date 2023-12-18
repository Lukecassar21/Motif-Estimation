from Bio import SeqIO
from collections import Counter
import numpy as np
import pandas as pd

def find_motif(file):
    # Initialize a Counter object
    kmer_counts = Counter()

    # Parse the fasta file
    sequences = SeqIO.parse(file, "fasta")

    # Loop over each sequence in the file
    for seq in sequences:
        # Convert the sequence to a string
        sequence = str(seq.seq)

        # Generate all k-mers for the sequence
        kmers = [sequence[i:i+10] for i in range(len(sequence) - 9)]

        # Update the counts for each k-mer
        kmer_counts.update(kmers)

    # Find the most common k-mer
    motif = kmer_counts.most_common(1)[0][0]

    return motif

# Use the function
file = "C:/Users/Luke/PycharmProjects/Motif Estimator/data/reads_motif.fa"
motif = find_motif(file)
print(f"The most common motif is: {motif}")

def find_best_kmers(file, motif):
    # Initialize a list to store the best k-mers
    best_kmers = []

    # Parse the fasta file
    sequences = SeqIO.parse(file, "fasta")

    # Loop over each sequence in the file
    for seq in sequences:
        # Convert the sequence to a string
        sequence = str(seq.seq)

        # Initialize the best score and best k-mer for this sequence
        best_score = -1
        best_kmer = None

        # Generate all k-mers for the sequence
        kmers = [sequence[i:i+10] for i in range(len(sequence) - 9)]

        # Loop over each k-mer
        for kmer in kmers:
            # Calculate the score for this k-mer
            score = sum(a == b for a, b in zip(kmer, motif))

            # If this score is better than the best score, update the best score and best k-mer
            if score > best_score:
                best_score = score
                best_kmer = kmer

        # Add the best k-mer for this sequence to the list
        best_kmers.append(best_kmer)

        # If we have found 10,000 best k-mers, stop
        if len(best_kmers) == 10000:
            break

    return best_kmers

# Use the function
file = "C:/Users/Luke/PycharmProjects/Motif Estimator/data/reads_motif.fa"
motif = "ATTGAGTTTC"
best_kmers = find_best_kmers(file, motif)
print(f"The best k-mers are: {best_kmers}")

def create_pfm(best_kmers, motif):
    # Initialize a matrix of zeros
    pfm = np.zeros((4, len(motif)))

    # Define a mapping from nucleotides to indices
    nuc_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    # Loop over each k-mer in the list
    for kmer in best_kmers:
        # Increment the count for each nucleotide at each position
        for i in range(len(kmer)):
            nuc = kmer[i]
            pfm[nuc_to_index[nuc], i] += 1

    # Normalize the PFM by dividing each count by the total number of k-mers
    pfm /= len(best_kmers)

    # Create a DataFrame from the PFM
    pfm_df = pd.DataFrame(pfm, index=['A', 'C', 'G', 'T'], columns=range(1, len(motif) + 1))

    #Write the DataFrame to a CSV File
    pfm_df.to_csv("pfm.csv")

    return pfm_df

# Use the function
pfm_df = create_pfm(best_kmers, motif)
print(f"The positional frequency matrix is:\n{pfm_df}")