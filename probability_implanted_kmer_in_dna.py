# Compute the probability that ten randomly selected 15-mers from the ten 600-nucleotide 
# long strings in the Subtle Motif Problem capture at least one implanted 15-mer. (Allowable error: 0.000001)

import math

dna_length = 600
kmer_length = 15
strands_dna = 10

def probability_implanted_kmer_in_dna(dna_length, kmer_length, strands_dna):
    # First you compute p1 - probability of not capturing the implanted k-mer (15-mer) in one string.
    # Then you notice for the entire problem we have to deal with ten similar cases, i.e. you have to
    # multiply p1 * p2... *p10, where p1 = p2 = ... = p10. So you just compute p1 to the 10th power:
    # Then you just compute the 'opposite' probability, i.e. the probability that from ten 600-length 
    # nucleotide string, we capture at least one implanted 15-mer! 
    
    
    #p = (dna_length - kmer_length) / ((dna_length - kmer_length) + 1)
    #return 1 - math.pow(p, strands_dna)
    
    return 1 -  ((600-15)/(600-15+1))**10 -  10*(600-15)**9/(600-15+1)**10
print(probability_implanted_kmer_in_dna(dna_length, kmer_length, strands_dna))
