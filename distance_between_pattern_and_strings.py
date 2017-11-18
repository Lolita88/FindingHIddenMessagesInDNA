""" Returns the hamming distance between a pattern passed in and all patterns
    of that length in dna strands passed in."""

import math

def distance_between_pattern_and_strings(pattern, dna):
    k = len(pattern)
    distance = 0
    for row in dna:
        this_hamming_distance = math.inf
        for i in range(len(row) - k + 1):
            temp_pattern = row[i:i+k]
            if this_hamming_distance > hamming_distance(pattern, temp_pattern):
                this_hamming_distance = hamming_distance(pattern, temp_pattern)
        distance = distance + this_hamming_distance
    return distance

def hamming_distance(pattern1, pattern2):
    distance = 0
    for x, y in zip(pattern1, pattern2):
        if x != y :
            distance += 1
    return distance

#got tired of copying in and formatting the data by hand. This is faster - putting into an array, dna.
dna_unparsed = "ACACGACAATTAGGATTGCCCAATTCATTGTATAGCAAAATATATTCCCCTAAGCTTCACCGAAGATATCCAATGTCTCATAGGCGGTACGCAGCG TAGAAATTCGGTTGCAAACGGACTTAGATTCATGCCCACTTCAGCCAATATCCTACGGTAGAGAAAGTCCCCCATTACGCTACAGCTGGCAATTGG CTAAGAGAGACCATAGATCTAGGTTTAAGGTCGGCGGTCTCTATCGCCGAATCCTGCTGAACCAGTCTAGAGTCTATTTGGAGAATTATAGACGAA ATGTGCGTTGGTCGGTTGGTGGAGTAACAACAGGCGGCAACTCCTCCGGATTCATGTCGGCTTGCCTCCATTTATAACGCGTGTTAGCCCTATGTC TGGTCTGACAAGTTATTGATCCTGAGACTAACTTCTCAGTCACGGAACATTCCCGCCCGCCCGAAGTCGGATTTATGATCACAGCTGAGACAGCGG CGTTCGCATCACTCTGCGTGTGACGTTCAGTTGTTACTAGTGCACACTGACTCTCCTCCGCTGCGTCATAAGAAACTGATCCTTGTGACTAGCTCT AATCGACAGCGTTAGAGCACCACGGACTATGACTAGGAGCTCTAGCGCATTCAATCGTACCATAGATGGCGTCAACTAGCGATTGGGTTCTATGGC AAACATTTTGGTAAGTACACAGCTCCCCTAATCCCAATCTGCTATTGGAGTTTACTAAGAGTCGTGAAACTCTTCATGGCCAAGCCCCATCGGACT CACACTTTGCTCGCTGTCGTAGCAACTCTGAACAGTGCGGTTGATCCGTCGCCAGCCAACCTATGGTTGGTATCTAGTCATTTATCCTATTGGGTT TGTGTTATAAGCTCTTGAATATTACCCGCTTTAGGACTCCCGTGAATTCTAAAGCATATCACGTAGGAGCAGCACGAGCGAAAAGGGCGTATATGA CATAGCCTTTTACCGAACGCATGCGCCAACGAACGTCCTCTGGGCCCTGCACTGAAACTATCTGAGAGTGTATTCTGGAACCGGGATATGTGATAA GGATTAAGTTAGAAGAAAAAGCCGAGGGCACATAAGAGCCTTCCAACTTAAAGCGTAGTATCCGACACGGCCGAGATCCTCATCTCCGCCCTTAGC TCAAGAGAACGCTTCTGCCTCTGGTGTTGTCCAGTACGTGAGAGACTTAAACAAATAGCGCTTGGGCATTCCGCCACGTTCTCTTAGTCACGGAGC TGTAAGGCGGAGGTTCTACTTTGCTATACCAGTATGGGCCCGGGCAGGCAGCCAAAAGCACCAGATTCTATCTGTAACGCCATGCACCGAAGCCGC TACACGTTCCAGCATTCGTTTCGGCGAGACTTTTCCTTTGGAGACATAAGTACCAAGCTCACCACTCACCCCGAGTGATTGGTCTGACGGCAGATC AAACGCTCGCTGCACTATGCTGTCCAGTATTTGCGTCCTTGCCAGATCCTGCGTTCGGCTACAGTACTCACCGAAACCGTATACGTGGTTGGATCC CACGTGGATCCACGGAAGGTACTCTTGGGAGCGGGCGCAAGGGGGTGGTATCACGGAGTCGGCTATGTGATTAAGAGTTGATCCTTCAGGAATTGA AAGTACCTCTGCGGCGCGAGACAATTAGAAGCGTAACTAACACAATAGGACTAAAATACTATCTTGCGCTGTCCCAACGCAGCTTCGACATTTGGT GTTGTCTCAAGAACGAACCTCGCAAAGTTAGGTCGATCCAACTCTCGGATGCTATGCTGATACGATGCATACGTTGCCAACGGGCTTTAGACATAG TCCCGGGGCAGAAGATCCAAGTAACGGACTGTAGCCGGATACGTCCTCTGTCGCGGATAGGCCTTGTATAGCGCCCCCGAATACCTACTAATCGGT CCATGACTACTTAGAAAGCCCCGATCAGTCGTGGCTCGCAAGCGTCACCAGGATGGGCACACATGCAGCAATTCTAAACCTAACCTTTCTGCAGAA GAAAGGACGATACTCTGTTGCCCCGCCTAAGTCTTTTAAACGAAGTTCCAGAGGTATATACTTGGAATAACTACTTAGCGCCGGTAATAATCAGGC TGACAAACACGCTCTTATAAGGGAGAGTGAGGAGCATTATTACTATGACGACAACACTCCAAGCTATGTCAAGATCGCCTTAACTCAAATGAAGGT GCCCAACCACTTGGCAAAGAGGCCCCTACTCGAAGGTTGAAACGATCGAACGTATAGCGCTGATCCCTCAAGCCTTACGTGTTAGTCAACCAGGCT CTTCTCGCGTTATTCCTGGTCTTTGTTTCTGCTCACATGAGACGTATACTCACACTGCGACGATGATTGAATCTATCATCAGCGGCCCTTACATGG ACGCCATGGGGACTGTGTTGAACCTTACGGAATCATTCCGGTTGGTTCGATATTGAGCCACACGTGACAGAATTTCTTGACGGTATAAGGCCGATA AGCATTTGGGCAATATCAAGCCTGGTGGCAAACCCACGCAGACGTTATTTGGAAATATCGAAGTGTGCCAAACGCGAGTGGCGTTTGCGGAATCCT TTCAAAACTACCGTCTATGATAGCATAAAAATGGTATACGCGGGAGCGTTCCTTTAAGTCTGTTTGTGGTACATCAGAATTCCAGCCGTTTTATCT GCAGGAAACATTCCGGCTTAGGCTGCCCATTTAGCAACAGACTACCCTGTCTTCCTGATATAGTACAACGAATCATAGGTCCCATGTTGGAATCGC TATCTACTACGCTTCACAAATCGAATTAAATCGAGCTTGTTTCAGCCCTTTTTAGTGCTGATGATGAAACTGTCTTCTGGTTGGTTAGCGAAGCAA AAATCTCCCTACCAAAGAACAGGATAGACTCTAAACCACAAACAGGACTATAGCGTTACCAGAGACGAGACTCAGTCAGAACAGCATGCTGGGCCA"
dna_unparsed = dna_unparsed.replace(" ", "\n")
dna = dna_unparsed.splitlines()

pattern = "ACAAA"

print(distance_between_pattern_and_strings(pattern, dna))