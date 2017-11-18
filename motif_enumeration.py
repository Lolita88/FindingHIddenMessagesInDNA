

import copy,sys

pattern_neighbors = {}

def motif_enumeration(dna, k, d):
    # A brute force algorithm for solving the implanted motif problem.
    # It is based on the observation that any (k, d)-motif must be at most 
    # d mismatches apart from some k-mer appearing in the first string in Dna.
    # We use neighborhoods and the neighbors method to generate the motifs.
    # What we are looking for are the best motifs
    dna_patterns = []
    dna_length = len(dna[0])
    
    #save all patterns from first dna row in pattern_neighbors
    for i in range(0, dna_length-k + 1):
        dna_patterns.append(dna[0][i:i+k])
    #get neighborhood for first strand and save in dictionary
    for each in dna_patterns:
        neighborhood = neighbors(each, d)
        for each in neighborhood:
            pattern_neighbors[each] = pattern_to_number(each)

    #loop through remaining dna strands over and over to see if neighbors are found
    #if a neighbor is not found, stop checking for it and remove it from the pattern_neighbors dict
    for i in range(1, len(dna)):
        #if hamming distance too great - not found, remove from pattern_neighbors
        check_neighborhood_dict(pattern_neighbors, dna[i], d, k)
    for key, value in pattern_neighbors.items() :
        print(key)

def neighbors(pattern, d):
    if d == 0:
        return ([pattern])
    if len(pattern) == 1:
        return (['A', 'C', 'G', 'T'])
    neighborhood = []
    suffix_neighbors = neighbors(suffix(pattern), d)

    for text in (suffix_neighbors):
        if hamming_distance(suffix(pattern), text) <d:
            for base in "ACGT":
                neighborhood.append(base + text)
        else:
            neighborhood.append(first_symbol(pattern) + text)
    return (neighborhood)

def suffix(pattern):
    suffix = pattern[1:len(pattern)]
    return (suffix)

def hamming_distance(pattern1, pattern2):
    distance = 0
    for x, y in zip(pattern1, pattern2):
        if x != y :
            distance += 1
    return distance

def first_symbol(pattern):
    return pattern[0]

def pattern_to_number(pattern):
    k = len(pattern)
    if k == 0:
        return 0
    symbol = pattern[k-1]
    prefix = pattern[:-1]
    return 4 * pattern_to_number(prefix) + symbol_to_number(symbol)

def symbol_to_number(symbol):
    pass
    if symbol == 'A':
        return 0
    elif symbol == 'C':
        return 1
    elif symbol == 'G':
        return 2
    elif symbol == 'T':
        return 3

def check_neighborhood_dict(pattern_neighbors, dna, d, k):
    #see which pattern_neighbors are in at least one of each_dna_patterns
    #for each k-mer in the strand
    count = 0
    dna_length = len(dna)
    temp_pattern_neighbors = copy.deepcopy(pattern_neighbors)
    for each in temp_pattern_neighbors:
        count = 0
        for j in range(0, dna_length - k + 1):
            each_dna_pattern = dna[j:(j+k)]
            #if not pattern found, set counter to 1
            #only remove those without counter
            if(hamming_distance(each_dna_pattern, each) <= d):
                count = 1
        if count == 0:
            #print("removing " + each)
            pattern_neighbors.pop(each)

dna = ["ACGT", "ACGT", "ACGT"]

motif_enumeration(dna, 3, 0)

