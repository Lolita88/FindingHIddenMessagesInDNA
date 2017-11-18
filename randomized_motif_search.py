import random

def randomized_motif_search(dna, k, t):
    m = random_motifs(dna, k, t)
    best_motifs = m
    #m = ["TGA","GTT","GAA","TGT"]
    #best_motifs = m
    while True:
        profile = profile_with_pseudocounts(m)
        m = motifs(profile, dna)
        if score(m) < score(best_motifs):
            best_motifs = m
        else:
            return best_motifs 

def random_motifs(dna, k, t):
    #uses random.randint to choose a random k-mer from each 
    #of t different strings Dna, and returns a list of t strings
    #random.randint(1, M)
    kmers = []
    for i in range(len(dna)):
        j = random.randint(0, len(dna[0])-k-1)
        kmers.append(dna[i][j:j+k])
    return kmers

def profile_with_pseudocounts(motifs):
    #takes a list of strings Motifs as input and returns the profile matrix of Motifs with pseudocounts as a dictionary of lists
    profile_matrix = count_with_pseudocounts(motifs)
    t = len(motifs)
    for i in profile_matrix:
        for j in range(len(profile_matrix[i])):
            profile_matrix[i][j] = float(profile_matrix[i][j]) / (float(t) + 4)
    return profile_matrix  

def motifs(profile, dna):
    kmers = []
    #takes a profile matrix Profile corresponding to a list of strings Dna as input 
    #and returns a list of the Profile-most probable k-mers in each string from Dna
    k = len(profile["A"]) 
    for i in range(len(dna)):
        kmers.append(profile_most_probable_pattern(dna[i], k, profile))
    return kmers

def profile_most_probable_pattern(text, k, profile):
    best_pr = -1
    best_pattern = ""
    k = int(k)
    for i in range(len(text)-k+1): #loop through text until almost end, no index error
        pattern = text[i:i+k] #loop through pattern
        if pr(pattern, profile) > best_pr: #call Pr on each Pattern and Profile, if > set
            best_pr=pr(pattern, profile) #best probability
            best_pattern=pattern #update pattern with new match
    return best_pattern

def pr(text, profile):
    # returns probability number of text 
    p = 1
    for i in range(len(text)):
        #At position i of Text, we set p equal to p times the value of Profile 
        # corresponding to symbol Text[i] and column i, which is just Profile[Text[i]][i].
        p *= profile[text[i]][i]
        #print(p)
    return p

def count_with_pseudocounts(motifs): 
    #returns the count matrix of Motifs with pseudocounts as a dictionary of lists
    #same as count but offset using pseudocounts to handle zero values
    count = {}
    k = len(motifs[0])
    # We then range over all nucleotides symbol and create a list of zeroes corresponding to count[symbol]
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)
    # six zeroes under each nucleotide symbol
    # Range over all elements symbol = Motifs[i][j] of the count matrix and add 1 to count[symbol][j].
    t = len(motifs)
    for i in range(t): #loops through each row of data
        for j in range(k): #loops through each column
            symbol = motifs[i][j]
            count[symbol][j] += 1 #saves count of each nucleotide per column
    return count

def consensus(motifs):
    k = len(motifs[0])
    count = count_with_pseudocounts(motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequent_symbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequent_symbol = symbol
        consensus += frequent_symbol
    return consensus

def score(motifs):
    this_consensus = consensus(motifs)
    score = 0
    for i in range (len(motifs)):
        for j in range (len(motifs[i])):
            if this_consensus[j] != motifs[i][j]:
                score += 1
    return score

"""
#dataset 0
dna = [
"CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", 
"GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", 
"TAGTACCGAGACCGAAAGAAGTATACAGGCGT", 
"TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", 
"AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]"""
"""
#dataset 1
dna = [
"AATTGGCACATCATTATCGATAACGATTCGCCGCATTGCC",
"GGTTAACATCGAATAACTGACACCTGCTCTGGCACCGCTC",
"AATTGGCGGCGGTATAGCCAGATAGTGCCAATAATTTCCT",
"GGTTAATGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG",
"AATTGGACGGCAACTACGGTTACAACGCAGCAAGAATATT",
"GGTTAACTGTTGTTGCTAACACCGTTAAGCGACGGCAACT",
"AATTGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG",
"GGTTAAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA"
]
"""
"""
#dataset 2
dna = [
"GCACATCATTAAACGATTCGCCGCATTGCCTCGATTAACC",
"TCATAACTGACACCTGCTCTGGCACCGCTCATCCAAGGCC",
"AAGCGGGTATAGCCAGATAGTGCCAATAATTTCCTTAACC",
"AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAAGGCC",
"AACCGGACGGCAACTACGGTTACAACGCAGCAAGTTAACC",
"AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGAAGGCC",
"AAGCTTCCAACATCGTCTTGGCATCTCGGTGTGTTTAACC",
"AATTGAACATCTTACTCTTTTCGCTTTCAAAAAAAAGGCC"
]
"""
"""dna = [
"TGACGTTC",
"TAAGAGTT",
"GGACGAAA",
"CTGTTCGC",
]
"""
dna = [
"CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
"GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
"TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
"TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
"AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
]
k = 8
t = 5

count = 0
best_motifs = randomized_motif_search(dna, k, t)
best_score = score(best_motifs)
while count < 1000:
    cur_motif = randomized_motif_search(dna, k, t)
    if best_score > score(cur_motif):
        best_motifs = cur_motif
        best_score = score(cur_motif)
    count += 1

print('\n'.join(best_motifs))

#print('\n'.join(randomized_motif_search(dna, k, t)))