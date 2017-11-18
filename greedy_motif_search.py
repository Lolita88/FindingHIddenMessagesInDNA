def greedy_motif_search(dna, k, t):
    """forming a motif matrix from arbitrarily selected k-mers in each string 
    from Dna (which in our specific implementation is the first k-mer in each 
    string). It then attempts to improve this initial motif matrix by trying 
    each of the k-mers in Dna1 as the first motif. For a given choice of k-mer 
    Motif1 in Dna1, it builds a profile matrix Profile for this lone k-mer, and 
    sets Motif2 equal to the Profile-most probable k-mer in Dna2. It then iterates 
    by updating Profile as the profile matrix formed from Motif1 and Motif2, and 
    sets Motif3 equal to the Profile-most probable k-mer in Dna3. In general, 
    after finding i − 1 k-mers Motifs in the first i − 1 strings of Dna, 
    GreedyMotifSearch constructs Profile(Motifs) and selects the Profile-most 
    probable k-mer from Dnai based on this profile matrix. After obtaining a 
    k-mer from each string to obtain a collection Motifs, GreedyMotifSearch 
    tests to see whether Motifs outscores the current best scoring collection 
    of motifs and then moves Motif1 one symbol over in Dna1, beginning the entire 
    process of generating Motifs again."""
    count = 1
    best_motifs = []
    for i in range(0, t):
        best_motifs.append(dna[i][0:k])
    n = len(dna[0])
    for i in range(n-k+1):
        motifs = []
        motifs.append(dna[0][i:i+k])
        for j in range(1, t):
            p = profile(motifs[0:j])
            motifs.append(profile_most_probable_pattern(dna[j], k, p))
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    temp = best_motifs[0]        
    for i in range(len(best_motifs)):
        if best_motifs[i] != temp:
            temp = best_motifs[i]
            count += 1  
    #return best_motifs, count
    return best_motifs

def profile(motifs):
    profile_matrix = count(motifs)
    t = len(motifs)
    for i in profile_matrix:
        for j in range(len(profile_matrix[i])):
            profile_matrix[i][j] = float(profile_matrix[i][j]) / float(t)
    return profile_matrix

def score(motifs):
    this_consensus = consensus(motifs)
    score = 0
    for i in range (len(motifs)):
        for j in range (len(motifs[i])):
            if this_consensus[j] != motifs[i][j]:
                score += 1
    return score

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
def count(motifs):
    count = {}
    k = len(motifs[0])
    # range over all nucleotides symbol and create a list of zeroes corresponding to count[symbol]
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    # six zeroes under each nucleotide symbol
    # Range over all elements symbol = Motifs[i][j] of the count matrix and add 1 to count[symbol][j].
    t = len(motifs)
    for i in range(t): #loops through each row of data
        for j in range(k): #loops through each column
            symbol = motifs[i][j]
            count[symbol][j] += 1 #saves count of each nucleotide per column 
    return count

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

def pr(Text, Profile):
    # returns probability number of text 
    p = 1
    for i in range(len(Text)):
        #At position i of Text, we set p equal to p times the value of Profile 
        # corresponding to symbol Text[i] and column i, which is just Profile[Text[i]][i].
        p *= Profile[Text[i]][i]
    return p

dna = ["TTACCTTAAC","AGGATCTGTC","CCGACGTTAG","CAGCAAGGTG","CACCTGAGCT"]
temp = greedy_motif_search(dna, 4, 5)
print('\n'.join(temp))
#print(greedy_motif_search(dna, 3, 5))

#sample url for test data
#https://stepik.org/api/attempts/25309005/file

#print(pr("CAGTGA", {'A': [0.4,0.3,0.0,0.1,0.0,0.9], 'C': [0.2,0.3,0.0,0.4,0.0,0.1], 'G': [0.1,0.3,1.0,0.1,0.5,0.0], 'T': [0.3,0.1,0.0,0.4,0.5,0.0]}))
