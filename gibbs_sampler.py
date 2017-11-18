import random

def gibbs_sampler(dna, k, t, n):
    # form motifs by randomly selecting a kmer in each seq of dna
    # dna is a collection of dna strands
    # k is the length of kmer we are searching for
    # t is the number of strands of dna passed in
    # n is the number of iterations we are looping through the motifs, removing one, creating
    # a new motif, etc - needs better explaining
    motifs = random_motifs(dna, k, t)
    best_motifs = motifs
    best_score = score(best_motifs)
    for j in range(n):
        i = random.randint(0,t-1)
        #randomly choose one of the selected kmers and remove from motifs
        motifs.pop(i)
        #create profile matrix from remaining kmers in motifs
        profile = profile_with_pseudocounts(motifs)
        #for the missing kmer, calculate Pr(kmer/Profile)resulting in n-k+1 probabilities
        motifs_i = profile_generated_string(dna[i], profile, k)
        #roll a die(with n-k+1 sides) where probability of ending up at side
        #i is proportional to p of i.
        #Choose a new starting position based on rolling the die.
        #Add the kmer starting at this position to motifs
        #insert and compare, this is the algorithm
        motifs.insert(i, motifs_i)
        curr_score = score(motifs)
        if curr_score < best_score:
            best_motifs = motifs
    return best_motifs

def random_motifs(dna, k, t):
    #uses random.randint to choose a random k-mer from each 
    #of t different strings dna, and returns a list of t strings
    #random.randint(1, M)
    kmers = []
    for i in range(len(dna)):
        j = random.randint(0, len(dna[0])-k-1)
        kmers.append(dna[i][j:j+k])
    return kmers 

def score(motifs):
    # returns score, higher is worse I think, check
    this_consensus = consensus(motifs)
    score = 0
    for i in range (len(motifs)):
        for j in range (len(motifs[i])):
            if this_consensus[j] != motifs[i][j]:
                score += 1
    return score

def profile_with_pseudocounts(motifs):
    #takes a list of strings Motifs as input and returns the profile matrix of Motifs with pseudocounts as a dictionary of lists
    profile_matrix = count_with_pseudocounts(motifs)
    t = len(motifs)
    for i in profile_matrix:
        for j in range(len(profile_matrix[i])):
            profile_matrix[i][j] = float(profile_matrix[i][j]) / (float(t) + 4)
    return profile_matrix 

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

def profile_generated_string(text, profile, k):
    # takes a string Text, a profile matrix profile , and an integer k as input.  
    # It should then return a randomly generated k-mer from Text whose probabilities are generated from profile 
    n = len(text)
    probabilities = {}
    for i in range(0,n-k+1):
        probabilities[text[i:i+k]] = pr(text[i:i+k], profile)
    probabilities = normalize(probabilities)
    return weighted_die(probabilities)

def normalize(probabilities):
    #This function takes a dictionary Probabilities whose keys are k-mers and whose values are the probabilities 
    # of these k-mers (which do not necessarily sum to 1). The function should divide each value in Probabilities 
    # by the sum of all values in  Probabilities, then return the resulting dictionary.
    norm_kmers = {}
    kmer_sum = 0
    for k,v in probabilities.items():
        kmer_sum += v
    for k,v in probabilities.items():
        norm_kmers[k] = v/kmer_sum
    return norm_kmers

def weighted_die(probabilities):
    #takes a dictionary Probabilities whose keys are k-mers and whose values are the probabilities of 
    # these k-mers. The function should return a randomly chosen k-mer key with respect to the values 
    # in Probabilities.
    ran_num = random.uniform(0,1)
    count = 0
    #for each face of the die, add Pr(face) to the "count" variable. If "count" > the random number, then return the face
    for k,v in probabilities.items():
        count += v
        if count > ran_num:
            return k

def consensus(motifs):
    k = len(motifs[0])
    count = count_with_pseudocounts(motifs)
    this_consensus = ""
    for j in range(k):
        m = 0
        frequent_symbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequent_symbol = symbol
        this_consensus += frequent_symbol
    return this_consensus

def pr(text, profile):
    # returns probability number of text 
    p = 1
    for i in range(len(text)):
        #At position i of Text, we set p equal to p times the value of Profile 
        # corresponding to symbol Text[i] and column i, which is just Profile[Text[i]][i].
        p *= profile[text[i]][i]
    return p

dna = [
"GGATCCTAGGCTGTACATTGAGAACCTATAGATTTAATTGGCCCAACCTCGCCGCCCTCACCACAGCGGTATACCGATCGCACGCGGGTCGCACATACAGGTAATCAAATTCATGAACACGCACTCGACATCTGCTCAAGCAGTCTCCCAGCCCCCCGGTGTATATGGGGTAACACCTTCGATCACCACATGGTAATGGGGCACGTAGTTCCCATTCGAGGTACGGAACTTTCCCACAGACAGCCATCGGAGATCGCGTCCCAGAACATTCTTCTATTCGTTATTACCGGCCACCCAGGTTTAAACGCCTGTATGAGGGATCCTAGGCTGTA",
"CATTGAGAACCTATAGATTTAATTGGCCCAACCTCGCCGCCCTCACCACAGCGGTATACCGATCGCACGCGGGTCGCACATACAGGTAATCAAATTCATGAACACGCACTCGACATCTGCTCAAGCAGTCTCCCAGCCCCCCGGCCCAATCAATGGATGTGTATATGGGGTAACACCTTCGATCACCACATGGTAATGGGGCACGTAGTTCCCATTCGAGGTACGGAACTTTCCCACAGACAGCCATCGGAGATCGCGTCCCAGAACATTCTTCTATTCGTTATTACCGGCCACCCAGGTTTAAACGCCTGTATGAGGGATCCTAGGCTGTA",
"TACGTTCCCCTGACTTGGATGACATAGTCGGCCGACAGGTGTCCACAAGAGATTGAACTTACGCCGGCGGTTTCCTCGCATCTTCAGAATAATCCACGATGGTGGTCGAGATGCTTCTGCTATCAGGGAGTTGTGATGCGTTCCTTCCGCGGCCGACCGGGTTAGGGGACTACTTGACCTGTTAATCCCCGGACCCCCAATTACAAGTGTAAACGACGTCAATCAGGAGAGTCATAACTCTGACCCCTTTGACGATCAGTGTCTCCGCGGCACACGGCGCCTCTGTGAGCGTGGTGTATCAACTCGACCGACTACATTGTAGCTGGGCTTAA",
"CTCCACGACGTGCCCGCTTCGTCTAAAAACTTAAGCTTTTCCGGCGGAAGGAGCGTGCTTCCAATGAGCCAAGTTAAGTTAGTTAAGAACAGTAACTAGATGGCCCCGCGAATGGATGAGCCCCCACGTTCATGTCGAGTTCGGAAGCAAACACGGAGTTCAATGCGTCCAGCAAATTTACGCTTCTTATAATCCTAGCAGGCGGGTGTGCGAGATGCAACGGGGCTCTAAGGCTTAGTGTTGGACTTACGAATGCACGCACACCGTAGGCATTGCATTTCGGTGACTACGGGGTTTTTTCGAGTCAGGGTCCTAGATGAATAGTCATAATT",
"GAAAATTATTTAACTGGTAAACATATAGGGGAGTAAGCGAGGTCATTTCAAGATGGCAACTCATGGAGATGATCTGTTGTTGAAGGATCTAGGGCAAGGGAAGAACGTGCGTCGCGTTGTTCTTCCGGCCTGGGAGTGAAACGTGCCGACTACCCACAAAAGAGGATTCGGGTTTCTACGCATACGCGAACATCCGGTTCAACGGTAATGCAGATAGCACATAGTGTTACGGAAACCCAACCACGACATTTTAGGTGTGTCGGATTGTCTCACACTCGGTGATATCAGGCTATAACGAGCTTTAGATAGACCCCCTGCCGCGGATGCTGTCC",
"TTCTTATAGCTTAGTTAAACAACCTGGCCGAACGTCCCATCGATTTAGGCTCGGAGCGCGCTCAAGGCAGGGGACTGTTTAGCAGCACGCTAAAATGACAAGTCTTTACGTCTTGTAGCATTCTGAGTTAAGATTCTCTGTAAACCTACCCTGTATCAGAGCTCGAGCCTATGACGTAAAGGTGGACACACGGTGCGATCCAAAATGTATTTAGGTCAGTGGTTTTCCTCCTGTGCGCCAGAGGTGAGGACTATCTTCTCGCGACGTTGGTGTCCACAGTTGAGCTGCAATGGATGGTGTCTTCAGATCTCGAGAGATTATGCGAGCGGTGA",
"TTAGCCCCTGCAATATTTGAACATGAGCTGCTAACGTTCATCCAATGCGCTTTTGCCGTTGAGTTTCAAGGCAGAAGGTACAGGAGGTGCGCTGGTTGTAGGAAATGTAGGGCTGTGTTATCTGACGGTGTCCTTGTCCTATGTGTCTAGACCAACGATTAACAGGGGAAGCGAGTTGGTATATCCATTGGAGCGAAAATGGGATAAGGAACATCCCTATCGTAGCAAAAGCCGCCTTGTTCGATGACGCTATAGTTCATGTACGACGAGGGCTACTATTCGAGCTTCCCATAATCATACTGCCTGGCCTTCGCGATCTGGCTATAGGTAAG",
"AAGGGAGGTATACCCAAGGCATACCTGCGGCTTGGGTTGAAGGGTCATGTGATTTCACTCGCCAGAGTAGGGCCGGAGGTGAGCGCTAGTAGGCACGGCAGCCGTCAGTGGCCCGATGTTAAACCATCATTCTCGGAAAGTCTCCGAAGTTTCTCCCTTAGAAGAAATACCATACGTGTTTGCGATCACCATCACGCGGTCGTTTTAATATAATCTCACTGAACCCTGAGGTATGGATACACCATATAGCACCCCGTCTCACTCTTTGGTGGAAATCTACTTTCCCGCCAGGGGGTGCCCCCCTAAATGGATGCCTGAGTTAACAATTGACG",
"TGATATTCGGGGCACGGGTGGCTGTGTTTGCGCTCTAAACAGGCCCAAGCCAAAGACAAATGCAATGGATGCTTTTGCTCAAAGAGAGTAACGATCAGAGGATGCAGCATTTCTGAATGTCAGTTTGTAAACATCTTGTTTTGCGAACTGGGCCATACACAGAGCTTGTTCTCTACAGGGCATGCTAGGCAGTGTACTTAAGGGTCGCAGGCTCTTGCCAAGGTAACCAACAGATTAGAGAAGTCACGCAAAGTCATGTATACCTCCACTGTAGACATTAAACCCTTCAGATCTAGCGTACTCTGAAGTACATCGGTTCAAGGAGGATCGAC",
"CGACAGGTATCCGTAAAAGTTGGATAACAAACACTTACGAACCGTCTAGAGTATGTAAACTGGGAACGGGTTCATCACGTCCGATTACTGGCTCCGTAATTGCATATCTAATTAGCACTTTAATCGCCTCAGGGTCCGGCCTTCTGCTCCCGATCAATGGATGTAATGGCGCCCACCCTACAGCCTCTGTGATCCGGGTTGGACAAGGAACTCTACCCCTCGCTATGAGTACATAGAATACCGTAGGAGTCCCTAGATCATCCAGCGCCACATGTAGGACCTTGAGGTCTATAGCCGCACGACGCGGGGTCTTGGGTACCCCTTCAGGGCAA",
"CCTTTCCCGTTATTGGATATCAGGTGGCAGTCAGGGTTGGCTAAAAGCCCGGAGGAAACACACGGTGTAAAAACGCCCATCACTAATCTGCAAGATGAGTCAACCCGTCGACACTAAGGGGTGTGCCGGTACGCCACGATAGGTTGGGGTAGGGGGCCAAGCCACTCTCCTAGCCGAAAGTATATGGGACATGTTAGCCGATTAAAGCGTTGTTGTGAGGTACGCGCAGAGTCGTCCGCGTCTAACCATGGCAATGGATGGACGTCGGGGCCCAACTTCGTACCTCTGGGTTAAAACTTCGACGAGGTTTTCCTCCTGGCTCGGGATCCGTT",
"ACGAGAGAAGAATGAATCGGCTTACCCTCCGCATAATCCGAAAATAAGCCATACATCTCACTTTGCGATACCAAAAATCGACGCATGTACATCGGGATGTGATACCCTACCTAGGTATCCACGCTCTGATACAAACAGCCGAGAGTATATGCACACATATATGACTTCGGCGCATTTTGACGTCAGTCCCCTGCAAAAAATGCTACAGATAGAATGGTTGGTGCACCCTCACATGTCGCTCAGTTGTCCGCTCTGTCCTGAATAGAAGAATTGAGCAATCCCTGTCCGCGGCTTGACACAACCAGACGGGGACCCAACCCCGGTAGGCACTG",
"CTGCGCCGCACATCCACTAAGTAAAATGAATTACGCACCGCGTTCGGCAGAATGCGAAAGGCTAGATTGAACAAGATTTTCAGGGGTGCAACGAGCCTTTCAGCGCAATCCTGGGAGTGCGCCCGGCGGCCGACCTCGAATAGTTCTGATCTGTCTGCCCGGAATCAGGCACCCGGGCAGCAAATAGCTGAATGTGATATTCTCATCGATATCCGTTGACTGGCCGCTACCTCTTCACTTTTTGGGTACGGTTTATATTTGTGTTCAGCCTATGCCCCTTTGATGGATGAGGGCCTTTAACTTGCAACAGTGGATGGATTGCCCCTTCGAGC",
"GAGGATGCCCACGGCGGATGATATGGCTCCGCGATTACAGATATTGCCACTAACTCGCCCCTGCAATGGCCTTTACAAGCTTTTAGTAGAATCATTCTGCCTCTATCAGGTAACGACCGCCGGTCTGAACATCGCCAGCACCCTAGAAATGCCCTGGTCTGCCTTTCAGGGCGCCGCGACGACACTCGAATTCTTACGCAACTACGGAACTCATGGTGACTTGCAAAGTATCGGGCCCGCAAAATTCCCTGTCTATGCTCTAGCCCCCCCAAGCCGGTAATGAGCTTGCACATTTACGATCTCTAAGCTTGAGTTTGCGACAAGGGCAAAAT",
"TCCACGGAGGCGTCTCTCCGAGACCTACTTACGGACGACTAGTCTCAATTGCATCGAGCGACCAGACCCGGGCTGAACGAATAGCGTCCTGTTTTGATAACCTTATCTAGTAGTCTGCTGCGCCCCTGCAATGTCGGCTCCTTTCAAAGGGAGTTAGAGCATACCCCTCGTTGCCCCTTACCCCCTGAGATATCGCGATGAGTAGCTAAGACACTCCCCCGCCCTTGTGCTTCCATAGTATTTATTGTCGAGGAGTATCTGTAACATGAGACCTAGAAGCGGACGCTTGCTGGCACGATACCAAGAGGTAGTATAACGTAACTTAGGCCTGC",
"CTTTACACCCAAACGGCCGGTACATCTAAAATAAGGACGCAAGCACAGGGAACAGCCCATGTTGACCCTGACGAGATCACCCCCTCAAGGGAATTCAACGAGAGGAAGGGCTTAATGAGAGAACTGTTCACGTGCACCGTCCGCTAGCTGTATGATTTGGACCATAGTCATAGCATGTTGATTGCTGAGTTACAATCCTTTCCCGTGCGTAATGAAGCGGGGAGTGCACATACTGAAAATGGGTTGAGCTGTGGTCATTATAAATCGCACGCTTGTAGAGAGCCCCCTAACCTGCAATGGATAGGCGGCGATTGGACCAAACCAGTACTTTC",
"GGCGATCTCTGGAGACAGTATCGATTTCTTTACTGGCCGTCCGAACAGCAGCGGGCGATGAAATTTCGGCTCTGTTTATATGGTTAGTAAGCCTCCAGTGCTATTGTACTAGAGCCCCCCCCCTGTCGTGGATGGTTAACAAATGTGCGTAGTTCGTTTGATCCCGTTACCCTAGAGCTGCCTTAGCTCCGAGATGGGAAATCCTGAGCCATTAACTAAATGGAGACGAAGTTTATTCGGTAGCAGACACTTTGACCAAAAGGCTGTCAACAGCCAAACACACGGAGGGACCCGTATGGAGGATGCGTTTTGAGCTGACCGTCTAGTTTCGA",
"CCCCTAGGATGGATGTAGCATCAGTACGGTGCGGGCATGCCAGCTCTAACGTTGTAACGATCTGCTGCGATCAGCACGCCTATGTTACTATTGTCATTCCAAACGACTTATTTTGTTAGGCCGCCGCAGCGACCAAGGATAATTGATTGCGGGACTGTCTAGCTCAAGATTATTGCGGACCCCTTGTATTAGAGGTGAAATGGTACGCTTCAACTCATCATCACTATGGGGAGCAGTGTGTTCATGGCAGCCGTACTCATGAAGCAGAGCGTTTAAGGGGCCCTTTCAGGAGCGTCATCGCTGACGTGCATCATACCGCCCTTATCTCAAAG",
"CACTGTAAGCAGGATCCCAGTACTGGTTCTAGTGGTTTCCCTGTAGACACACAGGGCTACCTCAATGCTATTACAGACATGCAATTGACGCCTTAATTAACCCGCCTCCCTCGTGTGCACCTCCCAAAGCGACACTTCGAAACAAGTAAGTCAGCCAAAGCCCCGTCCGTATCTGAAGCGTCGTGCACGCGGTACCCCCCAGTTGTGTTGGATCTAGATGCAGGAAACCGGGAGTGCCCATATTGCTCTTACAGCCTTCCTGATTAGTTCACTCTAACATTGTCCCAAGCCACCCTGCAATGGAGTCACGGGGTTAGTAGAGCGCCTCCGAC",
"CGATGGCTGTTGATATAACCACCCCTGCAGATGATGCAATTTAACTAATTAACTTGGTTGACGTCTTTTGGTCCGGGGGGATCTTGTTGGCGACTGTCGCCAAGTTGATCATGATTATGTTAGCCTTGAACCGTTTCCCTCTTACCGCGTGAATATTCACCTGGATGGCGCTGGCTAGGCCGTCCCATACCTGAGAGAACGTTGGCTCCAATATCACGGTAATCCTGGGGATTACAACTGAGCCCAATCGATACCGCACAAGATCCGTATTACTATGCTAACGAGTGAAGTTGAAACAGAATGGCGTCAGGTAATGGAGTCTCACTGCATCG"
]
k = 15
t = 20
n = 2000

#print(gibbs_sampler(dna, k, t, n))
#print('\n'.join(gibbs_sampler(dna, k, t, n)))

#Calling gibbs_samplers 20 times in a row and passing the data back in
count = 0
best_motifs = []
best_motifs = gibbs_sampler(dna, k, t, n)
best_score = score(best_motifs)
while count < 20:
    curr_motif = gibbs_sampler(dna, k, t, n) 
    if best_score > score(curr_motif):
        best_motifs = curr_motif
        best_score = score(curr_motif)
    #print(curr_best_motif)  
    count += 1

#print(best_score)
print('\n'.join(best_motifs))
