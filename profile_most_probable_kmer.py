import sys

def profile_most_probable_kmer(text, k, profile):
    # Given a matrix Profile, we can evaluate the probability of every k-mer 
    # in a string Text and find a Profile-most probable k-mer in Text, i.e., 
    # a k-mer that was most likely to have been generated by Profile among 
    # all k-mers in Text.
    # If there are multiple Profile-most probable k-mers in Text, then we 
    # select the first such k-mer occurring in Text.
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

def profile(motifs):
    profile_matrix = Count(Motifs)
    #print "profileMatrix " + str(profileMatrix)
    t = len(motifs)
    for i in profile_matrix:
        for j in range(len(profile_matrix[i])):
            profile_matrix[i][j] = float(profile_matrix[i][j]) / float(t)
    return profile_matrix


"""profile = {

    'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
    'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
    'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
    'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
}"""
"""profile = {

    'A': [0.2,0.2,0.3,0.2,0.3],
    'C': [0.4,0.3,0.1,0.5,0.1],
    'G': [0.3,0.3,0.5,0.2,0.4],
    'T': [0.1,0.2,0.1,0.1,0.2]
}"""
profile = {

    'A': [0.349,0.229,0.181,0.277,0.217,0.229,0.265,0.193,0.241,0.361,0.289,0.193],
    'C': [0.205,0.265,0.325,0.253,0.289,0.193,0.265,0.217,0.205,0.217,0.253,0.181],
    'G': [0.205,0.265,0.229,0.253,0.229,0.265,0.193,0.386,0.301,0.169,0.241,0.398],
    'T': [0.241,0.241,0.265,0.217,0.265,0.313,0.277,0.205,0.253,0.253,0.217,0.229]
}
#print(pr("TCGTGGATTTCC", profile))
print(profile_most_probable_kmer("TCTTACAACCCAACCACGGGAGTAACTCCGACGCTGCTCCGGTACCATGACCGAAGGTCGCCCATTCGAGGCGGCTGCGTGAGGACTGTGAAAACACACTAATTGGCGCTGCCATACTGCGCAAACTTTCCTTCAAGAAGTAAATTAACTTAAACTCGCGTCGCGCGTACAAGATAGTATCCGGAGCCTCAAGCGCTGTTCTCAGTTCTCTCTCGTTCGCGTGATAGGAAGAAGTCGAACATCCGTTCCTCGAGAGCCAACATTAGTAAGAGGGGTAAAAGTAGTACATGTGCTGACGACCGACCCTACGTCAATCTAACCCAAGACAGAGGCATTGCCAATACATGTCTGACGTGTACGGGCGCTTGGCCTAGGTTTCGCTCCGTCGATCTGTCGATCTGGTGTAAGGGCACCCGCCTTACGGTTAGTTTCTTCTCACATGCTTGGGCAACATCTGTACCCTGCGCGGAATACGAACATGGCTCTCTATACCCCTTCCCGAGATTTCTAACCATAAAGACAAGGCCATAATCGTAGACCTAGGAACCGGTACCCTAAATCACCTCCGGTAAGATGCCCTTCGCAGTTCCGATCGGTGTAAGTGGGTAACGAAGACAGCTATAAAGCAGCTACGACAGTAGTTAGTCAAGCTGGACCGGGAGGATGCGGCTGCTTCAGACTTTATGTAGTAGCTGAGTACCTTGCTCTGTGCCCCTGACAGTATGTAAACGTCGGCTAAGATTCGAACTAGCTTTATGCAGAGCTGGGGAACGCATGGAGCTCTCCCCAGCACGAGGCTTAGCTGTACCGCTGACGACCATCCCTTTTACCTCGTTACGCGCTTTGGGAACTCATATAGAACAAAATCATACCAGCAATTCAATATTGTGGCCAGTGACAGGTGACCAATCTATGACCAGATGAGTCTGCTTCTCAAGGGCGTTGCTAGCGAGCTGTGAGCGAACGATCTATTGCGTCTCGCTAAGGCCGGTTGACCGTA", 12, profile))