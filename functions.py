
import time, math
import urllib2

response0 = urllib2.urlopen('http://bioinformaticsalgorithms.com/data/realdatasets/Rearrangements/E_coli.txt')
pageData1 = str(response0.read())
response0.close()

number_to_symbol = {0:'A', 1:'C', 2:'G', 3:'T'}
symbol_to_number = {'A':0, 'C':1, 'G':2, 'T':3}


def pattern_count(text, pattern):
    count = 0
    for i in range((len(text) - len(pattern))+1):
        if(text[i:i + len(pattern)] == pattern):
            count += 1
    return count

def pattern_index(text, pattern):
    pattern_positions = []
    for i in range((len(text) - len(pattern))+1):
        if(text[i:i + len(pattern)] == pattern):
            pattern_positions.append(i)
    return pattern_positions

def frequent_words(text, k):
    frequent_patterns = []
    count = []
    for i in range(0, len(text) - k + 1):
        pattern = text[i:i + k]
        count.append(pattern_count(text, pattern))
    maxCount = max(count)
    for i in range(0, len(text) - k + 1):
        if count[i] == maxCount:
            frequent_patterns.append(text[i:i + k])
    return remove_duplicates(frequent_patterns)  

def clump_finding(genome, k, t, L):
        #L = 50 is window we are sliding down
        #t = 4 is min times we're looking for kmer
        #k = 5 is len of kmer
        frequent_patterns = []
        clump = []
        temp_num = 4 ** k - 1
        count = 0
        for i in range(0, temp_num):
            clump.append(0)
        #print len(genome)
        #print (len(genome)) - L
        #print 4 ** k - 1
        for i in range(0, (len(genome)) - L):
            count += 1
            #print count
            text = genome[i:i + L]
            #print text
            frequency_array = computing_frequencies(text, k)
            #print frequency_array
            for r in range(0, 4 ** k - 1):
                if frequency_array[r] >= t:
                    clump[r] = 1
            #print clump
        
        for i in range(0, 4 ** k - 1):
            if clump[i] == 1:
                #print "found clump!!"
                pattern = number_to_pattern(i, k)
                frequent_patterns.append(pattern)
        return frequent_patterns
    #L is the window of length we're looking at
    #Output: All distinct k-mers forming (L, t)-clumps in Genome.
    #frequent_patterns = []
    #frequency_array = computing_frequencies(text, k)
    #print frequency_array

def better_clump_finding(genome, k, t, L):
    #L is window we are sliding down genome
    #t is min times we're looking for kmer
    #k is len of kmer
    frequent_patterns = []
    clump = []
    for i in range(0, 4 ** k):
        clump.append(0)
    text = genome[0:L]
    #instead of calling every time. call once
    frequency_array = computing_frequencies(text, k)
    for i in range(0, 4 ** k):
    # runs clump for first freq array only
        if(frequency_array[i] >= t):
        # if >= t nucleotide found in freq, add to clump
            clump[i] = 1
    my_genome_length = len(genome) - L + 1
    for i in range(1, my_genome_length): 
        #starts moving window down, starting at 1 through 25
        first_pattern = genome[i - 1:k + i -1]
        index = pattern_to_number(first_pattern)
        frequency_array[index] = frequency_array[index] - 1
        #last_pattern = genome[i + L - k:k] # changed this because it was giving a neg number
        #pseudocode meant to give a range for second number. Annoying as FUCK
        last_pattern = genome[i + L - k:i + L]
        index = pattern_to_number(last_pattern)
        frequency_array[index] = frequency_array[index] + 1
        if(frequency_array[index]) >= t:
            clump[index] = 1
    for i in range(0, 4 ** k):
        if(clump[i] == 1):
            pattern = number_to_pattern(i, k)
            frequent_patterns.append(pattern)
    #return len(frequent_patterns) #used for e coli genome
    return frequent_patterns
    #look up sparse arrays, only need one array, use index as pattern, log 4, read about bloom arrays


def remove_duplicates(some_list):
    new_list = []
    new_list.append(some_list[0])
    for i in some_list:
        if i not in new_list:
            new_list.append(i)
    return new_list

def computing_frequencies(text, k):
    frequency_array = []
    for i in range(0, 4 ** k):
        frequency_array.append(0)
    for i in range(0, len(text) - k + 1):
        pattern = text[i:i+k]
        j = pattern_to_number(pattern)
        frequency_array[j] += 1
    return frequency_array

def symbol_to_number(symbol):
    #print symbol
    if symbol == 'A':
        return 0
    elif symbol == 'C':
        return 1
    elif symbol == 'G':
        return 2
    elif symbol == 'T':
        return 3
    else:
        #return -1
        print 'not a nucleotide'

def number_to_symbol(number):
    if number == 0:
        return 'A'
    elif number == 1:
        return 'C'
    elif number == 2:
        return 'G'
    elif number == 3:
        return 'T'
    else:
        print 'not a nucleotide'


def pattern_to_number(pattern):
    k = len(pattern)
    if k == 0:
        return 0
    symbol = pattern[k-1]
    prefix = pattern[:-1]
    return 4 * pattern_to_number(prefix) + symbol_to_number(symbol)

def number_to_pattern(index, k):
    if k == 1:
        return number_to_symbol(index)
    prefix_index = math.floor(index/4)
    r = index%4
    symbol = number_to_symbol(r)
    prefix_pattern = number_to_pattern(prefix_index, k - 1)
    return prefix_pattern + symbol

start = time.time()    

#temp = better_clump_finding('CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA', 5, 4, 50)
temp = better_clump_finding('AAAACGTCGAAAAA', 2, 2, 4) #works by book way
#temp = better_clump_finding('ACGTACGT', 1, 2, 5)
#temp = better_clump_finding('CCACGCGGTGTACGCTGCAAAAAGCCTTGCTGAATCAAATAAGGTTCCAGCACATCCTCAATGGTTTCACGTTCTTCGCCAATGGCTGCCGCCAGGTTATCCAGACCTACAGGTCCACCAAAGAACTTATCGATTACCGCCAGCAACAATTTGCGGTCCATATAATCGAAACCTTCAGCATCGACATTCAACATATCCAGCG', 3, 3, 25)
#print better_clump_finding(pageData1, 9, 3, 500)
#print (" ".join(str(x) for x in temp)) 
#print pattern_count("CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC", "CGCG")
#print frequent_words("TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT", 3)
print pattern_index("ATGACTTCGCTGTTACGCGC", "CGC")

#print number_to_pattern(7812, 7)
#print pattern_to_number('ATGCAA')
#print (pattern_to_number2('CCCTGCGCCCAGGGG'))
#temp = computing_frequencies('CGGCAAGTCTGCTTAGGCAACGCGCACCCTGCAGGGGCTTGCGAGACCGCCACAATCGTGGAGGTGACGATTGGTCTCGAACGTGTATCTCGTAACGGCGCTTCCCCTCGCAGCTACTCGATCCCATAAAGCCTATGAGGGGCCGTTCCGCGCCGGAAACGGCGTTTATACGTGGACCAACGGTCACGACACCAAAGAGTGTACAGTCGCGCCTCAGTGGATCCGGTTCCCGACCGCGGCTAATTTTACATGGTCGTGGGTCATACGGCGTCATACTCGCGAAACATCCAGATGCGAAGAAGTCTGACCATCAGAAACGCTTCTTGGGTCTCTTCTGCCAGCTCAAAGGACGGACCCAACACCCGTCGTGCGGCCTGCCCTCATTGGCCGACAATTCCTTGATCTTCCTTACCAAGTGCTTTTTACCTATAGTTCATCGACTCAGACCGATCTGCCCCCTGCACGTTGTTACCGCAGGGAAACCTTGATGGCCCATAATAGCACATACCTTTTCGTGCATGATTGGGGACTTGAGGACATTGTTTAAAATCATCATCAACTAGGTATGCAAGGGATGAACGCTGCTTATATTGTTTCAGAACGTGAAGTAGCTCCG', 5) 
#print (" ".join(str(x) for x in temp))   
#print frequent_words('ACGTTGCATGTCGCATGATGCATGAGAGCT',4)
print 'It took ', time.time()-start, 'seconds.'
