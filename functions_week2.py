import time, math, sys, urllib.request, string
import sys #if using stdin
#from string import rstrip

response = urllib.request.urlopen('http://bioinformaticsalgorithms.com/data/Salmonella_enterica.txt')

encoding = response.headers.get_content_charset('utf-8')

lines = []
for line in response:
    lines.append(str(line.decode(encoding).splitlines()))

lines.pop(0)

test_genome = (" ".join(str(x) for x in lines))
test_genome = test_genome.replace("']", "")
test_genome = test_genome.replace("['", "")
test_genome = test_genome.replace(" ", "")
#print(test_genome)
response.close()


def skew(genome):
    skew_array = {}
    skew_array[0] = 0
    minimun_skew = 0
    minimum_skew_array = []
    for i in range(1, len(genome) + 1):
        if genome[i-1] == 'G':
            skew_array[i] = skew_array[i-1] + 1
        elif genome[i-1] == 'C':
            skew_array[i] = skew_array[i-1] - 1
        else:
            skew_array[i] = skew_array[i-1]
    return skew_array
    #return ' '.join([str(skew_array[i]) for i in sorted(skew_array.keys())])

def minimum_skew(genome):
    temp_skew = skew(genome)
    min_value = min(temp_skew.values())
    min_keys = [k for k in temp_skew if temp_skew[k] == min_value]
    #return min_keys
    return ' '.join(str(x) for x in min_keys)

def hamming_distance(pattern1, pattern2):
    distance = 0
    for x, y in zip(pattern1, pattern2):
        if x != y :
            distance += 1
    return distance

def approximate_pattern_matching(pattern, genome, d):
    positions = []
    if not(d.isdigit()):
        d = 0
    else:
        d = int(d)
    for i in range(0, len(genome) - len(pattern) + 1):
        if(hamming_distance(pattern, genome[i:i + len(pattern)]) <= d):
            positions.append(i)
    #return positions
    return ' '.join(str(x) for x in positions)

def approximate_pattern_count(genome, pattern, d):
    count = 0
    for i in range(len(genome) - len(pattern) + 1):
        if(hamming_distance(pattern, genome[i:i + len(pattern)]) <= d):
            count += 1
    return count

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

def number_to_symbol(number):
    if number == 0:
        return 'A'
    elif number == 1:
        return 'C'
    elif number == 2:
        return 'G'
    elif number == 3:
        return 'T'

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

def reverse_complement(pattern):
    complement = []

    for nucleotide in pattern:
        if nucleotide == 'A':
            nucleotide = 'T'
        elif nucleotide == 'T':
            nucleotide = 'A'
        elif nucleotide == 'G':
            nucleotide = 'C'
        else:
            nucleotide = 'G'

        complement.append(str(nucleotide))

    complement.reverse()

    return ''.join(complement)

def get_window(genome, base, length):
    #new_base = int(base - length/2)
    #print(new_base)
    #print(int(base - length/2))
    #print(base + int(length/2))
    return(genome[int(base - length/2):base + int(length/2)])
    #return(genome[base-500:base+502])

def frequent_words_with_mismatches(genome, k, d):
    frequent_patterns = []
    frequency_array = {}
    close = {}
    for i in range(0, 4 ** k - 1):
        close[i] = 0
        frequency_array[i] = 0
    for i in range(0, len(genome) - k + 1):
        neighborhood = neighbors(genome[i:i +  k], d)
        for pattern in neighborhood:
            index = pattern_to_number(pattern)
            close[index] = 1
    for i in range(0, 4 ** k -1):
        if close[i] == 1:
            pattern = number_to_pattern(i, k)
            frequency_array[i] = approximate_pattern_count(genome, pattern, d)
    max_count = max(frequency_array.values())
    for i in range(0, 4 ** k -1):
        if frequency_array[i] == max_count:
            pattern = number_to_pattern(i, k)
            frequent_patterns.append(pattern)
    frequent_patterns.reverse()
    #return len(frequent_patterns)
    return frequent_patterns     


def frequent_words_with_mismatches_and_reverse_complements(genome, k, d):

    frequent_patterns = []
    close = {}
    frequency_array = {}

    for i in range(4 ** k - 1):
        close[i] = 0
        frequency_array[i] = 0

    for i in range(len(genome) - k + 1):
        neighborhood = neighbors(genome[i:i + k], d)
        for pattern in neighborhood:
            index = pattern_to_number(pattern)
            close[index] = 1
    for i in range(4 ** k - 1):
        if close[i] == 1:
            pattern = number_to_pattern(i, k)
            frequency_array[i] = approximate_pattern_count(genome, pattern, d) \
                + approximate_pattern_count(genome, reverse_complement(pattern), d)
    max_count = max(frequency_array.values())
    #print(max_count)

    for i in range(4 ** k - 1):
        if frequency_array[i] == max_count:
            pattern = number_to_pattern(i, k)
            frequent_patterns.append(pattern)
            if len(genome) == k:
                frequent_patterns.append(reverse_complement(pattern))
    return frequent_patterns

def neighbors(pattern, d):
    alt_bases = {'A': 'CGT', 'C': 'AGT', 'G': 'ACT', 'T': 'ACG'}
    suffix = pattern[1:]
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return {"A", "C", "G", "T"} 
    neighborhood = []
    suffix_neighbors = neighbors(suffix, d)
    for text in suffix_neighbors:
        if hamming_distance(suffix, text) < d:
            for nucleotide_x in alt_bases:
                neighborhood.append(nucleotide_x + text)
        else:
            neighborhood.append(pattern[0] + text)
    #neighborhood = neighborhood.sort()
    #print (neighborhood.sort())
    return neighborhood

start = time.time() 

#temp = frequent_words_with_mismatches_and_reverse_complements("TATAATTATGTTATGTGAATGTGTGGAAGTAATAATTATGTTGTGTGTGAATGTGAAGAAAATAATGAATATGTAATGAAAATAATGTGAAAATTGGAAGTTGTATTATTATGAATGTATGAATATTATGTTGTATGAATATAATAATGTGTGTTGAATTGAATTGTATGTTGGTTGAATTATAATGAAGAATATTATTATGAATGAATAATTGTGAATGAA", 6, 2)
#temp = frequent_words_with_mismatches("TGAGCTCAACCCTCAACCAACCTGAGCTCTGAGTGAGCTCCCCTCCCTCTCCCCTCTCGGTTCCCTCTCTGAGAACCAACCGGTTTGAGCTCCCCTGGTTCCCTAACCTGAGCTCCCCTCCCTCTCAACCCTCCTCGGTTGGTTGGTTCCCTCCCTCCCTCTCCCCTGGTTCCCTCTCTGAGGGTTGGTTGGTTAACCCCCTAACCAACCAACCTGAGCCCTCCCTCTCGGTTTGAGCTCCCCTGGTTGGTTGGTTCTCCTCCTCCTCCCCTAACCCTCAACCCCCTCCCTGGTTGGTTCCCTTGAGTGAGTGAGGGTTCTCAACCGGTTCTC", 7, 3)
#print (" ".join(str(x) for x in temp)) 
#print (approximate_pattern_count("ATA","ATA", 1))
#print(approximate_pattern_matching("ATTCTGGA", "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT", "3"))

#print(alist)
#print(hamming_distance("CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG", "ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT"))
#print(minimum_skew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"))
#print(minimum_skew(sys.stdin.readline()))
#print(approximate_pattern_count("CATGCCATTCGCATTGTCCCAGTGA", "CCC", 2))
#print(frequent_words_with_mismatches_and_reverse_complements("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1))
print(get_window(test_genome, 3764856, 1000))
#print(test_genome)
genome = get_window(test_genome, 3764856, 1000)
#print(genome)
print(frequent_words_with_mismatches_and_reverse_complements(genome, 9, 2))
#print(approximate_pattern_matching("CGGATCATC", genome, 2))
#TTATCCACA and TGTGGATAA for 9-mers
#temp = neighbors("TACCAAAAGATT", 2)
#print ("\n".join(str(x) for x in temp)) 
"""data = []
for line in sys.stdin:
    data.append(line)"""

print ('It took ', time.time()-start, 'seconds.')
