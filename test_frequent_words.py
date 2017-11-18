import sys

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


def quotient(index):
    return index // 4


def remainder(index):
    return index % 4


def pattern_to_number(pattern):
    if pattern == '':
        return 0

    symbol = pattern[len(pattern) - 1]
    prefix = pattern[:len(pattern) - 1]

    return 4 * pattern_to_number(prefix) + symbol_to_number(symbol)


def number_to_pattern(index, k):
    if k == 1:
        return number_to_symbol(index)

    prefix_index = quotient(index)
    r = remainder(index)

    symbol = number_to_symbol(r)

    prefix_pattern = number_to_pattern(prefix_index, k - 1)

    return prefix_pattern + symbol


def hamming_distance(p, q):
    dif = 0
    for i in range(len(p)):
        if not q[i] == p[i]:
            dif += 1
    return dif


def neighbors(pattern, d):
    alt_bases = {'A': 'CGT', 'C': 'AGT', 'G': 'ACT', 'T': 'ACG'}
    suffix = pattern[1:]

    if d == 0:
        return pattern

    if len(pattern) == 1:
        return {'A', 'C', 'G', 'T'}

    neighborhood = []
    suffix_neighbors = neighbors(suffix, d)
    for st in suffix_neighbors:
        if hamming_distance(suffix, st) < d:
            for nucleotide_x in alt_bases:
                neighborhood.append(nucleotide_x + st)
        else:
            neighborhood.append(pattern[0] + st)
    return neighborhood

def reverse_complement(pattern):
    complement = []

    for acid in pattern:
        if acid == 'A':
            acid = 'T'
        elif acid == 'T':
            acid = 'A'
        elif acid == 'G':
            acid = 'C'
        else:
            acid = 'G'

        complement.append(str(acid))

    complement.reverse()

    return ''.join(complement)

def computing_frequencies(text, k):
    frequency_array = {}
    keys = {}

    for i in range(len(text)-k+1):
            pattern = text[i:i+k]
            j = pattern_to_number(pattern)
            if keys.get(j):
                frequency_array[j] += 1
            else:
                frequency_array[j] = 1
                keys[j]=j
    return frequency_array

def faster_frequent_words(text, k):
        frequent_patterns = []
        frequency_array = computing_frequencies(text, k)
        max_count = max(frequency_array.values())

        for i in frequency_array:
                if frequency_array[i] == max_count:
                    pattern = number_to_pattern(i, k)
                    frequent_patterns.append(pattern)
        return frequent_patterns

def pattern_count(text, pattern):
    count = 0
    for i in range(len(text)-len(pattern)+1):
        if text[i:i+len(pattern)] == pattern:
            count += 1
    return count

def approximate_pattern_count(text, pattern, d):
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        pattern_prima = text[i:i + len(pattern)]
        if hamming_distance(pattern, pattern_prima) <= d:
            count += 1
    return count


def reverse_complement(pattern):
    complement = []

    for acid in pattern:
        if acid == 'A':
            acid = 'T'
        elif acid == 'T':
            acid = 'A'
        elif acid == 'G':
            acid = 'C'
        else:
            acid = 'G'

        complement.append(str(acid))

    complement.reverse()

    return ''.join(complement)


def frequent_words_with_mismatches(text, k, d):
    if d == 0:
        return faster_frequent_words(text, k)
    frequent_patterns = []
    close = {}
    frequency_array = {}
    for i in range(4 ** k - 1):
        close[i] = 0
        frequency_array[i] = 0
    for i in range(len(text) - k+1):
        neighborhood = neighbors(text[i:i + k], d)
        for pattern in neighborhood:
            index = pattern_to_number(pattern)
            close[index] = 1
    for i in range(4 ** k - 1):
        if close[i] == 1:
            pattern = number_to_pattern(i, k)
            frequency_array[i] = approximate_pattern_count(text, pattern, d)
    max_count = max(frequency_array.values())
    for i in range(4 ** k - 1):
        if frequency_array[i] == max_count:
            pattern = number_to_pattern(i, k)
            frequent_patterns.append(pattern)
    frequent_patterns.reverse()
    return frequent_patterns

def frequent_words_with_mismatches_and_reverse_complements(text, k, d):
    if d == 0 and len(text) != k:
        return faster_frequent_words(reverse_complement(text)+text, k)

    frequent_patterns = []
    close = {}
    frequency_array = {}

    for i in range(4 ** k - 1):
        close[i] = 0
        frequency_array[i] = 0

    for i in range(len(text) - k + 1):
        neighborhood = neighbors(text[i:i + k], d)
        for pattern in neighborhood:
            index = pattern_to_number(pattern)
            close[index] = 1
    for i in range(4 ** k - 1):
        if close[i] == 1:
            pattern = number_to_pattern(i, k)
            frequency_array[i] = approximate_pattern_count(text, pattern, d) \
                + approximate_pattern_count(text, reverse_complement(pattern), d)
    max_count = max(frequency_array.values())

    for i in range(4 ** k - 1):
        if frequency_array[i] == max_count:
            pattern = number_to_pattern(i, k)
            frequent_patterns.append(pattern)
            if len(text) == k:
                frequent_patterns.append(reverse_complement(pattern))
    return frequent_patterns

print(frequent_words_with_mismatches_and_reverse_complements("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1))