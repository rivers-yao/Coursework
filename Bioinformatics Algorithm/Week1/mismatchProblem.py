import sys
import itertools


def reverse_complement(pattern):
    """Return reverse complement of pattern"""
    dna_revcomps = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    rev_comp = ""
    for base in pattern[::-1]:
        rev_comp += dna_revcomps[base]
    return rev_comp


def SymbolToNumber(symbol):
    num = -9
    if symbol == "A":
        num = 0
    elif symbol == "C":
        num = 1
    elif symbol == "G":
        num = 2
    elif symbol == "T":
        num = 3
    return num


def PatternToNumber(seq):
    if seq == "":
        return 0
    lastSymbol = seq[-1]
    seq = seq[:-1]
    return 4 * PatternToNumber(seq) + SymbolToNumber(lastSymbol)


def NumberToSymbol(index):
    if index == 0:
        return "A"
    elif index == 1:
        return "C"
    elif index == 2:
        return "G"
    elif index == 3:
        return "T"


def NumberToPattern(index, k):
    if k == 1:
        return NumberToSymbol(index)
    prefixIndex = index / 4
    r = index % 4
    prefixPattern = NumberToPattern(prefixIndex, k-1)
    symbol = NumberToSymbol(r)
    return prefixPattern + symbol


def hammingDistance(seq1, seq2):
    if len(seq1) != len(seq2):
        print "Not equal length of two strings"
        sys.exit()
    dist = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist += 1
    return dist


def ImmediateNeighors(pattern):
    res = []
    base = ['A', 'T', 'C', 'G']
    for i in range(0, len(pattern)):
        seq = list(pattern)
        symbol = pattern[i]
        for j in base:
            if symbol == j:
                continue
            seq[i] = j
            res.append("".join(seq))
    return res


def findNeighors(kmer, hammingDist, letters='AGCT'):
    """
       Helper generator to find all kmers within hamming_dist of the starting
       kmer.
    """
    for indices in itertools.combinations(range(len(kmer)), hammingDist):
        """
           This loop is to extract mismatch count base in kmer, for example:
           a sequence 'ACG' and find two mismatch: it will extract (0,1),(0,2),
           (1,2)
        """

        this_kmer = [[char] for char in kmer]
        # print this_kmer
        for index in indices:
            # print indices, index
            orgi_char = kmer[index]
            # print "orgi_char =", orgi_char
            this_kmer[index] = [l for l in letters if l != orgi_char]
            # print "this_kmer = ", this_kmer
        for poss in itertools.product(*this_kmer):
            # *syntax unpacks this_kmer into a list of individual letters
            #  and one list of all four possibilities
            yield ''.join(poss)


def neighbors(pattern, d):
    """
        Helper function to recursively find all kmers within d of the
        starting kmer.Recursive alternative to find_neighbors.
    """
    if d == 0:
        return set([pattern])
    if len(pattern) == 1:
        return set(['A', 'C', 'G', 'T'])
    neighborhood = set([])
    first, suffix = pattern[0], pattern[1::]
    suffix_neighbors = neighbors(suffix, d)
    for text in suffix_neighbors:
        if hammingDistance(text, suffix) < d:
            for base in 'ACGT':
                neighborhood = neighborhood.union([base + text])
                # print [base + text], neighborhood
        else:
            neighborhood = neighborhood.union([first + text])
    return neighborhood


def mismacthPatternFind(text, pattern, hammingDist):
    """
        Return the position of occurrances of pattern within text, including
        variants of pattern where the hamming distance from the variant to
        the pattern is less than hamming_dist
    """
    # print "pattern =", pattern
    index = []
    for i in range(len(text)-len(pattern)+1):
        if text[i:i+len(pattern)] == pattern:
            index.append(i)
        elif text[i:i+len(pattern)] == reverse_complement(pattern):
            index.append(i)
        elif hammingDistance(pattern, text[i:i+len(pattern)]) <= hammingDist:
            index.append(i)
    return index


def freqWordsWithMismatch(text, k, d):
    freqPattern = []
    mismatchSeqIndex = list()
    frequentArray = list()
    for i in range(4**k):
        mismatchSeqIndex.append(0)
        frequentArray.append(0)
    for i in range(len(text) - k):
        neighborhood = neighbors(text[i:i+k], d)
        # print neighborhood
        # Trans all neighbors to index
        for p in neighborhood:
            index = PatternToNumber(p)
            # if mismatchSeqIndex[index] == 1:
            #    print mismatchSeqIndex.count(1)
            mismatchSeqIndex[index] = 1
    for i in range(4**k):
        # Search all the genome by slide on the genome

        if mismatchSeqIndex[i] == 1:
            pattern = NumberToPattern(i, k)
            #print "Pattern =", pattern
            frequentArray[i] = len(mismacthPatternFind(text, pattern, d))
    # print frequentArray,len(frequentArray)
    maxCount = max(frequentArray)
    for i in range(4**k - 1):
        if frequentArray[i] == maxCount:
            pattern = NumberToPattern(i, k)
            freqPattern.append(pattern)
    return (freqPattern, maxCount)
    #print freqPattern







if __name__ == "__main__":
    # input = open("/Users/Heyao/Courses/bioinformatics/dataset_9_7.txt", "r")
    # text = input.readline()
    # text = text.strip()
    # pattern = input.readline()
    # pattern = pattern.strip()
    # d = input.readline()
    # d = int(d)
    # print len(pattern)
    # print text
    text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
    res = freqWordsWithMismatch(text, 4, 1)
    print res


    #f1 = open("test.txt","w+")
    #for i in res:
     #   print >> f1, i,

