"""PatternCount Algorithm for DNA Sequence"""

import sys

def PatternToNumber(text):
    table = {'A':0, 'C':1, 'G':2, 'T':3}
    if len(text) == 0:
        return 0
    symbol = text[-1:]
    text = text[:-1]
    return 4*PatternToNumber(text) + table[symbol]


def NumberToPattern(index, k):
    table = {0:'A', 1:'C', 2:'G', 3:'T'}
    if k == 1:
        return table[index]
    prefixIndex = index / 4
    r = index % 4
    prefixPattern = NumberToPattern(prefixIndex, k-1)
    symbol = table[r]
    return prefixPattern + symbol


def PatternCount(sequence, kmer):
    """
       Input Sequence and kmer(the pattern you want to check),
       return the number of times of the kmer appears in the sequence
    """
    count = 0
    for i in range(len(sequence) - len(kmer) ):
        if sequence[i:i+len(kmer)] == kmer:
            count += 1
    return count


def ReverseComplement(sequence):
    table = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    ReverseSeq = ""
    for base in sequence[::-1]:
        ReverseSeq += table[base]
    return ReverseSeq


def FindFrequentPattern(sequence, k):
    """
       A straightforward algorithm for finding the most frequent words in
       a string Text checks all k-mers appearing in this string and then
       computes how many times each k-mer appears in Text
       ---------------------------------------------------------------------
       Notice: I(O) == len(sequence)^2 as PatternCount() need to loop sequence
       again
    """
    PatternArray = set()
    CountArray = []
    for i in range(len(sequence) - k):
        CountArray.append(0)
    for i in range(len(sequence) - k):
        kmer = sequence[i:i+k]
        CountArray[i] = PatternCount(sequence, kmer)
    maxCount = max(CountArray)
    for i in range(len(sequence) - k):
        if CountArray[i] == maxCount:
            PatternArray.add(sequence[i:i+k])
    return PatternArray

def ComputingFrequenceArray(sequence, k):
    """
        Use frequency array to compute most frequent k-mers in sequence
        ---------------------------------------------------------------
        Notice:I(O) == 4**k
    """
    FrequencyArray = []
    for i in range(4**k - 1):
        FrequencyArray.append(0)
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        index = PatternToNumber(kmer) #This is for hashing the kmer sequence
        FrequencyArray[index] = FrequencyArray[index] + 1
    return FrequencyArray


def FindFrequentPatternByHashing(sequence, k):
    """
       A Hashing algorithm for finding the most frequent words in
       a string Text checks all k-mers appearing in this string and then
       computes how many times each k-mer appears in Text
       ---------------------------------------------------------------------
       Notice: I(O) == 4**k
    """
    FrequenceArray = ComputingFrequenceArray(sequence, k)
    maxCount = max(FrequencyArray)
    PatternArray = set()
    for i in range(4**k - 1):
        if FrequenceArray[i] == maxCount:
            pattern = NumberToPattern(i, k)
            PatternArray.add(pattern)
    return (PatternArray, maxCount)


def FindFrequentPatternBySorting(sequence, k):
    """
        This algorithm do not initialize all kmer hashing array.
        Instead , it generates a hashing arry only when a kmer exists
        Then sorted this array and calculate the max value in this array
        so it reduce the I(O)
    """
    PatternArray = []
    IndexArray = []
    CountArray = []
    for i in range(len(sequence) - k + 1):
        pattern = sequence[i:i+k]
        IndexArray.append(PatternToNumber(pattern))
        CountArray.append(1)
    SortedIndexArray = sorted(IndexArray)
    for i in range(1, len(sequence) - k + 1):
        if SortedIndexArray[i] == SortedIndexArray[i-1]:
            CountArray[i] = CountArray[i-1]
    maxCount = max(CountArray)
    for i in range(len(sequence) - k + 1):
        if CountArray[i] == maxCount:
            pattern = NumberToPattern(SortedIndexArray[i], k)
            PatternArray.append(pattern)
    return PatternArray


def ClumpFinding(Genome, k, t, L):
    """
        Find patterns forming clumps in a string.
        Input: A string Genome, and integers k, L, and t.
        Output: All distinct k-mers forming (L, t)-clumps in Genome.
        ---------------------------------------------------------------
        After computing the frequency array for the current window,
        it identifies (L, t)-clumps simply by finding which k-mers occur
        at least t times within the window.
    """
    clumps = []
    FrequencePattern = set()
    for i in range(4**k - 1):
        clumps.append(0)
    for i in range(len(Genome) - L):
        ClumpText = Genome[i:i+L]
        FrequenceArray = ComputingFrequenceArray(ClumpText, k)
        for j in range(4**k - 1):
            if FrequenceArray[j] >= t:
                clumps[j] = 1
    for i in range(4**k - 1):
        if clumps[i] == 1:
            pattern = NumberToPattern(i, k)
            FrequencePattern.add(pattern)
    return FrequencePattern


def HammingDistance(seq1, seq2):
    """Input two strings of equal length and count mismatches between them"""
    assert len(seq1) == len(seq2)
    dist = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist += 1
    return dist


def ApproximatePatternFind(seq, pattern, HammingDist):
    PatternOccurrence = []
    for i in range(len(seq)-len(pattern) + 1):
        text = seq[i:i+len(pattern)]
        if text == pattern:
            PatternOccurrence.append(i)
        elif HammingDistance(text, pattern) <= HammingDist:
            PatternOccurrence.append(i)
    return PatternOccurrence


def ApproximatePatternCount(seq, pattern, HammingDist):
    count = 0
    # print len(seq), len(pattern)
    # print len(seq) - len(pattern) + 1
    for i in range(len(seq)-len(pattern) + 1): # Remember +1 important
        text = seq[i:i+len(pattern)]
        if HammingDistance(text, pattern) <= HammingDist:
            count = count + 1
    return count

def GenerateNeighborhood(pattern, d):
    """
        Try to add each nucleotide in one loop of sequence length
        For example: ACG
        Firstly, test HammingDistance(G,set([A,C,G,T])), then the
        HammingDistance means that whether G == A or not
        if G != A, the HammingDistance == 1, we can add the orginal base
        before G, such as 'CA'
        Recursive this process, then we get all Neighborhood of the sequence
    """
    if d == 0:
        return set([pattern])
    if len(pattern) == 1:
        return set(['A', 'C', 'G', 'T'])
    Neighborhood = set([])
    first, suffix = pattern[0], pattern[1:]
    # print "first =", first
    # print "suffix =", suffix
    NeighborhoodSuffix = GenerateNeighborhood(suffix, d)
    for text in NeighborhoodSuffix:
        # print text
        # print HammingDistance(text, suffix)
        if HammingDistance(text, suffix) < d:
            for base in 'ACGT':
                Neighborhood = Neighborhood.union([base + text])
                # print Neighborhood
        else:
            Neighborhood = Neighborhood.union([first + text])
            # print Neighborhood
    return Neighborhood

def FindFrequentPatternWithMismatch(sequence, k, d, countRC = False):
    """
        Frequent Words with Mismatches Problem: Find the most frequent k-mers
        with mismatches in a string.
        Input: A string Text as well as integers k and d.
        Output: All most frequent k-mers with up to d mismatches in Text.
    """
    FrequenceArray = []
    Close = []
    for i in range(4**k ):
        FrequenceArray.append(0)
        Close.append(0)
    for i in range(len(sequence)-k):
        pattern = sequence[i:i+k]
        NeighborPattern = GenerateNeighborhood(pattern, d)

        for string in NeighborPattern:
            # print pattern, string, PatternToNumber(string), max(FrequenceArray)
            index = PatternToNumber(string)
            Close[index] = 1
            # FrequenceArray[index] = FrequenceArray[index] + 1
    # print "Loop 2 done"
    for i in range(4**k):
        if Close[i] == 1:
            pattern = NumberToPattern(i, k)
            FrequenceArray[i] = ApproximatePatternCount(sequence, pattern, d)
            if countRC == True:
                patternRC = ReverseComplement(pattern)
                # print FrequenceArray[i]
                FrequenceArray[i] = FrequenceArray[i] + ApproximatePatternCount(sequence, patternRC, d)

            # print FrequenceArray[i]
    maxCount = max(FrequenceArray)
    print maxCount
    FrequencePattern = set([])
    for i in range(4**k):
        if FrequenceArray[i] == maxCount:
            pattern = NumberToPattern(i, k)
            FrequencePattern.add(pattern)
    return FrequencePattern

def FindFrequentPatternWithMismatchBySorting(sequence, k, d, countRC=False):
    FrequentPattern = set([])
    NeighborPatternArray = []
    IndexArray = []
    CountArray = []
    for i in range(len(sequence)-k+1):
        pattern = sequence[i:i+k]
        NeighborPatternArray.extend([string for string in GenerateNeighborhood(pattern, d)])
        if countRC == True:
            patternRC = ReverseComplement(pattern)
            NeighborPatternArray.extend([string for string in GenerateNeighborhood(patternRC, d)])
    print "Loop 1 done"
    for i in range(len(NeighborPatternArray)):
        pattern = NeighborPatternArray[i]
        IndexArray.append(PatternToNumber(pattern))
        CountArray.append(1)
    SortedIndexArray = sorted(IndexArray)
    for i in range(len(NeighborPatternArray)-1):
        if SortedIndexArray[i] == SortedIndexArray[i+1]:
            CountArray[i+1] = CountArray[i] + 1
    maxCount = max(CountArray)
    print maxCount
    for i in range(len(NeighborPatternArray)-1):
        if CountArray[i] == maxCount:
            FrequentPattern.add(NumberToPattern(SortedIndexArray[i], k))
    return FrequentPattern





# test = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
# print HammingDistance(s1, s2)
# print FindFrequentPattern(test, 4)
# print PatternToNumber("AGT")
# print NumberToPattern(45, 4)
# f = open("/Users/Heyao/Courses/Coursework/Bioinformatics Algorithm/Week1/dataset_9_8.txt", "r")
# pattern = f.readline()
# pattern = pattern.strip()
# print pattern
# print pattern
# pattern = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
pattern = "CTTGCCGGCGCCGATTATACGATCGCGGCCGCTTGCCTTCTTTATAATGCATCGGCGCCGCGATCTTGCTATATACGTACGCTTCGCTTGCATCTTGCGCGCATTACGTACTTATCGATTACTTATCTTCGATGCCGGCCGGCATATGCCGCTTTAGCATCGATCGATCGTACTTTACGCGTATAGCCGCTTCGCTTGCCGTACGCGATGCTAGCATATGCTAGCGCTAATTACTTAT"
# print ApproximatePatternCount(pattern, "GCACACAGAC", 2)
print pattern
# print FindFrequentPatternWithMismatch(pattern, 7, 2, True)
print FindFrequentPatternWithMismatchBySorting(pattern, 9, 3, True)
