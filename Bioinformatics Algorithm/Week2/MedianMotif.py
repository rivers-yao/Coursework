def NumberToPattern(index, k):
    table = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    if k == 1:
        return table[index]
    prefixIndex = index / 4
    r = index % 4
    prefixPattern = NumberToPattern(prefixIndex, k-1)
    symbol = table[r]
    return prefixPattern + symbol


def CalHammingDistance(seq1, seq2):
    assert len(seq1) == len(seq2)
    dist = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist += 1
    return dist


def DistanceBetweenPatternAndStrings(pattern, dna):
    res = 0
    for seq in dna:
        distance = len(pattern)
        seq = seq.upper()
        # print seq
        for j in range(len(seq)-len(pattern)+1):
            text = seq[j:j+len(pattern)]
            tempDist = CalHammingDistance(pattern, text)
            if tempDist <= distance:
                distance = tempDist
            # print distance
        res += distance
    return res

def MedianMotif(dna, k):
    distance = 0
    for i in dna:
        distance += len(i)
    print distance
    for i in range(4**k - 1):
        pattern = NumberToPattern(i, k)
        tempDist = DistanceBetweenPatternAndStrings(pattern, dna)
        if tempDist <= distance:
            distance = tempDist
            motif = pattern
    return motif








    # pattern = f.readline()
    # print pattern
    # kmer = int(f.readline())
    # mat = [map(float, line.strip("\n").split(" ")) for line in f]
    # print mat
#
# print MedianMotif(dnaSet, k)
# FindPortableKmer(pattern, kmer, mat)
# GenerateProfileMat(dnaSet)

# strings = ["GCG","AAG", "AAG", "ACG", "CAA"]
# print score(strings)
