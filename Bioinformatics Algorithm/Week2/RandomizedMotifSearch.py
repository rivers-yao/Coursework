from random import randint


def CalHammingDistance(seq1, seq2):
    assert len(seq1) == len(seq2)
    dist = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist += 1
    return dist


def GenerateProfileMat(dna):
    seqNum = float(len(dna))
    nucleotide = ['A', 'C', 'G', 'T']
    matrix = []
    for i in range(len(dna[0])):
        base_i = [seq[i] for seq in dna]
        # print base_i
        colProfile = [float(base_i.count(n))/seqNum for n in nucleotide]
        matrix.append(colProfile)
    return [list(i) for i in zip(*matrix)]

def FindPortableKmer(text, k, matrix):
    # print text
    maxPr = 0
    prList = []
    for i in range(len(text)-k+1):
        KmerPr = 1
        pattern = text[i:i+k]
        # print pattern
        for j in range(len(pattern)):
            profile = GenerateProfileDict(matrix, j)
            # print profile
            KmerPr *= profile[pattern[j]]
        # print "KmerPr =", KmerPr
        prList.append(KmerPr)
    i = prList.index(max(prList))
    maxKmer = text[i:i+k]
    return maxKmer


def GenerateProfileDict(matrix, j):
    return {'A': matrix[0][j], 'C': matrix[1][j], 'G': matrix[2][j], 'T': matrix[3][j]}


def score(motifs):
    '''Returns the score of the dna list motifs.'''
    score = 0
    for i in range(len(motifs[0])):
        motif = ''.join([motifs[j][i] for j in range(len(motifs))])
        # print motif
        # print [homogeneous*len(motif) for homogeneous in 'ACGT']
        # Calculate the min score between motif and [AAAAA, CCCCC, GGGGG, TTTTT]
        # avoiding find the consensus strings
        score += min([CalHammingDistance(motif, homogeneous*len(motif)) for homogeneous in 'ACGT'])
    return score


def GenerateProfileMatWithPesudoCount(dna):
    seqNum = float(len(dna))
    nucleotide = ['A', 'C', 'G', 'T']
    matrix = []
    for i in range(len(dna[0])):
        base_i = [seq[i] for seq in dna]
        # print base_i
        # print "seqNum = ", seqNum + 4
        colProfile = [float(base_i.count(n) + 1)/float(seqNum + 4) for n in nucleotide]
        # print "colProfile = ", colProfile
        matrix.append(colProfile)
    return [list(i) for i in zip(*matrix)]


def RandomizedMotifSearch(dna, k):
    # Random select initial motifs
    kmerIndex = [randint(0, len(dna[0]) - k) for i in range(len(dna))]
    Motifs = [dna[i][j:j+k] for i, j in enumerate(kmerIndex)]
    print "First =" , Motifs
    BestMotifs = Motifs
    while True:
        ProfileMatrix = GenerateProfileMatWithPesudoCount(Motifs)
        # print ProfileMatrix
        Motifs = [FindPortableKmer(dna[i], k, ProfileMatrix) for i in range(len(dna))]
        # print "Motifs = ", Motifs
        # print score(Motifs), score(BestMotifs)
        if score(Motifs) < score(BestMotifs):
            BestMotifs = Motifs
        else:
            return BestMotifs

# dnaSet = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]
dnaSet = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", \
            "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", \
            "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", \
            "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
# a = ["TCTCGGGG", "CCAAGGTG", "TACAGGCG", "TTCAGGTG", "TCCACGTG"]
# print score(a)
with open("/Users/Heyao/Courses/Coursework/Bioinformatics Algorithm/Week2/dataset_161_5.txt", "r") as f:
     k = f.readline().strip()
     dnaSet = [line.strip() for line in f]
     print dnaSet
# print score(res)
# for i in range(0,10):
res = []
for i in range(0,100):
    res.append(RandomizedMotifSearch(dnaSet, 15))
s = []
for i in res:
    s.append(score(i))
print s.index(min(s))
for i in res[s.index(min(s))]:
    print i

