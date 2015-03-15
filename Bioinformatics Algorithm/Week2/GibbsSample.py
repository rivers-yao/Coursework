from random import randint

def CalHammingDistance(seq1, seq2):
    assert len(seq1) == len(seq2)
    dist = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist += 1
    return dist

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


def GibbsSampler(dna, k, t, N):
    randomIndex = [randint(0, len(dna[0])-k) for i in range(len(dna))]
    Motifs = [seq[i:i+k] for seq, i in zip(dna, randomIndex)]
    # print randomIndex
    # print Motifs
    BestMotifs = Motifs
    for j in range(N):
        i = randint(0, len(dna) - 1)
        # print "i = ", i
        SampleMotifs = [dna[index] for index in range(len(dna)) if index != i]
        # print SampleMotifs
        ProfileMatrix = GenerateProfileMatWithPesudoCount(Motifs)
        NewMotif = FindPortableKmer(dna[i], k, ProfileMatrix)
        # print Motifs
        Motifs[i] = NewMotif
        # print Motifs
        if score(Motifs) < score(BestMotifs):
            BestMotifs = Motifs
        # print score(BestMotifs)
    # print score(BestMotifs)
    return BestMotifs
        # print ProfileMatrix

dnaSet = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", \
          "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", \
          "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]


# output = ["TCTCGGGG", "CCAAGGTG", "TACAGGCG", "TTCAGGTG", "TCCACGTG"]
# print score(output)


# res = []
# for i in range(200):
    # res.append(GibbsSampler(dnaSet, 8, 5, 100))

res = []
with open("/Users/Heyao/Courses/Coursework/Bioinformatics Algorithm/Week2/dataset_163_4.txt", "r") as f:
     k = f.readline().strip()
     dnaSet = [line.strip() for line in f]
     print dnaSet
for i in range(0,20):
    res.append(GibbsSampler(dnaSet, 15, 20, 2000))
    # print "***************"
s = []
for i in res:
    s.append(score(i))
print min(s)
print s.index(min(s))

for i in res[s.index(min(s))]:
    print i

# res = GibbsSampler(dnaSet, 8, 5, 100)

# for i in res:
    # print i
