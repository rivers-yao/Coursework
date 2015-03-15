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

def GenerateProfileMatWithPesudoCount(dna):
    seqNum = float(len(dna))
    nucleotide = ['A', 'C', 'G', 'T']
    matrix = []
    for i in range(len(dna[0])):
        base_i = [seq[i] for seq in dna]
        # print base_i
        # print "seqNum = ", seqNum + 1
        colProfile = [float(base_i.count(n) + 1)/seqNum + 4 for n in nucleotide]
        matrix.append(colProfile)
    return [list(i) for i in zip(*matrix)]


def GenerateProfileDict(matrix, j):
    return {'A': matrix[0][j], 'C': matrix[1][j], 'G': matrix[2][j], 'T': matrix[3][j]}


# def ScoreMotifs(motifs):
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



def GreedyMotifSearch(dna, k, t, PesudoCount = False):
    BestMotifs = [seq[0:k] for seq in dna]
    # print BestMotifs
    # profileMat = GenerateProfileMat(BestMotifs)
    seq = dna[0]
    # print seq
    for i in range(len(seq) - k):
        kmers = seq[i:i+k]
        # print "Kmer = ", kmers
        motifs = [kmers]
        for j in range(1, len(dna)):
            if PesudoCount == True:
                mat = GenerateProfileMatWithPesudoCount(motifs)
            else:
                mat = GenerateProfileMat(motifs)
            # print "Dna =", dna[j]
            # print mat
            # print "Find Portable Kmer:"
            tempMotif = FindPortableKmer(dna[j], k, mat)
            motifs.append(tempMotif)
        if score(motifs) < score(BestMotifs):
            BestMotifs = motifs
    return BestMotifs


with open("/Users/Heyao/Courses/Coursework/Bioinformatics Algorithm/Week2/dataset_160_9.txt", "r") as f:
     k = f.readline().strip()
     dnaSet = [line.strip() for line in f]
     print dnaSet
res = GreedyMotifSearch(dnaSet, 12, 25, True)
for i in res:
    print i
