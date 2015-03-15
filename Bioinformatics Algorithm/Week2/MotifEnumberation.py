

def HammingDistance(seq1, seq2):
    """Input two strings of equal length and count mismatches between them"""
    assert len(seq1) == len(seq2)
    dist = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist += 1
    return dist


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


def ApproximatePatternFindInAllSeq(seq, pattern, HammingDist):
    res = []
    for subseq in seq:
        inOneSeq = []
        for j in range(len(subseq)-len(pattern)+1):
            text = subseq[j:j+len(pattern)]
            if HammingDistance(text, pattern) <= HammingDist:
                inOneSeq.append(True)
                break
            else:
                inOneSeq.append(False)
        # print inOneSeq
        if any(inOneSeq):
            res.append(True)
        else:
            res.append(False)
        # print "res = ", res
    if all(res):
        return True
    else:
        return False


def GetKmer(seq, k):
    Kmer = set([])
    for i in range(len(seq)-k+1):
        pattern = seq[i:i+k]
        Kmer.add(pattern)
    return Kmer


def MotifEnumeration(dna, k, d):
    """
        Brute Force Algorithm
        It generate all kmer and all kmer Neighborhood and check whether exists
        in all DNA seq you gave
    """
    pattern = set([])
    kmer = set([])
    for seq in dna:
        kmer = kmer.union(GetKmer(seq, k))
    for mer in kmer:
        # print "mer = ", mer
        for m in GenerateNeighborhood(mer, k):
            # print m
            if ApproximatePatternFindInAllSeq(dna, m, d):
                # print "Here"
                pattern.add(m)
    return pattern


def GenerateProfileMatrix(dna):
    l = len(dna[0])
    colMatrix = []
    base = ['A', 'C', 'G', 'T']
    for i in range(l):
        colMatrix.append([seq[i] for seq in dna])
    for col in colMatrix:
    print [col.count(b) / len(dna) for b in base]


def GreedyMotifSearch(dna, k, t):
    firstkmer = [seq[0:k] for seq in dna]
    print firstkmer
    BestMotifs = GenerateProfileMatrix(firstkmer)

dnaSet = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]
GenerateProfileMatrix(dnaSet)
# GreedyMotifSearch(dnaSet, 3, 5)
