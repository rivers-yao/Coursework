
def MinimumSkew(Genome):
    """
        Minimum Skew Problem: Find a position in a genome minimizing the skew.
        Input: A DNA string Genome.
        Output: All integer(s) i minimizing Skewi (Genome) among all values of
                i (from 0 to |Genome|).
    """
    SkewCodes = {'A' : 0, 'T' : 0, 'G' : 1, 'C' : -1}
    SkewValueList = []
    CurrentSkewValue = 0
    MinSkewPosition = []
    for base in Genome:
        SkewValueList.append(CurrentSkewValue)
        CurrentSkewValue += SkewCodes[base]
    print SkewValueList



seq = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
MinimumSkew(seq)




