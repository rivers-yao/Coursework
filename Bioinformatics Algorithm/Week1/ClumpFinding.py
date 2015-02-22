import sys
import ComputingFrequency
import NumberToPattern
import PatternToNumber


def ClumpFinding(Genome, k, t, L):
    clump = []
    freqPattern = []
    for i in range(4**k - 1):
        clump.append(0)
    for i in range(len(Genome)-L):
        text = Genome[i:i+L]
        FrequencyArray = ComputingFrequency.ComputingFreq(text, k)
        for j in range(4 ** k - 1):
            if FrequencyArray[j] >= t:
                clump[j] = 1

    for i in range(4**k - 1):
        if clump[i] == 1:
            pattern = NumberToPattern.NumberToPattern(i, k)
            freqPattern.append(pattern)
    return freqPattern


def BetterClumpingFinding(Genome, k, t, L):
    freqPattern = []
    clump = []
    freqArray = []
    for i in range(4**k - 1):
        clump.append(0)
    text = Genome[0:L]
    freqArray = ComputingFrequency.ComputingFreq(text, k)
    for i in range(4**k - 1):
        if freqArray[i] >= t:
            clump[i] = 1
    for i in range(1, len(Genome) - L):
        FirstPattern = Genome[i-1:k]
        j = PatternToNumber.PatternToNumber(FirstPattern)
        freqArray[j] = freqArray[j] - 1
        LastPattern = Genome[i+L-k:k]
        j = PatternToNumber.PatternToNumber(LastPattern)
        freqArray[j] = freqArray[j] + 1
        if freqArray[j] >= t:
            clump[j] = 1
    for i in range(4**k - 1):
        if clump[i] == 1:
            pattern = NumberToPattern.NumberToPattern(i, k)
            freqPattern.append(pattern)
    return freqPattern


if __name__ == "__main__":
    input = open(sys.argv[1])
    seq = input.readline()
    res = BetterClumpingFinding(seq, 9, 1, 500)
    print res
