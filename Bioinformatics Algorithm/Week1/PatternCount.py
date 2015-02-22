
import sys


def PatternCountMatch(seq, pattern):
    count = 0
    position = []
#    for is to split kmer
#    seq = seq.strip("\n");
    pattern = pattern.strip("\n")
#    print pattern;
    seqLen = len(seq)
    patternLen = len(pattern)
    print seqLen, patternLen
    for i in range(0, (seqLen - patternLen)):
        substring = seq[i:i+patternLen]
# print substring;
        if substring == pattern:
            count += 1
            position.append(i)
    print count

    for i in position:
        print i,

if __name__ == "__main__":
    file = open(sys.argv[1], "r")
    pattern = file.readline()
    seq = file.readline()
#   print seq;
    PatternCountMatch(seq, pattern)
