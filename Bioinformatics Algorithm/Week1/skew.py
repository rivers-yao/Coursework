import sys


def skewValue(seq):
    skew = []
    for i in range(len(seq)+1):
        skew.append(0)
    for i in range(0, len(seq)):
        if seq[i] == "C":
            skew[i+1] = skew[i] - 1
        elif seq[i] == "G":
            skew[i+1] = skew[i] + 1
        else:
            skew[i+1] = skew[i]
    return skew




if __name__ == "__main__":
    f = open(sys.argv[1], "r")
    input = f.readline()

    res = skewValue(input)
    minValue = min(res)
    print minValue
    minSkew = [i for i in range(len(res)) if res[i] == minValue]
    print minSkew,
