#!/bin/bash


import PatternToNumber
import NumberToPattern
import sys


def ComputingFreq(text, k):
#    print text;
    FrequenceArray = []
    for i in range(4**k):
        FrequenceArray.append(0)
    for i in range(len(text) - k):
#        print len(text),i;
        pattern = text[i:i+k]
#        print "Pattern: ", pattern;
        j = PatternToNumber.PatternToNumber(pattern)
#        print "FrequenceArray[" , j ,"] = ",FrequenceArray[j];
        FrequenceArray[j] = FrequenceArray[j] + 1
    return FrequenceArray
#        print "Count =",FrequenceArray[j];


def FasterFrequentWords(text, k):
    FrequentPattern = {}
    FreqArray = ComputingFreq(text, k)
    maxCount = max(FreqArray)
    print maxCount
    for i in range(4**k - 1):
        if FreqArray[i] == maxCount:
            pattern = NumberToPattern.NumberToPattern(i, k)
            FrequentPattern[i] = pattern
    return FrequentPattern


def FindFrequentWordsBySorting(text, k):
    FrequentPattern = {}
    Index = []
    Count = []  #A marker for counting
    for i in range(len(text)-k):
        pattern = text[i:i+k]
        Index.append(PatternToNumber.PatternToNumber(pattern))
        Count.append(1)
    Index.sort()
    for i in range(1, len(text)-k):
        if Index[i] == Index[i-1]:
            Count[i] = Count[i-1] + 1
    maxCount = max(Count)
    print maxCount
    for i in range(len(text)-k):
        if Count[i] == maxCount:
            pattern = NumberToPattern.NumberToPattern(Index[i], k)
            FrequentPattern[i] = pattern
    return FrequentPattern




if __name__ == "__main__" :
    file = open(sys.argv[1],"r");
    seq = file.readline();
#    print seq;
    k = file.readline();
#    print k;
    k = int(k);
    res1 = FasterFrequentWords(seq,k);
    print res1;
    res2 = FindFrequentWordsBySorting(seq,k);
    print res2;

#    for i in res:
#        print i,;


