#!/bin/python 
import sys
def PatternCount(seq,k):
    res = {};
    count = 0;
    k = int(k);
#for is to split kmer
    for i in range(0,(len(seq) - k)):
        pattern = seq[i:k+i];        
        if pattern in res:
            res[pattern] += 1;
        else:
            res[pattern] = 1;
    return res; 

if __name__ == "__main__":
    lineCount = 0;
    k =  sys.argv[2];
    with open(sys.argv[1],"r") as f:
        for line in f:
            lineCount = lineCount + 1;
            if lineCount % 2 != 0 : 
                continue;
            else : 
                tt = PatternCount(line , k);
    maxFreq = max(tt.values());
    
    res = [k for k,v in tt.items() if v == maxFreq];
    print res;
