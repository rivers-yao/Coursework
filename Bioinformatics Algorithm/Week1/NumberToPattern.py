#!/bin/bash/

def NumberToSymbol(index):
    if index == 0:
        return "A";
    elif index == 1:
        return "C";
    elif index == 2:
        return "G";
    elif index == 3:
        return "T";

def NumberToPattern(index,k):
    if k == 1:
        return NumberToSymbol(index);# Here is the end 
    prefixIndex = index / 4;
#    print "prefixIndex : ",prefixIndex;
    r = index % 4;
    prefixPattern = NumberToPattern(prefixIndex,k-1);
#    print prefixPattern;
    symbol = NumberToSymbol(r);
    return prefixPattern + symbol;


if __name__ == "__main__":
    res =  NumberToPattern(5437,8);
    print res;

