#!/bin/bash

def SymbolToNumber(symbol):
    num = -9;
    if symbol == "A":
        num = 0
    elif symbol == "C":
        num = 1;
    elif symbol == "G":
        num = 2;
    elif symbol == "T":
        num = 3;
    return num;

def PatternToNumber(seq):
    if seq == "" :
        return 0;
    lastSymbol = seq[-1];
    seq = seq[:-1];
#print "seq = ", seq;
#    print "lastSymbol = ",lastSymbol;
    return 4 * PatternToNumber(seq) + SymbolToNumber(lastSymbol);


if __name__ == "__main__" :
    seq = "ATGCAA";
    num = PatternToNumber(seq);
    print "Result:",num;
# print test;
