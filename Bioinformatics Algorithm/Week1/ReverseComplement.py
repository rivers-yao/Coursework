#!/bin/bash

import sys


file = open(sys.argv[1], "r")
rev_seq = ""
for seq in file:
    for i in seq:
        print i
        if i == "A":
            rev_seq += "T"
        elif i == "T":
            rev_seq += "A"
        elif i == "C":
            rev_seq += "G"
        elif i == "G":
            rev_seq += "C"
res = open("result.txt", "w")
res.write(rev_seq[::-1])
