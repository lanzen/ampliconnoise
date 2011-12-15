#!/usr/bin/python
import sys

fastafile = open(sys.argv[1], 'r')
cropPos=int(sys.argv[2])

for line in fastafile:
    if line[0]=='>':
        print line[:-1]
        fNext=True
    elif fNext:
        print(line[cropPos:-1])
        fNext=False
    else: print(line[:-1])

fastafile.close()
