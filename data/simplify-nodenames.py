#!/usr/bin/env python

import re
import sys

def inttostring(N):
    ## Required: N>=0
    r = N%26
    string = chr(ord("A")+r)
    N /= 26
    while(N>0):
        r = N%26
        string += chr(ord("A") + r-1)
        N /= 26
    return string
    

## this processing will simply kill any unattached nodes

linematch = re.compile('([A-Za-z0-9_]+) -> ([A-Za-z0-9_]+)')
rename_tbl = {}
count = 0

for line in sys.stdin.readlines():
    m = linematch.search(line)
    if(not m):
        sys.stdout.write(line)
    else:
        n1 = m.group(1)
        n2 = m.group(2)

        if(not rename_tbl.has_key(n1)):
            rename_tbl[n1] = inttostring(count)
            count += 1
        if(not rename_tbl.has_key(n2)):
            rename_tbl[n2] = inttostring(count)
            count += 1

        sys.stdout.write("\t%s -> %s;\n" % (rename_tbl[n1], rename_tbl[n2]))

