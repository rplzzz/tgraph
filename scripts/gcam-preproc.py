#!/usr/bin/python

from sys import argv
import re


## Some patterns we will need
header = re.compile(r'digraph\s+(\w+)\s*{')
footer = re.compile(r'}')
labelx = re.compile(r'(\w+)\[label\s*=\s*"([\w\s\-]+)"[^\]]*\]')
edge   = re.compile(r'(\w+)\s*->\s*(\w+)')
# This one is to replace whitespace with underscores
fixws  = re.compile(r'[\s\-]+')

## translation table
transtbl = {}
## name table
nametbl = {}
## edge list
edges = []
## graph name
gname = "graph"


infile = open(argv[1],"r")
outfile = open(argv[2],"w")

for line in infile.readlines():
    m = header.search(line)
    if m:
        gname = m.group(1)
        continue
    
    m = footer.search(line)
    if m:
        break

    m = labelx.search(line)
    if m:
        idnum = m.group(1)
        label = m.group(2)
        labelfix = fixws.sub("_",label)
        i = 1
        # ensure a unique label
        while labelfix in nametbl:
            i += 1
            labelfix = "%s%d" % (labelfix,i)
        transtbl[idnum] = labelfix
        nametbl[labelfix] = 1
        #print "%s label= %s\n" % (idnum,labelfix)
        
        continue

    m = edge.search(line)
    if m:
        e1 = m.group(1)
        e2 = m.group(2)
        edges.append((e1,e2))
        continue

    print "Discarded: ", line

## finished all the lines, now output the graph
outfile.write("digraph %s {\n" % gname)
for edge in edges:
    outfile.write("%s -> %s;\n" % (transtbl[edge[0]],transtbl[edge[1]]))

outfile.write("}\n")

