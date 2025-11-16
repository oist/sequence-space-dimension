#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Function to split the fasta output of mmseqs cluster into multiple fastas (1 per cluster/orthologuous group).

Arguments:
- first argument is mmseqs fasta file, where separate orthogroups are separated by repeating sequence header/id 
of the first sequence
- second argument (MINSEQS) is the minimum number of sequences the cluster should have for a separate fasta file to be 
generated for it
- third argument is a folder name to record the output fasta files (optional, if not provided the program will save them 
in the directory of execution/current directory)

Note on output fastas indexing: indexes will depend on the value of MINSEQS because only the recorded clusters are counterd
(e.g. the file EntOG0003 will likely correspond to different proteins when executing this function with MINSEQS=10 and
MINSEQS=100) 
"""

import sys
from tqdm import tqdm


if len(sys.argv)==4:
    outfolder=sys.argv[3]
else:
    outfolder=''

minseqs=int(sys.argv[2])

cl=[]
with open(sys.argv[1], 'r') as infasta:
    lines=[]
    ognumid=1
    numseqs=0
    for line in tqdm(infasta):
        lines.append(line)
        if line[:5]=='>GAMA':
            numseqs=numseqs+1
            if len(line)<30:
                ogid=line[:-1]
            elif line[:len(ogid)]==ogid:
                #write sequences to file only if there are enough (>=than the threshold)
                if numseqs>=minseqs+3:
                    with open(outfolder+'EntOG'+str(ognumid).zfill(4)+'.fasta', 'w+') as out:
                        out.write(''.join(lines[1:-2]))
                    #increase index of OG
                    ognumid=ognumid+1
                #reset lines list to the cluster ID and fisrt seq ID
                lines=lines[-2:]
                numseqs=2
    #record final OG
    if numseqs>=minseqs+1:
        with open(outfolder+'EntOG'+str(ognumid).zfill(4)+'.fasta', 'w+') as out:
            out.write(''.join(lines[1:]))