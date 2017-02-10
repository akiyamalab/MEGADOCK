#!/usr/bin/python
#
# Masahito Ohue 
# 
# Usage: Contact.py pdb1 pdb2
#

import os
import sys
import math
import Bio
from Bio.PDB import *
from math import *

dist = 4 # contact distance

PDB1 = sys.argv[1]
PDB2 = sys.argv[2]

parser = PDBParser()
str1 = parser.get_structure("1", PDB1)
str2 = parser.get_structure("2", PDB2)
rescon = 0
#chain 1 -> chain 2
for model1 in str1.get_list():
    for chain1 in model1.get_list():
        for r1 in chain1.get_list():
            for a1 in r1.get_list():
                for model2 in str2.get_list():
                    for chain2 in model2.get_list():
                        for r2 in chain2.get_list():
                            for a2 in r2.get_list():
                                if rescon == 1:
                                    continue
                                d = a1 - a2
                                if d < dist:
                                    rescon = 1
                                    
                            if rescon == 1:
                                print r1, r2
                            rescon = 0

