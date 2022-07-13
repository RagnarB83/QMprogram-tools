#!/bin/env python3
import sys
import os
import numpy as np
import re
import math

#Elements of interest
#syselements=['Pd']
syselements=[sys.argv[2]]
# Threshold for including orbitals
#popthreshold=0.1
popthreshold=float(sys.argv[3])

print("Finding localized orbitals for element:", syselements, "with threshold:", popthreshold)


##################

dict_alpha={}
dict_beta={}
#MAIN program
try:
    args=sys.argv
    file=args[1]
except IndexError:
    print("Script usage:")
    print("orbloc-find.py ORCA-outputfile (for simple usage)")
    quit()


stronggrab=False
bondgrab=False
delocgrab=False
with open(file) as f:
    for line in f:
        if 'More delocalized orbitals:' in line:
            #print(line)
            #print("deloc switch")
            stronggrab=False
            bondgrab=False
            delocgrab=True
        if 'ORCA ORBITAL LOCALIZATION' in line:
            loc=True
        if 'Operator                                 ... 0' in line:
            operator='alpha'
        if 'Operator                                 ... 1' in line:
            operator='beta'
        if 'Rather strongly localized orbitals:' in line:
            stronggrab=True
        if stronggrab == True:
            for el in syselements:
                if el in line:
                    #print(line)
                    if operator == 'alpha':
                        atom=line.split()[2]
                        monumber=line.split()[1][:-1]
                        dict_alpha.setdefault(atom, []).append(int(monumber))
                    elif operator == 'beta':
                        atom=line.split()[2]
                        monumber=line.split()[1][:-1]
                        dict_beta.setdefault(atom, []).append(int(monumber))
                    else:
                        print("neither wtf")
                        exit()
        if bondgrab == True:
            for el in syselements:
                if el in line:
                    #print(line)
                    for i in line.split():
                        if el in i:
                            at=i
                            #print("at is", at)
                            pos=line.split().index(at)
                            #print("pos is", pos)
                            #print("float(line.split()[pos+2]) is", float(line.split()[pos+2]))
                            if float(line.split()[pos+2]) > popthreshold:
                                if operator == 'alpha':
                                    atom=line.split()[pos]
                                    monumber=line.split()[1][:-1]
                                    dict_alpha.setdefault(atom, []).append(int(monumber))
                                elif operator == 'beta':
                                    atom=line.split()[pos]
                                    monumber=line.split()[1][:-1]
                                    #print("atom is", atom, "and monumber is", monumber)
                                    dict_beta.setdefault(atom, []).append(int(monumber))
                                else:
                                    print("neither wtf")
                                    exit()
        if delocgrab == True:
            if 'More delocalized orbita' not in line:
                for el in syselements:
                    if el in line:
                        #print(line)
                        linechanged=line.replace("-","")
                        #print(linechanged)
                        for i in linechanged.split():
                            if el   in i:
                                at=i
                                #print("at is", at)
                                pos=linechanged.split().index(at)
                                #print("pos is", pos)
                                #print("float(line.split()[pos+2]) is", float(line.split()[pos+2]))
                                if float(linechanged.split()[pos+1]) > popthreshold:
                                    if operator == 'alpha':
                                        atom=linechanged.split()[pos]
                                        monumber=linechanged.split()[1][:-1]
                                        dict_alpha.setdefault(atom, []).append(int(monumber))
                                    elif operator == 'beta':
                                        atom=linechanged.split()[pos]
                                        monumber=linechanged.split()[1][:-1]
                                        dict_beta.setdefault(atom, []).append(int(monumber))
                                    else:
                                        print("neither wtf")
                                        exit()
        if 'Bond-like localized orbitals:' in line:
            stronggrab=False
            bondgrab=True
        if 'Localized MO\'s were stored in:' in line:
            stronggrab=False
            bondgrab=False
            delocgrab=False

#print("Alpha orbitals")
#print(dict_alpha)

#print("Beta orbitals")
#print(dict_beta)

alphalist=[]
betalist=[]
for avals in dict_alpha.items():
    j=sorted(avals[1])
    #Deleting metal s and p orbitals
    #j.pop(0);j.pop(0);j.pop(0);j.pop(0)
    for i in j:
        alphalist.append(i)

for bvals in dict_beta.items():
    j=sorted(bvals[1])
    #Deleting metal s and p orbitals
    #j.pop(0);j.pop(0);j.pop(0);j.pop(0)
    for i in j:
        betalist.append(i)

alphalist=sorted(list(set(alphalist)))
betalist=sorted(list(set(betalist)))
print("")

print("Alpha orbitals to be plotted with getorbitals:")
print(*alphalist, sep=' ')
print("Beta orbitals to be plotted with getorbitals:")
print(*betalist, sep=' ')
