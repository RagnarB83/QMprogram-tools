#!/usr/bin/env python3
# Make 

#Pointsize
pointsize=4000

#Linewidth of bars
linewidth=2

##################################################
# probably unnecessary to change anything below
##################################################

#Lists and defintions
occorbsgrab="unset"
virtorbsgrab="unset"
endocc="unset"
tddftgrab="unset"
tddft="unset"

bands_alpha=[]
bands_beta=[]
virtbands=[]
f=[]
virtf=[]
tddftstates=[]

spinflag="unset"
hftyp="unset"
#Only using alpha orbitals (regardless of RHF or UHF). Handled by break.
import sys
with open(sys.argv[1]) as file:
    for line in file:
        if '%tddft' in line:
            tddft="yes"
        if 'Hartree-Fock type      HFTyp' in line:
            hftyp=line.split()[4]
            print("HF type is", hftyp)
            #if hftyp=="UHF":
        if hftyp == "RHF":
            spinflag="alpha"
        if 'SPIN UP ORBITALS' in line:
            spinflag="alpha"
        if 'SPIN DOWN ORBITALS' in line:
            spinflag="beta"
        if occorbsgrab=="yes":
            endocc=line.split()[1]
            if endocc == "0.0000" :
                occorbsgrab="no"
                #virtorbsgrab="yes"
            else:
                if spinflag=="alpha":
                    #print("adding alpha:", line.split()[3])
                    bands_alpha.append(float(line.split()[3]))
                if spinflag=="beta":
                    bands_beta.append(float(line.split()[3]))
                    #print("adding beta:", line.split()[3])
        if 'NO   OCC          E(Eh)            E(eV)' in line:
            occorbsgrab="yes"

#Grabbing TDDFT sates
if tddft=="yes":
    with open(sys.argv[1]) as file:
        for line in file:
            if tddftgrab=="yes":
                if 'STATE' in line:
                    tddftstates.append(float(line.split()[5]))
                tddftgrab="yes"
            if 'the weight of the individual excitations' in line:
                tddftgrab="yes"

print("Occupied MOs are (eV):")
print(bands_alpha)
print(bands_beta)

# Check for numpy and matplotlib, try to exit gracefully if not found
import imp
import matplotlib.pyplot

try:
    imp.find_module('numpy')
    foundnp = True
except ImportError:
    foundnp = False
try:
    imp.find_module('matplotlib')
    foundplot = True
except ImportError:
    foundplot = False
if not foundnp:
    print("Numpy is required. Exiting")
    sys.exit()
if not foundplot:
    print("Matplotlib is required. Exiting")
    sys.exit()
import numpy as np
import matplotlib.pyplot as plt
import math


fig, ax = plt.subplots()
#Alpha MOs
ax.scatter([1]*len(bands_alpha), bands_alpha, color='blue', marker = '_',  s=pointsize, linewidth=linewidth )

#Beta MOs
if hftyp == "UHF":
    ax.scatter([1.05]*len(bands_beta), bands_beta, color='red', marker = '_',  s=pointsize, linewidth=linewidth )

plt.xlim(0.98,1.07)
plt.ylim(-5,3)
plt.xticks([])
plt.ylabel('MO energy (eV)')

#Vertical line
plt.axhline(y=0.0, color='black', linestyle='--')
plt.savefig(str(sys.argv[1])+'.png', format='png', dpi=200)

plt.show()

