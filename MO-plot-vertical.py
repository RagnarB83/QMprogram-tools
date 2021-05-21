#!/usr/bin/env python3
# MO-plot-vertical

# Option to print only occupied or occupied and verticals.
#OccOnly=True  only occ. OccOnly=False  both coc and vit
OccOnly=False

#Pointsize
pointsize=4000

#Linewidth of bars
linewidth=2

##################################################
# probably unnecessary to change anything below
##################################################


#Lists and defintions
occorbsgrab=False
virtorbsgrab=False
endocc="unset"
tddftgrab="unset"
tddft="unset"

bands_alpha=[]
bands_beta=[]
virtbands_a=[]
virtbands_b=[]
f=[]
virtf=[]
tddftstates=[]

spinflag="unset"
hftyp="unset"
#Only using alpha orbitals (regardless of RHF or UHF). Handled by break.
import sys
with open(sys.argv[1]) as file:
    for line in file:
        #print("spinflag:",spinflag)
        #print("occorbsgrab:", occorbsgrab, "virtorbsgrab:", virtorbsgrab)
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
        if occorbsgrab==True:
            endocc=line.split()[1]
            if endocc == "0.0000" :
                occorbsgrab=False
                virtorbsgrab=True
            else:
                if spinflag=="alpha":
                    bands_alpha.append(float(line.split()[3]))
                if spinflag=="beta":
                    bands_beta.append(float(line.split()[3]))
        if virtorbsgrab==True:
            print("line2:", line)
            if '------------------' in line:
                break
            if line == '\n':
                print("Setting virtorbs to false", line)
                virtorbsgrab=False
                spinflag="unset"
                continue
            if spinflag=="alpha":
                virtbands_a.append(float(line.split()[3]))
            if spinflag=="beta":
                virtbands_b.append(float(line.split()[3]))
            endvirt=line.split()[1]
        if 'NO   OCC          E(Eh)            E(eV)' in line:
            occorbsgrab=True


print("Occupied MOs-alpha are (eV):")
print(bands_alpha)
print("")
print("Occupied MOs-beta are (eV):")
print(bands_beta)
print("")
print("virtual alpha")
print(virtbands_a)
print("")
print("virtual beta")
print(virtbands_b)
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


if OccOnly==True:
    print("OccOnly True! Plotting only occupied orbital levels")
else:
    print("OccOnly False! Plotting occupied and virtual orbital levels")

fig, ax = plt.subplots()
#Alpha MOs
ax.scatter([1]*len(bands_alpha), bands_alpha, color='blue', marker = '_',  s=pointsize, linewidth=linewidth )
if OccOnly!=True:
    ax.scatter([1]*len(virtbands_a), virtbands_a, color='cyan', marker = '_',  s=pointsize, linewidth=linewidth )

#Beta MOs
if hftyp == "UHF":
    ax.scatter([1.05]*len(bands_beta), bands_beta, color='red', marker = '_',  s=pointsize, linewidth=linewidth )
    if OccOnly!=True:
        ax.scatter([1.05]*len(virtbands_b), virtbands_b, color='pink', marker = '_',  s=pointsize, linewidth=linewidth )

plt.xlim(0.98,1.07)
plt.ylim(-12,3)
plt.xticks([])
plt.ylabel('MO energy (eV)')

#Vertical line
plt.axhline(y=0.0, color='black', linestyle='--')
plt.savefig(str(sys.argv[1])+'.png', format='png', dpi=200)



plt.show()

