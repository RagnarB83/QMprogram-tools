#!/usr/bin/env python3

#Hessian Tool: Analyze the Hessian in terms of normal mode compositions etc.
# Read ORCA Hessian file, grabs Hessian, diagonalizes and prints out normal mode compositions for all atoms, elements or specific atom
import numpy as np
import sys
from numpy import linalg as la

#Assuming nonlinear molecule. Change to 5 if linear
TRmodenum=6


class bcolors:
    HEADER = '\033[95m' ; OKBLUE = '\033[94m'; OKGREEN = '\033[92m'; WARNING = '\033[93m'; FAIL = '\033[91m'; ENDC = '\033[0m'; BOLD = '\033[1m'; UNDERLINE = '\033[4m'

try:
    hessfile=sys.argv[1]
    option=sys.argv[2]
except IndexError:
    print(bcolors.OKGREEN,"Hess-tool: Normal mode analysis of ORCA-Hessian",bcolors.ENDC)
    print(bcolors.OKGREEN,"Script usage: Hess-tool.py ORCA-Hessianfile Atomgroup", bcolors.ENDC)
    print(bcolors.WARNING,"Atomgroup option: 'all' , 'elements' or  'X' (where X is specific atomnumber)   (numbering starts from 0)",bcolors.ENDC)
    print(bcolors.WARNING,"Example: ./Hess-tool.py file.hess elements", bcolors.ENDC)
    quit()

###############
# FUNCTIONS
##############
#Is integer
def isint(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

#Function to grab masses and elements
def masselemgrab(hessfile):
    grab=False
    elems=[]; masses=[]
    with open(hessfile) as hfile:
        for line in hfile:
            if '$actual_temperature' in line:
                grab=False
            if grab==True and len(line.split()) == 1:
                numatoms=int(line.split()[0])
            if grab==True and len(line.split()) == 5 :
                elems.append(line.split()[0])
                masses.append(float(line.split()[1]))
            if '$atoms' in line:
                grab=True
    return masses, elems,numatoms

#Function to grab Hessian from ORCA-Hessian file
def Hessgrab(hessfile):
    hesstake=False
    j=0
    orcacoldim=5
    shiftpar=0
    lastchunk=False
    with open(hessfile) as hfile:
        for line in hfile:
            #print(len(line))
            if '$vibrational_frequencies' in line:
                hesstake=False
                continue
            if hesstake==True and len(line.split()) == 1:
                hessdim=int(line.split()[0])
                hessarray2d=np.zeros((hessdim, hessdim))
            if hesstake==True and len(line.split()) == 5:
                continue
                #Headerline
            if hesstake==True and lastchunk==True:
                if len(line.split()) == hessdim - shiftpar +1:
                    for i in range(0,hessdim - shiftpar):
                        hessarray2d[j,i+shiftpar]=line.split()[i+1]
                    j+=1
            if hesstake==True and len(line.split()) == 6:
                #Hessianline
                for i in range(0,orcacoldim):
                    hessarray2d[j,i+shiftpar]=line.split()[i+1]
                j+=1
                if j==hessdim:
                    shiftpar+=orcacoldim
                    j=0
                    if hessdim - shiftpar < orcacoldim:
                        lastchunk=True
            if '$hessian' in line:
                hesstake=True
    return hessarray2d

def massweight(matrix):
    mass_mat = np.zeros( (3*numatoms,3*numatoms), dtype = float )
    molwt = [ masses[int(i)] for i in range(numatoms) for j in range(3) ]
    for i in range(len(molwt)):
        mass_mat[i,i] = molwt[i] ** -0.5
    mwhessian = np.dot((np.dot(mass_mat,matrix)),mass_mat)
    return mwhessian,mass_mat

def calcfreq(evalues):
    hartree2j = 4.3597438e-18
    bohr2m = 5.29177208e-11
    amu2kg = 1.66054e-27
    c = 2.99792458e10
    pi = np.pi
    evalues_si = [val*hartree2j/bohr2m/bohr2m/amu2kg for val in evalues]
    vfreq_hz = [1/(2*pi)*np.sqrt(np.complex_(val)) for val in evalues_si]
    vfreq = [val/c for val in vfreq_hz]
    return vfreq

#Give normal mode composition factors for mode j and atom a
def normalmodecomp(evectors,j,a):
    #square elements of mode j
    esq_j=[i ** 2 for i in evectors[j]]
    #Squared elements of atom a in mode j
    esq_ja=[]
    esq_ja.append(esq_j[a*3+0]);esq_ja.append(esq_j[a*3+1]);esq_ja.append(esq_j[a*3+2])
    return sum(esq_ja) 

# Redefines complex number as normal-looking
# Todo: deal with complex
def clean_number(number):
    return np.real_if_close(number)



##########################################
# Main program
##########################################
print(bcolors.OKGREEN,"Hess-tool: Normal mode analysis of ORCA-Hessian",bcolors.ENDC)
if option=="all":
    print(bcolors.OKBLUE,"Option: All atom composition factors",bcolors.ENDC)
    print("")
    print("Normal modes and their composition factors (e**2_ja)")
elif option=="elements":
    print(bcolors.OKBLUE,"Option: Element composition factors",bcolors.ENDC)
    print("")
    print("Normal modes and their element composition factors")
elif isint(option)==True:
    print(bcolors.OKBLUE,"Option: atom", option, "composition only.",bcolors.ENDC)
    print("")
    print("Normal modes and the composition factor for atom:", int(option))
else:
    print(bcolors.FAIL,"Unknown option. Doing all atoms.",bcolors.ENDC)
    option="all"

print("")

#Grab masses, elements and numatoms from Hessianfile
masses,elems,numatoms=masselemgrab(hessfile)
atomlist=[]
for i,j in enumerate(elems):
    atomlist.append(str(j)+'-'+str(i))

#Grab Hessian from Hessianfile
hessian=Hessgrab(hessfile)

#Massweight Hessian
mwhessian,massmatrix =massweight(hessian)

#Diagonalize mass-weighted Hessian
evalues,evectors = la.eigh(mwhessian) 
evectors = np.transpose(evectors)

#Calculate frequencies from eigenvalues
vfreq = calcfreq(evalues)

#Unweight eigenvectors to get normal modes
nmodes = np.dot(evectors,massmatrix) 

#Normalmodecomposition factors for mode j and atom a
if option=="all":
    line = "{}   {}    {}".format("Mode", "Freq(cm**-1)", '     '.join(atomlist))
    print(line)
    for mode in range(0,3*numatoms):
        normcomplist=[]
        if mode < TRmodenum:
            line = "{:2d}      {:2.4f}".format(mode,0.000)
            print(line)
        else:
            vib=clean_number(vfreq[mode])
            for n in range(0,numatoms):
                normcomp=normalmodecomp(evectors,mode,n)
                normcomplist.append(normcomp)
            normcomplist=['{:.6f}'.format(x) for x in normcomplist]
            line = "{:2d}     {:4.4f}        {}".format(mode, vib, '   '.join(normcomplist))
            print(line)
elif option=="elements":
    uniqelems=[]
    for i in elems:
        if i not in uniqelems:
            uniqelems.append(i)
    line = "{}   {}      {}".format("Mode", "Freq(cm**-1)", '        '.join(uniqelems))
    print(line)
    for mode in range(0,3*numatoms):
        normcomplist=[]
        if mode < TRmodenum:
            line = "{:2d}      {:2.4f}".format(mode,0.000)
            print(line)
        else:
            vib=clean_number(vfreq[mode])
            for n in range(0,numatoms):
                normcomp=normalmodecomp(evectors,mode,n)
                normcomplist.append(normcomp)
            elementnormcomplist=[]
            #Sum components together
            for u in uniqelems:
                elcompsum=0.0
                elindices=[i for i, j in enumerate(elems) if j == u]
                for h in elindices:
                    elcompsum=float(elcompsum+float(normcomplist[h]))
                elementnormcomplist.append(elcompsum)
            elementnormcomplist=['{:.6f}'.format(x) for x in elementnormcomplist]
            line = "{:2d}     {:4.4f}        {}".format(mode, vib, '   '.join(elementnormcomplist))
            print(line)
elif isint(option)==True:
    atom=int(option)
    if atom > numatoms-1:
        print(bcolors.FAIL,"Atom index does not exist. Note: Numbering starts from 0",bcolors.ENDC)
        exit()
    line = "{}   {}    {}".format("Mode", "Freq(cm**-1)", atomlist[atom])
    print(line)
    for mode in range(0,3*numatoms):
        normcomplist=[]
        if mode < TRmodenum:
            line = "{:2d}      {:2.4f}".format(mode,0.000)
            print(line)
        else:
            vib=clean_number(vfreq[mode])
            for n in range(0,numatoms):
                normcomp=normalmodecomp(evectors,mode,n)
                normcomplist.append(normcomp)
            normcomplist=['{:.6f}'.format(x) for x in normcomplist]
            line = "{:2d}     {:4.4f}        {}".format(mode, vib, normcomplist[atom])
            print(line)

