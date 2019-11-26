#!/usr/bin/env python3

#Hessian Tool: Analyze the Hessian in terms of normal mode compositions etc.
# Reads ORCA Hessian file, grabs Hessian, diagonalizes and prints out normal mode compositions for all atoms, elements or specific atom or atomgroup.
# Experimental: printing of VDOS and PVDOS spectra (dat and stk files).
import numpy as np
import os
import sys
from numpy import linalg as la

#Assuming nonlinear molecule. Will be changed to 5 if diatomic. Other linear detection not present
TRmodenum=6


class bcolors:
    HEADER = '\033[95m' ; OKBLUE = '\033[94m'; OKGREEN = '\033[92m'; WARNING = '\033[93m'; FAIL = '\033[91m'; ENDC = '\033[0m'; BOLD = '\033[1m'; UNDERLINE = '\033[4m'

#Check which option chosen
try:
    # Check if doing modecompare on 2 Hessians or if doing normal-mode-comparison analysis of 1 Hessian
    if sys.argv[1]=="Modecompare":
        mainoption="Modecompare"
        hessfileA = sys.argv[2]
        hessfileB = sys.argv[3]
    else:
        mainoption="regular"
        hessfile=sys.argv[1]
        option=sys.argv[2]
except IndexError:
    print(bcolors.OKGREEN,"Hess-tool: Normal mode analysis of ORCA-Hessian",bcolors.ENDC)
    print(bcolors.OKGREEN,"Script usage: Hess-tool.py ORCA-Hessianfile Atomgroup", bcolors.ENDC)
    print(bcolors.OKGREEN,"Other option: Hess-tool.py Modecompare Hessfile1.hess Hessfile2.hess [optional: threshold]", bcolors.ENDC)
    print(bcolors.WARNING,"Atomgroup option: 'all' , 'elements', 'X' (where X is specific atomnumber) or \'X,Y,M\' (where X, Y, M etc. are atom numbers)", bcolors.ENDC)
    print(bcolors.WARNING,"Atom numbering starts from 0.",bcolors.ENDC)
    print(bcolors.WARNING,"Example: ./Hess-tool.py file.hess all", bcolors.ENDC)
    print(bcolors.WARNING,"Example: ./Hess-tool.py file.hess elements", bcolors.ENDC)
    print(bcolors.WARNING,"Example: ./Hess-tool.py file.hess 2", bcolors.ENDC)
    print(bcolors.WARNING,"Example: ./Hess-tool.py file.hess 2,3,8,9", bcolors.ENDC)
    print("")
    print(bcolors.WARNING,"Add \"-VDOS\" as last flag for printing of VDOS spectra",bcolors.ENDC)
    quit()




VDOS=False
try:
    if "VDOS" in sys.argv[3]:
        VDOS=True
except IndexError:
    VDOS=False


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
    grabsize=False
    with open(hessfile) as hfile:
        for line in hfile:
            if '$vibrational_frequencies' in line:
                hesstake=False
                continue
            if hesstake==True and len(line.split()) == 1 and grabsize==True:
                grabsize=False
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
                grabsize=True
    return hessarray2d

def massweight(matrix,masses,numatoms):
    mass_mat = np.zeros( (3*numatoms,3*numatoms), dtype = float )
    molwt = [ masses[int(i)] for i in range(numatoms) for j in range(3) ]
    for i in range(len(molwt)):
        mass_mat[i,i] = molwt[i] ** -0.5
    mwhessian = np.dot((np.dot(mass_mat,matrix)),mass_mat)
    return mwhessian,mass_mat

#
def calcfreq(evalues):
    hartree2j = 4.3597438e-18
    bohr2m = 5.29177210903e-11
    #amu2kg = 1.66054e-27
    amu2kg = 1.66053906660e-27
    #speed of light in cm/s
    c = 2.99792458e10
    pi = np.pi
    evalues_si = [val*hartree2j/bohr2m/bohr2m/amu2kg for val in evalues]
    vfreq_hz = [1/(2*pi)*np.sqrt(np.complex_(val)) for val in evalues_si]
    vfreq = [val/c for val in vfreq_hz]
    return vfreq


#using constants by ORCA
def calcfreqORCA(evalues):
    vfreq = [np.sqrt(eval*0.1602186765*27.2107/0.529177210903/0.529177210903)*1302.78 for eval in evalues]
    return vfreq

#using constants by ORCA
def calcfreqORCAimproved(evalues):
    vfreq = [np.sqrt(eval*0.1602176634*27.211386245988/0.529177210903/0.529177210903)*1302.78 for eval in evalues]
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
    #print("number:", number)
    #print(type(number))
    #print("np.real_if_close(number)", np.real_if_close(number))
    return np.real_if_close(number)

#Can use intensity info. For examples oscillator strengths
def badgaussBand(x, band, strength, broad):
    "Produces a Gaussian curve"
    #print("x is", x)
    #print("band is", band)
    #print("strength is", strength)
    #print("stdev is", stdev)
    bandshape = strength / (broad)  * np.exp(-(((x)-(band))/(broad))**2)
    return bandshape

def Gaussian(x, mu, strength, sigma):
    "Produces a Gaussian curve"
    #print("x:", x,)
    #print( "mu:", mu, "strength:", strength, "sigma:", sigma)
    #print(strength / (sigma*np.sqrt(2*np.pi)))
    #bandshape = (strength / (sigma*np.sqrt(2*np.pi)))  * np.exp(-1*((x-mu))**2/(2*(sigma**2)))
    bandshape = (strength)  * np.exp(-1*((x-mu))**2/(2*(sigma**2)))
    return bandshape


def Lorentzian(x, mu, strength, sigma):
    "Produces a Lorentz curve"
    print("Lorentzian function is not ready yet!")
    exit()
    bandshape = (strength)  * sigma / (((x-u)**2)+(sigma**2))
    return bandshape

def lorentzian(x,c,w):
    """ Analytic Lorentzian function with amplitude 'a', center 'c', width 'w'.
        The FWHM of this fn is 2*w
        NOT NORMALISED """
    L = w**2 / ( (x-c)**2 + w**2 )
    L /= L.max()
    return L

def Voight(x, mu, strength, sigma):
    "Produces a Voight curve"
    print("Voight function is not ready yet")
    exit()
    gamma=sigma
    z=x-mu+i*gamma/(sigma*np.sqrt(2))
    w=np.exp(-1*(z**2))*erfc(-1*i*z)
    bandshape = (strength)  * Re(w)/(sigma*np.sqrt(2*np.pi))
#    return bandshape

#Function to print normal mode composition factors for all atoms, element-groups, specific atom groups or specific atoms
def printnormalmodecompositions(option,TRmodenum):
    #Normalmodecomposition factors for mode j and atom a
    freqs=[]
    #If one set of normal atom compositions (1 atom or 1 group)
    comps=[]
    #If multiple (case: all or elements)
    allcomps=[]

    #Change TRmodenum to 5 if diatomic molecule since linear case
    if numatoms==2:
        TRmodenum=5

    if option=="all":
        #Case: All atoms
        line = "{:>4}{:>14}      {:}".format("Mode", "Freq(cm**-1)", '       '.join(atomlist))
        print(line)
        for mode in range(0,3*numatoms):
            normcomplist=[]
            if mode < TRmodenum:
                line = "{:>3d}   {:>9.4f}".format(mode,0.000)
                print(line)
            else:
                vib=clean_number(vfreq[mode])
                freqs.append(float(vib))
                for n in range(0,numatoms):
                    normcomp=normalmodecomp(evectors,mode,n)
                    normcomplist.append(normcomp)
                allcomps.append(normcomplist)
                normcomplist=['{:.6f}'.format(x) for x in normcomplist]
                line = "{:>3d}   {:>9.4f}        {}".format(mode, vib, '   '.join(normcomplist))
                print(line)
    elif option=="elements":
        #Case: By elements
        uniqelems=[]
        for i in elems:
            if i not in uniqelems:
                uniqelems.append(i)
        line = "{:>4}{:>14}      {:45}".format("Mode", "Freq(cm**-1)", '         '.join(uniqelems))
        print(line)
        for mode in range(0,3*numatoms):
            normcomplist=[]
            if mode < TRmodenum:
                line = "{:>3d}   {:>9.4f}".format(mode,0.000)
                print(line)
            else:
                vib=clean_number(vfreq[mode])
                freqs.append(float(vib))
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
                #print(elementnormcomplist)
                allcomps.append(elementnormcomplist)
                elementnormcomplist=['{:.6f}'.format(x) for x in elementnormcomplist]
                line = "{:>3d}   {:>9.4f}        {}".format(mode, vib, '   '.join(elementnormcomplist))
                print(line)
    elif isint(option)==True:
        #Case: Specific atom
        atom=int(option)
        if atom > numatoms-1:
            print(bcolors.FAIL,"Atom index does not exist. Note: Numbering starts from 0",bcolors.ENDC)
            exit()
        line = "{:>4}{:>14}      {:45}".format("Mode", "Freq(cm**-1)", atomlist[atom])
        print(line)
        for mode in range(0,3*numatoms):
            normcomplist=[]
            if mode < TRmodenum:
                line = "{:>3d}   {:>9.4f}".format(mode,0.000)
                print(line)
            else:
                vib=clean_number(vfreq[mode])
                freqs.append(float(vib))
                for n in range(0,numatoms):
                    normcomp=normalmodecomp(evectors,mode,n)
                    normcomplist.append(normcomp)
                comps.append(normcomplist[atom])
                normcomplist=['{:.6f}'.format(x) for x in normcomplist]
                line = "{:>3d}   {:>9.4f}        {}".format(mode, vib, normcomplist[atom])
                print(line)
    elif len(option.split(",")) > 1:
        #Case: Chemical group defined as list of atoms
        selatoms = option.split(",")
        selatoms=[int(i) for i in selatoms]
        grouplist=[]
        for at in selatoms:
            if at > numatoms-1:
                print(bcolors.FAIL,"Atom index does not exist. Note: Numbering starts from 0",bcolors.ENDC)
                exit()
            grouplist.append(atomlist[at])
        simpgrouplist='_'.join(grouplist)
        grouplist=', '.join(grouplist)
        line = "{}   {}    {}".format("Mode", "Freq(cm**-1)", "Group("+grouplist+")")
        print(line)
        for mode in range(0,3*numatoms):
            normcomplist=[]
            if mode < TRmodenum:
                line = "{:>3d}   {:>9.4f}".format(mode,0.000)
                print(line)
            else:
                vib=clean_number(vfreq[mode])
                freqs.append(float(vib))
                for n in range(0,numatoms):
                    normcomp=normalmodecomp(evectors,mode,n)
                    normcomplist.append(normcomp)
                #normcomplist=['{:.6f}'.format(x) for x in normcomplist]
                groupnormcomplist=[]
                for q in selatoms:
                    groupnormcomplist.append(normcomplist[q])
                comps.append(sum(groupnormcomplist))
                sumgroupnormcomplist='{:.6f}'.format(sum(groupnormcomplist))
                line = "{:>3d}   {:9.4f}        {}".format(mode, vib, sumgroupnormcomplist)
                print(line)
    else:
        print("Something went wrong")

    return allcomps,comps,freqs

#Print vibrational density of states plots based on calculated normal mode composition factors

def printVDOS(option,allcomps,comps,freqs):
    print("");print("")
    print(bcolors.OKGREEN,"VDOS printing ON!",bcolors.ENDC)
    print(bcolors.WARNING,"Note: VDOS feature is experimental.",bcolors.ENDC)
    # Adjust the following three variables to change which area of the spectrum is plotted and number of points used
    # in plotting the curves. Units in cm**-1.
    start=0
    finish=4000.0
    points=10000
    import math

    # broadening in cm**-1
    FWHM = 10
    #Lineshape function options: Gaussian, Lorentzian, Voight
    lineshape=Gaussian
    try:
        start=int(sys.argv[4])
        finish=int(sys.argv[5])
        lineshape=sys.argv[6]
        points=int(sys.argv[7])
        FWHM=int(sys.argv[8])
    except IndexError:
        print(bcolors.OKBLUE,"Using default parameters.",bcolors.ENDC)
        print(bcolors.OKBLUE,"To configure:", "Hess-tool.py Hessianfile Atomgroup -VDOS startxvalue endxvalue Lineshapefunction Numpoints Broadening ",bcolors.ENDC)
        print(bcolors.OKBLUE,"Example:", "Hess-tool.py", hessfile, "all -VDOS 0 4000 Gaussian 10000 10 ",bcolors.ENDC)
        lineshapefunction=Gaussian
    print("")
    if lineshape=="Gaussian":
        lineshapefunction=Gaussian
    elif lineshape=="Lorentzian":
        lineshapefunction=Lorentzian
    elif lineshape=="Voight":
        lineshapefunction=Voight
    print(bcolors.OKBLUE,"Plotting VDOS from", start, "to", finish, "cm**-1",bcolors.ENDC)
    print(bcolors.OKBLUE,"Lineshape function:", lineshapefunction.__name__, "(Options: Gaussian, Lorentzian, Voight)", bcolors.ENDC)
    print(bcolors.OKBLUE,"Number of points:", points, bcolors.ENDC)
    print(bcolors.OKBLUE,"FWHM broadening:", FWHM, "cm**-1",bcolors.ENDC)
    print("")
    #Height of stick for full composition
    stkheight=1.0
    broad=FWHM/(2*np.sqrt(2*np.log(2)))
    superfinalcomps=[]

    if len(allcomps) >0:
        finalcomps=[]
        #Getting list of normalmodecompositions for each mode per atom or element
        for i in range(0,len(allcomps[0])):
            for modecomplist in allcomps:
                finalcomps.append(modecomplist[i])
            superfinalcomps.append(finalcomps)
            finalcomps=[]
    else:
        superfinalcomps.append(comps)

    #X-range for VDOS
    x = np.linspace(start,finish,points)
    pi = np.pi

    #Total VDOS
    totdospeak = 0
    datfile=open(str(basename)+'-TotalVDOS.dat', 'w')
    stkfile=open(str(basename)+'-TotalVDOS.stk', 'w')
    #Not scaling strength here
    strength=stkheight
    for count,f in enumerate(freqs):
        stkfile.write(str(f)+', '+str(stkheight)+'\n')
        dospeak = lineshapefunction(x, f, strength, broad)
        totdospeak += dospeak
    for i in range(0,len(x)):
        datfile.write(str(x[i])+", ")
        datfile.write(str(totdospeak[i])+" \n")
    print(bcolors.OKGREEN,"Total VDOS files created:", str(basename+'-TotalVDOS.dat'), "and", str(basename+'-TotalVDOS.stk'),bcolors.ENDC)
    datfile.close()
    stkfile.close()
    uniqelems = []
    for i in elems:
        if i not in uniqelems:
            uniqelems.append(i)
    #Partial VDOS, per every-atom, every-element, specific atom or group-of-atoms
    for count,subset in enumerate(superfinalcomps):
        if optionmode=="all":
            label=atomlist[count]
        elif optionmode=="elements":
            label=uniqelems[count]
        elif optionmode=="singleatom":
            atom = int(option)
            label=atomlist[atom]
        elif optionmode=="group":
            selatoms = option.split(",")
            selatoms = [int(i) for i in selatoms]
            grouplist = []
            for at in selatoms:
                if at > numatoms - 1:
                    print(bcolors.FAIL, "Atom index does not exist. Note: Numbering starts from 0", bcolors.ENDC)
                    exit()
                grouplist.append(atomlist[at])
            simpgrouplist = '_'.join(grouplist)
            label="Group_"+simpgrouplist
        totdospeak = 0
        datfile=open(str(basename)+'-PVDOS_'+str(label)+'.dat', 'w')
        stkfile=open(str(basename)+'-PVDOS_'+str(label)+'.stk', 'w')
        for f,g in zip(freqs,subset):
            stkfile.write(str(f)+', '+str(g)+'\n')
            strength=g
            dospeak = lineshapefunction(x, f, strength, broad)
            totdospeak += dospeak
        for i in range(0,len(x)):
            datfile.write(str(x[i])+", ")
            datfile.write(str(totdospeak[i])+" \n")
        print(bcolors.OKGREEN,"Partial VDOS files created:", str(basename+'-PVDOS_'+str(label)+'.dat'), "and", str(basename+'-PVDOS_'+str(label)+'.stk'),bcolors.ENDC)
        datfile.close()
        stkfile.close()


##########################################
# Main program
##########################################

if mainoption=="regular":
    print(bcolors.OKGREEN,"Hess-tool: Normal mode analysis of ORCA-Hessian",bcolors.ENDC)
    if option=="all":
        optionmode="all"
        print(bcolors.OKBLUE,"Option: All atom composition factors",bcolors.ENDC)
        print("")
        print("Normal modes and their composition factors (e**2_ja)")
    elif option=="elements":
        optionmode="elements"
        print(bcolors.OKBLUE,"Option: Element composition factors",bcolors.ENDC)
        print("")
        print("Normal modes and their element composition factors")
    elif isint(option)==True:
        optionmode="singleatom"
        print(bcolors.OKBLUE,"Option: atom", option, "composition only.",bcolors.ENDC)
        print("")
        print("Normal modes and the composition factor for atom:", int(option))
    elif len(option.split(",")) > 1:
        optionmode="group"
        selatoms=option.split(",")
        #selatoms=[int(i) for i in selatoms]
        print(bcolors.OKBLUE,"Option: Atom-list:", ', '.join(selatoms), " group composition.",bcolors.ENDC)
        print("")
        print("Normal modes and the composition factor for atom-group:", ', '.join(selatoms))
    else:
        print(bcolors.FAIL,"Unknown option. Doing all atoms.",bcolors.ENDC)
        option="all"
elif mainoption=="Modecompare":
    print(bcolors.OKBLUE,"Modecomparison: Will compare modes from Hessian files for similarity (cosine similarity, i.e. normalized dot product):")
    print(bcolors.OKBLUE, "Hessian-A:", hessfileA, "and Hessian-B:", hessfileB,bcolors.ENDC)


#REGULAR OPTION
if mainoption=="regular":
    basename = os.path.splitext(hessfile)[0]
    #Grab masses, elements and numatoms from Hessianfile
    masses,elems,numatoms=masselemgrab(hessfile)
    atomlist=[]
    for i,j in enumerate(elems):
        atomlist.append(str(j)+'-'+str(i))

    #print(atomlist)
    #Grab Hessian from Hessianfile
    hessian=Hessgrab(hessfile)

    #Massweight Hessian
    mwhessian,massmatrix =massweight(hessian,masses,numatoms)

    #print(massmatrix)
    #Diagonalize mass-weighted Hessian
    evalues,evectors = la.eigh(mwhessian)
    evectors = np.transpose(evectors)


    #Calculate frequencies from eigenvalues

    vfreq = calcfreq(evalues)
    #print("vfreq:", vfreq)
    #vfreqORCA = calcfreqORCA(evalues)
    #print("vfreq from ORCA:", vfreqORCA)
    #vfreqORCAimproved = calcfreqORCAimproved(evalues)
    #print("vfreq from ORCA-improved:", vfreqORCAimproved)
    print("")
    #Unweight eigenvectors to get normal modes
    nmodes = np.dot(evectors,massmatrix)

    #Now print normalmodecompositions depending on which atomgroup was chosen. Return allcomps,comps and freqs for VDOS
    allcomps,comps,freqs=printnormalmodecompositions(option,TRmodenum)
elif mainoption=="Modecompare":
    #Grab masses, elements and numatoms from Hessianfile
    massesA,elemsA,numatomsA=masselemgrab(hessfileA)
    massesB,elemsB,numatomsB=masselemgrab(hessfileB)

    if elemsA != elemsB:
        print("Elemental compositions of Hessian-files differ! That makes no sense!")
        exit()
    numatoms=numatomsA

    print("Masses of Hessian-A:",massesA)
    print("Masses of Hessian-B:",massesB)

    atomlistA=[]
    for i,j in enumerate(elemsA):
        atomlistA.append(str(j)+'-'+str(i))

    print(atomlistA)

    #Grab Hessians from Hessianfiles
    hessianA=Hessgrab(hessfileA)
    hessianB=Hessgrab(hessfileB)

    #Massweight Hessians
    mwhessianA,massmatrixA =massweight(hessianA,massesA,numatoms)
    mwhessianB,massmatrixB =massweight(hessianB,massesB,numatoms)


    #Diagonalize mass-weighted Hessian
    evaluesA,evectorsA = la.eigh(mwhessianA)
    evaluesB,evectorsB = la.eigh(mwhessianB)
    evectorsA = np.transpose(evectorsA)
    evectorsB = np.transpose(evectorsB)

    #Calculate frequencies from eigenvalues
    vfreqA = calcfreq(evaluesA)
    vfreqB = calcfreq(evaluesB)


    print("")
    #Unweight eigenvectors to get normal modes
    nmodesA = np.dot(evectorsA,massmatrixA)
    nmodesB = np.dot(evectorsB, massmatrixB)


    cos_sims=[]

    #Dict to keep track of which mode corresponds to which (key: modeA, value=[modeB,cosinesimilarity])
    # modeA : [modeB, cosine-similarity]
    # Does not work for the TR modes but who cares

dict_modes = {}
for countA,A in enumerate(nmodesA):
    for countB,B in enumerate(nmodesB):
        #Now doing absolute cosine similarity to get rid of sign
        cos_sim = abs(np.dot(A, B) / (np.linalg.norm(A) * np.linalg.norm(B)))
        cos_sims.append(cos_sim)
    #print(cos_sims)
    #max_cos=max(cos_sims)
    #max_cos_pos=cos_sims.index(max_cos)
    #largest_sim_indices=[cos_sims.index(x) for x in sorted(cos_sims, reverse=True)[:3]]
    #largest_sims=sorted( [x for x in cos_sims], reverse=True )[:3]
    #Getting 3 largest cos_sims and their indices from whole mode-comparison
    largest_sims_and_indices=sorted( [(i,x) for (x,i) in enumerate(cos_sims)], reverse=True )[:3]
    #print("largest_sims_and_indices:", largest_sims_and_indices)
    #print(largest_sims_and_indices[0][1])
    #print("largest_sims:", largest_sims)

    #dict_modes[countA] = [max_cos_pos,max_cos]
    dict_modes[countA] = largest_sims_and_indices
    cos_sims = []

#print(dict_modes)
printing="Modemapping"

#Regular printing of modes in each Hessian from 0 to 3N with cos_sim. Assumes 1-to-1 mapping.
if printing=="Regular":
    line = "{:>4}".format("Mode  Freq-A(cm**-1)  Freq-B(cm**-1)    Cosine-similarity")
    print(line)
    for mode in range(0,3*numatoms):
        if mode < TRmodenum:
            line = "{:>3d}   {:>9.4f}       {:>9.4f}".format(mode,0.000,0.000)
            print(line)
        else:
            vibA=clean_number(vfreqA[mode])
            vibB = clean_number(vfreqB[mode])
            cos_sim = np.dot(nmodesA[mode], nmodesB[mode]) / (np.linalg.norm(nmodesA[mode]) * np.linalg.norm(nmodesB[mode]))
            if abs(cos_sim) < 0.9:
                line = "{:>3d}   {:>9.4f}       {:>9.4f}          {:.3f} {}".format(mode, vibA, vibB, cos_sim, "<------" )
            line = "{:>3d}   {:>9.4f}       {:>9.4f}          {:.3f}".format(mode, vibA, vibB, cos_sim )
            print(line)
elif printing=="Modemapping":
    print("For each mode in Hessian A, printing best corresponding mode in Hessian B.")
    print("")
    try:
        threshold=float(sys.argv[4])
    except IndexError:
        threshold=0.9
    print("Printing threshold for alternative matches:",  threshold)
    print ("To change: Hess-tool.py Modecompare file1.hess file2.hess threshold")
    print("")
    line = "{:>4}".format("Mode(A)  Freq-A           Mode(B)   Freq-B    Cosine similarity")
    print(line)
    for mode in range(0,3*numatoms):
        if mode < TRmodenum:
            line = "{:>3d}   {:>9.4f}                   {:>9.4f}".format(mode,0.000,0.000)
            print(line)
        else:
            vibA=clean_number(vfreqA[mode])
            if np.iscomplex(vibA):
                imagA="i"
                vibA = float(vibA.imag)
            else:
                imagA=""

            #For current mode in Hessian-A, determine the corresponding mode in Hessian-B:
            #Best matching mode
            modeB=dict_modes[mode][0][1]
            vibB=clean_number(vfreqB[modeB])
            if np.iscomplex(vibB):
                imagB="i"
                vibB = float(vibB.imag)
            else:
                imagB=""

            cos_simab=dict_modes[mode][0][0]
            #2nd and 3rd best matching modes
            modeB2=dict_modes[mode][1][1]
            #print("----")
            vibB2=clean_number(vfreqB[modeB2])
            #print("vfreqB[modeB2]:", vfreqB[modeB2])
            #print("vibB2:", vibB2)
            #print(type(vibB))
            cos_simab2 = dict_modes[mode][1][0]
            #exit()
            modeB3=dict_modes[mode][2][1]
            vibB3=clean_number(vfreqB[modeB3])
            cos_simab3 = dict_modes[mode][2][0]

            #print("cos_simab:", cos_simab)
            #print("cos_simab2:", cos_simab2)
            #print("cos_simab3:", cos_simab3)
            #print("modeB:",modeB)
            #print("vibB",vibB)
            #vibB = clean_number(vfreqB[mode])

            if cos_simab3 > threshold:
                if mode == modeB:
                    line = "{:>3d}   {:>9.3f}{}   {}  {} {:>3d} {}  {:>9.3f}{}          {:.3f}     {} {}  {:.3f}   {} {:>}   {:>} {:.3f} {} {}  {:.3f}   {} {:>}   {:>} {:.3f} {}".format(
                        mode, vibA, imagA, " ===>", bcolors.OKBLUE, modeB, bcolors.ENDC, vibB, imagB, cos_simab, "Other matches:", "( Cos.sim:", cos_simab2,
                        "Mode:", modeB2, "Freq:", vibB2, "cm**-1)", "( Cos.sim:", cos_simab3, "Mode:", modeB3, "Freq:",
                        vibB3, "cm**-1)")

                else:
                    line = "{:>3d}   {:>9.3f}{}   {}  {} {:>3d} {}  {:>9.3f}{}          {:.3f}     {} {}  {:.3f}   {} {:>}   {:>} {:.3f} {} {}  {:.3f}   {} {:>}   {:>} {:.3f} {}".format(
                        mode, vibA, imagA, " ===>", bcolors.WARNING, modeB, bcolors.ENDC, vibB, imagB, cos_simab, "Other matches:", "( Cos.sim:", cos_simab2,
                        "Mode:", modeB2, "Freq:", vibB2, "cm**-1)", "( Cos.sim:", cos_simab3, "Mode:", modeB3, "Freq:",
                        vibB3, "cm**-1)")

            elif cos_simab2 > threshold:
                if mode == modeB:
                    line = "{:>3d}   {:>9.3f}{}   {}  {} {:>3d} {}  {:>9.3f}{}          {:.3f}     {} {}  {:.3f}   {} {:>}   {:>} {:.3f} {}".format(
                        mode, vibA, imagA, " ===>", bcolors.OKBLUE, modeB, bcolors.ENDC, vibB, imagB, cos_simab, "Other matches:", "( Cos.sim:", cos_simab2,
                        "Mode:", modeB2, "Freq:", vibB2, "cm**-1)")

                else:
                    line = "{:>3d}   {:>9.3f}{}   {}  {} {:>3d} {}  {:>9.3f}{}          {:.3f}     {} {}  {:.3f}   {} {:>}   {:>} {:.3f} {}".format(
                        mode, vibA, imagA, " ===>", bcolors.WARNING, modeB, bcolors.ENDC, vibB, imagB, cos_simab, "Other matches:", "( Cos.sim:", cos_simab2,
                        "Mode:", modeB2, "Freq:", vibB2, "cm**-1)")
            else:
                if mode == modeB:
                    line = "{:>3d}   {:>9.3f}{}   {}  {} {:>3d} {}  {:>9.3f}{}          {:.3f}".format(mode, vibA, imagA, " ===>", bcolors.OKBLUE, modeB, bcolors.ENDC, vibB, imagB, cos_simab )
                else:
                    line = "{:>3d}   {:>9.3f}{}   {}  {} {:>3d} {}  {:>9.3f}{}          {:.3f}".format(mode, vibA, imagA, " ===>", bcolors.WARNING, modeB, bcolors.ENDC, vibB, imagB, cos_simab )
            print(line)
            #exit()
print("")

if VDOS==True:
    printVDOS(option,allcomps,comps,freqs)
else:
    print(bcolors.WARNING, "No VDOS spectra requested. Use \"-VDOS\" option if wanted.", bcolors.ENDC)
