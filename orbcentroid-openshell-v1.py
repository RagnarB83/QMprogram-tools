#!/bin/env python3
import sys
import os
import numpy as np
import re
import math

#Oxidation state definition: fullvalence or dvalence
#If set to fullvalence then localized orbital cubefiles of all valence electrons have to be created, e.g. all 3s3p3d for Fe.
#If set to dvalence then only metal d-orbitals are assumed to be present and oxidation state is determined by counting d-electrons.
oxstatedef='dvalence'

#List of transition metals with d-electrons
tmlist=['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg']

bohrang=0.529177

elements=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']

oxnumbers={-1:'-I',0:'0',-2:'-II',-3:'-III',-4:'-IV',-5:'-V',-6:'-VI',-7:'-VII',-8:'-VIII',-9:'-IX',-10:'-X',1:'I',2:'II',3:'III',4:'IV',5:'V',6:'VI',7:'VII',8:'VIII',9:'IX',10:'X'}

#Element dictionary: Key: atomic number. Value: List of [Atomsymbol,covalent radius,valence electrons]
#Valence electrons: For p-block including whole electron shell (s,p). For TM: including s,p,d (and next s if applicable)
#Covalent radii from Alvarez
#Note: Carbon requires special treatment: sp3, sp2, sp. Sp3 in dict
#Note: Mn ls and hs options. hs in dict
#Note: Fe ls and hs options. hs in dict
#Note Co ls and hs options: hs in dict
#Elements added: H-Cd. Rest to be done
eldict={1:['H',0.31,1],2:['He',0.28,2],3:['Li',1.28,1],4:['Be',0.96,2],5:['B',0.84,3],6:['C',0.76,4],7:['N',0.71,5],8:['O',0.66,6],9:['F',0.57,7],10:['Ne',0.58,8],11:['Na',1.66],12:['Mg',1.41],13:['Al',1.21],14:['Si',1.11],15:['P',1.07],16:['S',1.05],17:['Cl',1.02,7],18:['Ar',1.06],19:['K',2.03],20:['Ca',1.76],21:['Sc',1.70],22:['Ti',1.6],23:['V',1.53],24:['Cr',1.39],25:['Mn',1.61],26:['Fe',1.52,16],27:['Co',1.50],28:['Ni',1.24],29:['Cu',1.32],30:['Zn',1.22],31:['Ga',1.22],32:['Ge',1.20],33:['As',1.19],34:['Se',1.20],35:['Br',1.20],36:['Kr',1.16],37:['Rb',2.2],38:['Sr',1.95],39:['Y',1.9],40:['Zr',1.75],41:['Nb',1.64],42:['Mo',1.54,14],43:['Tc',1.47],44:['Ru',1.46],45:['Rh',1.42],46:['Pd',1.39,18],47:['Ag',1.45],48:['Cd',1.44]}

#D-valence electrons for some elements (counting atomic s-electrons as d)
delectrondict={26:['Fe',8],42:['Mo',6],23:['V',5]}

def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

def read_cube (cubefile):
    global bohrang
    #Opening orbital cube file
    try:
        filename = cubefile
        a = open(filename,"r")
        print("Reading orbital file:", filename)
        filebase=os.path.splitext(filename)[0]
    except IndexError:
        print("error")
        quit()
    if denswrite==True:
        #Write orbital density cube file
        output = open(filebase+'-dens.cube', "w")
    #Read cube file and get all data. Square values
    count = 0
    X = False
    d = []
    densvals = []
    orbvals=[]
    elems=[]
    molcoords=[]
    molcoords_ang=[]
    numatoms=0
    for line in a:
        count += 1
        words = line.split()
        numwords=len(words)
        #Grabbing origin
        if count < 3:
            if denswrite==True:
                output.write(line) 
        if count == 3:
            numatoms=abs(int(line.split()[0]))
            orgx=float(line.split()[1])
            orgy=float(line.split()[2])
            orgz=float(line.split()[3])
            rlowx=orgx;rlowy=orgy;rlowz=orgz
            if denswrite==True:
                output.write(line)
        if count == 4:
            nx=int(line.split()[0])
            dx=float(line.split()[1])
            if denswrite==True:
                output.write(line)
        if count == 5:
            ny=int(line.split()[0])
            dy=float(line.split()[2])
            if denswrite==True:
                output.write(line)
        if count == 6:
            nz=int(line.split()[0])
            dz=float(line.split()[3])
            if denswrite==True:
                output.write(line)
        #Grabbing molecular coordinates
        if count > 6 and count <= 6+numatoms:
            elems.append(int(line.split()[0]))
            molcoord=[float(line.split()[2]),float(line.split()[3]),float(line.split()[4])]
            molcoord_ang=[bohrang*float(line.split()[2]),bohrang*float(line.split()[3]),bohrang*float(line.split()[4])]
            molcoords.append(molcoord)
            molcoords_ang.append(molcoord_ang)
            if denswrite==True:
                output.write(line)
        # reading gridpoints
        if X == True:
            b = line.rstrip('\n').replace('  ', ' ').replace('  ', ' ').split(' ')
            b=list(filter(None, b))
            c =[float(i) for i in b]
            #print("c is", c)
            #Squaring orbital values to get density
            csq = [q** 2 for q in c]
            dsq = [float('%.5e' % i) for i in csq]
            densvals.append(dsq)
            dbq = [float('%.5e' % i) for i in c]
            orbvals.append(dbq)
        # when to begin reading gridpoints
        if (count > 6 and numwords == 2 and X==False):
            X = True
            if denswrite==True:
                output.write(line)

    # Go through orb and dens list and print out density file
    alldensvalues=[]
    allorbvalues=[]
    for line in densvals:
        columns = ["%13s" % cell for cell in line]
        for val in columns:
            alldensvalues.append(float(val))
        if denswrite==True:
            linep=' '.join( columns)
            output.write(linep+'\n')

    for line in orbvals:
        dolumns = ["%13s" % cell for cell in line]
        for oval in dolumns:
            allorbvalues.append(float(oval))
    if denswrite==True:
        output.close()
        print("Wrote orbital density file as:", filebase+'-dens.cube')
        print("")
    sumdensvalues=sum(i for i in alldensvalues)
    if LargePrint==True:
        print("Sum of density values is:", sumdensvalues)
        print("Number of density values is", len(alldensvalues))
        print("Number of orb values is", len(allorbvalues))
    return rlowx,dx,nx,orgx,rlowy,dy,ny,orgy,rlowz,dz,nz,orgz,alldensvalues,elems,molcoords_ang,numatoms,filebase


def centroid_calc (rlowx,dx,nx,orgx,rlowy,dy,ny,orgy,rlowz,dz,nz,orgz,alldensvalues ):
    #########################################################
    # Calculate centroid. Based on Multiwfn
    ############################################################

    #Largest x,y,z coordinates
    rhighx=rlowx+(dx*(nx-1))
    rhighy=rlowy+(dy*(ny-1))
    rhighz=rlowz+(dz*(nz-1))
    #Lowest and highest density values
    rlowv = min(float(s) for s in alldensvalues)
    rhighv = max(float(s) for s in alldensvalues)

    sumuppos=0.0
    cenxpos=0.0
    cenypos=0.0
    cenzpos=0.0
    vcount=0

    #print ("dx, dy, dz is", dx, dy, dz)
    #print("range of x:", rlowx, rhighx)
    #print("range of y:", rlowy, rhighy)
    #print("range of z:", rlowz, rhighz)

    for i in range(1,nx+1):
        if (orgx+(i-1)*dx)<rlowx or (orgx+(i-1)*dx)>rhighx:
            print("If statement. Look into. x")
            exit()
            continue
        for j in range(1,ny+1):
            if (orgy+(j-1)*dy)<rlowy or (orgy+(j-1)*dy)>rhighy:
                print("If statement. Look into. y")
                exit()
                continue
            for k in range(1,nz+1):
                if (orgz+(k-1)*dz)<rlowz or (orgz+(k-1)*dz)>rhighz:
                    print("If statement. Look into. z")
                    exit()
                    continue
                #print("i,j,k is", i,j,k)
                valtmp=alldensvalues[vcount]
                if valtmp<rlowv or valtmp>rhighv:
                    print("If statement. Look into. v")
                    exit()
                    continue
                if valtmp>0:
                    sumuppos=sumuppos+valtmp
                    #print("sumuppos is", sumuppos)
                    cenxpos=cenxpos+(orgx+(i-1)*dx)*valtmp
                    cenypos=cenypos+(orgy+(j-1)*dy)*valtmp
                    cenzpos=cenzpos+(orgz+(k-1)*dz)*valtmp
                    #print("valtmp is", valtmp)
                    #print("-----------------------")
                vcount+=1

    #Final values
    cenxpos=cenxpos/sumuppos
    cenypos=cenypos/sumuppos
    cenzpos=cenzpos/sumuppos
    return cenxpos,cenypos,cenzpos


def centroid_assign (centroids,molcoords_ang,scalingpar,elemlist,mospinlist):
    #print("centroids is", centroids)
    #print("molcoords_ang is", molcoords_ang)
    #print("elemlist is", elemlist)
    rlist=[]
    centroiddict={};atnums=[]
    X=[];Y=[]
    numcentroidsassigned=0
    for ccounter,centroid in enumerate(centroids):
        #print("ccounter is", ccounter)
        #print("centroid is", centroid)
        for atom in molcoords_ang:
            r=math.sqrt((centroid[0]-atom[0])**2+(centroid[1]-atom[1])**2+(centroid[2]-atom[2])**2)
            rlist.append(r)

        #Multi-atom assignment test
        #In this for loop we need to pick the correct covalent radius
        for n,rval in enumerate(rlist):
            #print("rlist is", rlist)
            #index=rlist.index(rval)
            currel=elemlist[n]
            curratomnum=elements.index(currel)+1
            covalradius=eldict[curratomnum][1]
            threshold=covalradius*scalingpar
            #print("rval is", rval)
            if rval < threshold :
                #print("Yes.........")
                X.append(rval)
                rminindexB=rlist.index(rval)
                Y.append(rminindexB)
                if LargePrint==True:
                    print("Centroid distance to nearest atom (", str(rminindex)+currel,") is:", rmin, "and current threshold (based on covalent radius) is: ", threshold)
                moldict.setdefault(rminindexB, []).append(ccounter)
                numcentroidsassigned+=1

        rmin=0
        rminindex=0;currel='';curratomnum=0;covalradius=0;threshold=0
        rminindexB
        index=None
        rlist=[]
        X=[]; Y=[]
        #print("---------------_")
        #print(moldict)

    electronlist=[0] * numatoms
    orbshare=[]
    if 'b' in mospinlist:
        elperorb=1.0
    else:
        elperorb=2.0        
    #Assign electrons to atoms based on if multiple atoms own a centroid
    for centroid in range(0, len(centroids)):
        for atom,k in sorted(moldict.items()):
            if centroid in moldict[atom]:
                #print("yes. Atom :", atom, "has centroid:", centroid)
                orbshare.append(atom)
                atnums.append(atom)
        for bla in orbshare:
            electronlist[bla]=electronlist[bla]+elperorb/len(orbshare)
        orbshare=[]
        centroiddict[centroid] = atnums
        atnums=[]


    return moldict,numcentroidsassigned,electronlist,centroiddict


def writexyz_single (cenxpos, cenypos, cenzpos, numatoms,molcoords_ang,filebase):
    global bohrang
    #Printing xyzfile with molecule coordinates and centroid
    xyzfile = open(filebase+'-centroid.xyz', "w")
    xyzfile.write(str(numatoms+1)+'\n')
    xyzfile.write('molecule + centroid coordinate\n')
    for el,co in zip(elems,molcoords_ang):
        m=' '.join(map(str,co))
        ele=elements[el-1]
        xyzfile.write(ele+'  '+m+'\n')
    xyzfile.write('He  '+str(bohrang*cenxpos)+' '+str(bohrang*cenypos)+' '+str(bohrang*cenzpos)+'\n')
    xyzfile.close()
    print("Wrote molecule xyzfile with orbital-density centroid (as X):", filebase+'-centroid.xyz')

def writexyz_mult (numatoms,molcoords_ang,centroids,elemlist,monumlist,mospinlist):
    global bohrang
    #Printing xyzfile with molecule coordinates and centroid
    xyzfile = open('Mol-centroids.xyz', "w")
    xyzfile.write(str(numatoms+len(centroids))+'\n')
    xyzfile.write('molecule + centroid coordinates\n')
    for el,co in zip(elemlist,molcoords_ang):
        m=' '.join(map(str,co))
        #ele=elements[el-1]
        #elemlist.append(ele)
        xyzfile.write(el+'  '+m+'\n')
    for cent,spin in zip(centroids,mospinlist):
        if spin == 'a':
            xyzfile.write('X  '+str(cent[0])+' '+str(cent[1])+' '+str(cent[2])+'\n')
        else:
            xyzfile.write('He  '+str(cent[0])+' '+str(cent[1])+' '+str(cent[2])+'\n')
    xyzfile.close()
    molinfo = open('Mol-info.txt', "w")
    molinfo.write('elemlist='+str(elemlist)+'\n')
    molinfo.write('centroids='+str(centroids)+'\n')
    molinfo.write('molcoords_ang='+str(molcoords_ang)+'\n')
    molinfo.write('monumlist='+str(monumlist)+'\n')
    molinfo.write('mospinlist='+str(mospinlist)+'\n')
    molinfo.close()
    print("Wrote Mol-centroids.xyz with orbital-density centroids for all cubefiles (as X symbols):")
    print("Wrote Mol-info.txt file")
    print("")




#MAIN program
try:
    args=sys.argv
    print("Arguments:", args[1:])
    arg1=args[1]
    if '-denswrite' in args:
        denswrite=True
    else:
        denswrite=False
    if '-largeprint' in args:
        LargePrint=True
    else:
        LargePrint=False
except IndexError:
    print("Script usage:") 
    print("orbcentroid.py -simple (for simple usage)")
    print("orbcentroid.py -denswrite  (if cube density files are wanted)")
    print("orbcentroids.py -largeprint  (for more printing)")
    quit()



#Checking if Mol-info.txt already exists. Skip orb density procedure if so and read Mol-info.xt instead and proceed.
if os.path.isfile("Mol-info.txt") == False:
    # Finding centroid for each cubefile in current dir.
    centroids=[]
    cubefiles=[]
    monumlist=[]
    mospinlist=[]
#Reading each cubefile, calculating orbital density and then centroid position
    for cubefile in natural_sort(os.listdir()):
        if cubefile.endswith(".cube"):
            cubefiles.append(cubefile)
            monumlist.append(int(cubefile.split('.')[-2][2:-1]))
            mospinlist.append(cubefile.split('.')[-2][-1])
            rlowx,dx,nx,orgx,rlowy,dy,ny,orgy,rlowz,dz,nz,orgz,alldensvalues,elems,molcoords_ang,numatoms,filebase = read_cube(cubefile)
            cenxpos,cenypos,cenzpos = centroid_calc(rlowx,dx,nx,orgx,rlowy,dy,ny,orgy,rlowz,dz,nz,orgz,alldensvalues)
            print("")
            if LargePrint==True:
                print("Orbital-density centroid (Bohrs):", cenxpos, cenypos, cenzpos)
                print("Orbital-density centroid (Angstrom):", bohrang*cenxpos, bohrang*cenypos, bohrang*cenzpos)
            centroids.append([bohrang*cenxpos,bohrang*cenypos,bohrang*cenzpos])
            print("----------------------------------------")
            elemlist=[elements[el-1] for el in elems]
    print("")
    print("Centroid (Angstrom)for each Cubefile")
    for i,j in zip(cubefiles,centroids):
        print("File: ", i, "  Centroid:  ", j)
else:
    print("Mol-info.txt exists. Skipping cube analysis... (delete file to redo analysis)")
    with open('Mol-info.txt') as molinfo:
        exec(molinfo.read())
    numatoms=len(elemlist)



#Write xyzfile containing molecule and all centroids
writexyz_mult(numatoms,molcoords_ang,centroids,elemlist,monumlist,mospinlist)

############################
#Oxidation state assignment#
############################
#This number extends or shortens covalent radius by scaling 
#1.26 
scalingpar=0.9
molnumbegin=monumlist[0]

#print("molnumbegin is", molnumbegin)
#print("monumlist0 is", monumlist[0])
#print("monumlist is", monumlist)

#This molecule dictionary will contain the molecule atom index as key and then a list (as value) of centroid indices.
moldict={}


print("Assigning localized orbital centroids to atoms in molecule...")
print("Using list of covalent radii and scaling parameter", scalingpar) 
moldict,numcentroidsassigned,electronlist,centroiddict=centroid_assign(centroids,molcoords_ang,scalingpar,elemlist,mospinlist)
print("")
print("moldict is", moldict)
print("Centroiddict is", centroiddict)
print("")
print("Number of centroids calculated: ", len(centroids))
print("Number of centroids assignments:", numcentroidsassigned) 
print("Number of shared centroids, a.k.a. bonding orbitals:", numcentroidsassigned-len(centroids))
print("")


for mold,val in sorted(moldict.items()) :
    el=elemlist[mold]
    print("Atom", mold, "(", el,")"," owns centroids:", moldict[mold], " or MOs:", ' , '.join([ str(monumlist[j])+str(mospinlist[j]) for j in moldict[mold]]) )
    spinlist=[str(mospinlist[j]) for j in moldict[mold]]
    alphanum=spinlist.count('a')
    betanum=spinlist.count('b')
    
    print("Alpha-electrons:", alphanum)
    print("Beta-electrons:", betanum)
    print("Crude valence electron-count (fractional numbers indicate bonding):", electronlist[mold])
    if oxstatedef=='dvalence':
        if el in tmlist:
            delectrons=alphanum+betanum
            print("Whole d-electrons:",delectrons)
            print(el, str(alphanum)+'α', str(betanum)+'β')
            oxstate=delectrondict[elements.index(el)+1][1] - delectrons
            print("Oxidation state: ", el, oxnumbers[oxstate])
    else:
        print("Oxidation state: ", el, eldict[elements.index(el)+1][2] - electronlist[mold])
    print("-------------------")


print("")
print("")
print("List of shared centroids (chemical bonds, delocalized electrons etc.):")

print("")
for centr in centroiddict:
    if len(centroiddict[centr]) > 1:
        #print("Centroid", centr, "belongs to atoms:", centroiddict[centr])
        atomsoncentroid=[str(elemlist[j])+' '+str(j) for j in centroiddict[centr]]
        print("MO", str(monumlist[centr])+str(mospinlist[centr]), "belongs to atoms:", ' and '.join(atomsoncentroid))

print("")
