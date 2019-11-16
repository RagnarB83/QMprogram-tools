#!/bin/env python3
#import os.path
import os
import subprocess
import math
import glob
import sys
import numpy as np

# Read ORCA-CPCM file: filename.cpcm and create xyz file for visualization of both molecule and the surface charge coordinates

molcartgrab=False
surfacecartgrab=False
molcart=[]
surfacecart=[]
bohrang=0.529177
filename=sys.argv[1]
filebase=filename.strip('.cpcm')
with open(filename) as file:
    print("Reading CPCM file,  getting molecular and surface charge coordinates ", sys.argv[1])
    for line in file:
        if molcartgrab==True:
            if '#------' not in line:
                if '# SURFACE POINTS' not in line:
                    #print(len(line))
                    if len(line) > 1:
                        #print(line)
                        cart_x=bohrang*float(line.split()[0]);cart_y=bohrang*float(line.split()[1]);cart_z=bohrang*float(line.split()[2])
                        molcart.append([cart_x,cart_y,cart_z])
        if surfacecartgrab==True:
            if '#------' not in line:
                if '# SURFACE POINTS' not in line:
                    #print(len(line))
                    if 'X' not in line:
                        if len(line) > 1:
                            #print(line)
                            cart_x=bohrang*float(line.split()[0]);cart_y=bohrang*float(line.split()[1]);cart_z=bohrang*float(line.split()[2])
                            surfacecart.append([cart_x,cart_y,cart_z])
        if '# Number of atoms' in line:
            numatoms=int(line.split()[0])
        if '# Number of surface points' in line:
            numsurfacecharges=int(line.split()[0])
        if 'CARTESIAN COORDINATES' in line:
            molcartgrab=True

        if '# SURFACE POINT' in line:
            molcartgrab=False
            surfacecartgrab=True

#
print()
print("Getting molecule element information")
elementlist=[]
ocartgrab=False
CartFlag=2
print("Trying to find XYZfile:", str(filebase)+".xyz")
try:
    xyzfile=str(filebase)+".xyz"
    with open(xyzfile) as xfile:
        for i,line in enumerate(xfile):
            if i >1:
                elementlist.append(line.split()[0])
except FileNotFoundError:
    print("Not found. Trying outputfile")
    outfilelist=glob.glob('new-cpcm-snapA-78000-A-chargeB-rmin1_0.*out*')
    print("outfilelist is", outfilelist)
    if len(outfilelist) > 1:
        print("problem. Many outputfiles in dir matching basename:", outfilelist, "Remove extra ones from dir")
        exit()
    else:
        outfile=outfilelist[0]
        print("Reading element list from outputfile: ", outfile) 
        with open(outfile) as ofile:
            for line in ofile:
                if ocartgrab == True:
                    if len(line) < 2:
                        break
                    if '---' not in line:
                        #print(line)
                        elementlist.append(line.split()[0])
                if 'CARTESIAN COORDINATES' in line:
                    ocartgrab=True; CartFlag=2
                

if len(elementlist) != numatoms:
    print("Something wrong with grabbed elementlist")
    print("Using dummy atom C as element for all atoms")
    elementlist=numatoms*['C']

#Element list backup
if len(elementlist) < 2:
    print("No element list found. Using dummy atom C as element for all atoms")
    elementlist=numatoms*['C']

#outfilelist=glob.glob('new-cpcm-snapA-78000-A-chargeB-rmin1_0.*out*')
#if len(outfilelist) > 1:
#    print("problem. Many outputfiles in dir matching basename. Remove extra ones")
#    exit()
#else:
#    outfile=outfilelist[0]
#    with open(outfile) as ofile:
        


totparticles=numatoms+numsurfacecharges
#print(molcart)
#print(surfacecart)

outfile=open(sys.argv[1]+".surface.xyz",'w')

outfile.write(str(totparticles)+'\n')
outfile.write("title\n")
for el,at in zip(elementlist,molcart):
    outfile.write(el+' '+str(' '.join(map(str,at)))+'\n')
for ch in surfacecart:
    outfile.write('Li '+str(' '.join(map(str,ch)))+'\n')



outfile.close()
print("Created new XYZ file containing molecule and surface charge positions (labelled as Li atoms) (all coords in Angstrom):\n", sys.argv[1]+".surface.xyz")
