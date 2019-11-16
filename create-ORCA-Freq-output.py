#!/usr/bin/env python3

#Create fake ORCA-outputfile containing geometry and Freq output from ORCA-style Hessian file and orca_vib output.
# For visualization in Chemcraft

#Read in Hessian file as first argument and outputfile from orca_vib as second

import os
import sys

hessfile=sys.argv[1]
orcaviboutputfile=sys.argv[2]



orca_header="""                                 *****************
                                 * O   R   C   A *
                                 *****************

           --- An Ab Initio, DFT and Semiempirical electronic structure package ---

                  #######################################################
                  #                        -***-                        #
                  #  Department of molecular theory and spectroscopy    #
                  #              Directorship: Frank Neese              #
                  # Max Planck Institute for Chemical Energy Conversion #
                  #                  D-45470 Muelheim/Ruhr              #
                  #                       Germany                       #
                  #                                                     #
                  #                  All rights reserved                #
                  #                        -***-                        #
                  #######################################################


                         Program Version 3.0.3 - RELEASE   -




                       *****************************
                       * Geometry Optimization Run *
                       *****************************

         *************************************************************
         *                GEOMETRY OPTIMIZATION CYCLE   1            *
         *************************************************************
---------------------------------
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------"""

#Grab coordinates from hessfile
numatomgrab=False
cartgrab=False
elements=[]
x_coord=[]
y_coord=[]
z_coord=[]
count=0
bohrang=0.529177249
with open(hessfile) as hfile:
    for line in hfile:
        if cartgrab==True:
            count=count+1
            #print(line)
            #print(count)
            elem=line.split()[0]; x_c=bohrang*float(line.split()[2]);y_c=bohrang*float(line.split()[3]);z_c=bohrang*float(line.split()[4])
            elements.append(elem)
            x_coord.append(x_c);y_coord.append(y_c);z_coord.append(z_c)
            if count == numatoms:
                break
        if numatomgrab==True:
            numatoms=int(line.split()[0])
            numatomgrab=False
            cartgrab=True
        if "$atoms" in line:
            numatomgrab=True
            
#Grab Normal mode output from vibout file

vibout=[]
grabstuff=False
with open(orcaviboutputfile) as vfile:
    for line in vfile:
        if grabstuff==True:
            if '--------------------------------------' in line:
                grabstuff=False; break
            else:
                vibout.append(line)
        if 'VIBRATIONAL FREQUENCIES' in line:
            vibout.append(line)
            grabstuff=True


filename='bla'
outfile = open(filename+'fake.out', 'w')

outfile.write(orca_header+'\n')
for el,x,y,z in zip(elements,x_coord,y_coord,z_coord):
     line = "  {0:2s} {1:11.6f} {2:12.6f} {3:13.6f}".format(el, x, y, z)
     outfile.write(line+'\n')


outfile.write('\n')
outfile.write('-----------------------\n')
for l in vibout:
    outfile.write(l)

outfile.close()

