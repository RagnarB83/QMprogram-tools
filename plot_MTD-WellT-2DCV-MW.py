#!/bin/env python3

import os
import sys
import matplotlib.pyplot as plt
import re
import numpy as np
import glob

def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

#Script settings:

#Run Plumed script to get fes.dat from HILLS
#Setting PATH and LD_LIBRARY_PATH for PLUMED. LD-lib path to C library may also be required
os.environ['PATH'] = '/home/bjornsson/plumed-install/bin:/usr/bin:$PATH'
os.environ['LD_LIBRARY_PATH'] = '/home/bjornsson/plumed-install/lib:/opt/gcc-4.9.1/lib64:$LD_LIBRARY_PATH'

# If WellTempered MetaDynamics (generally recommended). For regular MTD a pointless plot would be generated for Welltemp=True
WellTemp=True

#Note: Script currently assumes the CV to be a torsion in radians. Gets converted to degrees.
#Note: TODO: Support other CVs
#PLUMED uses kJ/mol be default. This is assumed but is converted to kcal/mol here.

#################################
#END OF USER-REQUIRED SETTINGS
################################
print("Metadynamics Analysis Script by Ragnar Bjornsson")
print("")

try:
    path=sys.argv[1]
    os.chdir(path)
    print("Changing dir to path: ", path)
except:
    print("Assuming current dir contains COLVAR/HILLS files")


#Checking if Multiple Walker MTD or not (by analyzing whether HILLS.X files are present or not)
try:
    f = open("HILLS.0")
    f.close()
    print("Found numbered HILLS.X file. This is a multiple-walker run.")
    MultipleWalker=True
except IOError:
    try:
        f = open("HILLS")
        f.close()
        print("This is single-walker run")
        MultipleWalker=False
    except FileNotFoundError:
        print("Found no HILLS.X or HILLS file. Exiting...")
        exit()


#The plumed sum_hills command that is run.
print("")
if MultipleWalker==True:
    #Removing old HILLS.ALL if present
    try:
        os.remove('HILLS.ALL')
    except:
        pass
    #Gathering HILLS files
    #HILLSFILELIST=sorted(glob.glob("HILLS*"))
    HILLSFILELIST=natural_sort(glob.glob("HILLS*"))
    #Which COLVAR file to look at
    #COLVARFILELIST=sorted(glob.glob("COLVAR*"))
    COLVARFILELIST=natural_sort(glob.glob("COLVAR*"))
    print("MW= True. Concatenating files to HILLS.ALL")
    #os.system('cat HILLS.* > HILLS.ALL')
    print("HILLSFILELIST:", HILLSFILELIST)
    with open('HILLS.ALL', 'w') as outfile:
        for hfile in HILLSFILELIST:
            with open(hfile) as infile:
                for line in infile:
                    outfile.write(line)

    print("Running plumed to sum hills...")
    print("")
    os.system('plumed sum_hills --hills HILLS.ALL')
else:
    os.system('plumed sum_hills --hills HILLS')
    #HILLSFILE="HILLS"
    HILLSFILELIST=['HILLS']
    #Single COLVAR file
    COLVARFILELIST=['COLVAR']

print("")
print("COLVAR files:", COLVARFILELIST)
print("HILLS files:", HILLSFILELIST)
###########################################
# 0 K PES curve for comparison on plot
##########################################
PotCurve=True
if PotCurve==True:
    #Getting 0 Kelvin potential energy curve from file
    #File should be :  X-value: Torsion in Deg  Y-value: Energy in hartree
    potcurve_degs=[]
    potcurve_energy_au=[]
    try:
        with open("potcurve") as potfile:
            for line in potfile:
                if '#' not in line:
                    potcurve_degs.append(float(line.split()[0]))
                    potcurve_energy_au.append(float(line.split()[1]))
        potcurve_energy_kcal=np.array(potcurve_energy_au)*627.509
        potcurve_Relenergy_kcal=potcurve_energy_kcal-min(potcurve_energy_kcal)
    except FileNotFoundError:
        PotCurve=False
        print("File potcurve not found. Add file if desired.")
########################################
pi=3.14159265359
#Get temperature from plumed.in in dir or dir above
dihed1atoms=[]
dihed2atoms=[]
print("")
try:
    with open("plumed.in") as pluminpfile:
        print("Found plumed.in file. Reading variables")
        for line in pluminpfile:
            if '#' not in line:
                if 'TORSION' in line:
                    if len(dihed1atoms) > 0:
                        x=line.split()[-1]
                        y=line.split('=')[-1]
                        for z in y.split(','):
                            dihed2atoms.append(int(z))
                    else:
                        x=line.split()[-1]
                        y=line.split('=')[-1]
                        for z in y.split(','):
                            dihed1atoms.append(int(z))
                if 'TEMP' in line:
                    for x in line.split():
                        if 'TEMP' in x:
                            temperature=float(x.split('=')[1])
    print("Found temperature:", temperature)
except:
    print("Found no plumed.in in dir")
    print("Trying dir above...")
    try:
        with open("../plumed.in") as pluminpfile:
            print("Found plumed.in file. Reading variables")
            for line in pluminpfile:
                if '#' not in line:
                    if 'TORSION' in line:
                        if len(dihed1atoms) > 0:
                            x=line.split()[-1]
                            y=line.split('=')[-1]
                            for z in y.split(','):
                                dihed2atoms.append(int(z))
                        else:
                            x=line.split()[-1]
                            y=line.split('=')[-1]
                            for z in y.split(','):
                                dihed1atoms.append(int(z))
                    if 'TEMP' in line:
                        for x in line.split():
                            if 'TEMP' in x:
                                temperature=float(x.split('=')[1])
        print("Found temperature:", temperature)
    except:
        print("Unknown exception occurred when reading plumed.in")
        print("Setting temp to unknown")
        temperature="Unknown"

print("Atoms in CV1:", dihed1atoms)
if len(dihed2atoms)>0:
    print("Atoms in CV2:", dihed2atoms)
#READ HILLS. Only necessary for Well-Tempered Metadynamics and plotting of Gaussian height
if WellTemp==True:
    time_hills=[]
    gaussheight=[]
    time_hills_list=[]
    gaussheightkcal_list=[]
    for hillsfile in HILLSFILELIST:
        with open(hillsfile) as hillsf:
            for line in hillsf:
                if 'FIELDS' in line:
                    biasfcolnum=int(line.split().index('biasf'))
                if '#' not in line:
                    if biasfcolnum==6:
                        time_hills.append(float(line.split()[0]))
                        gaussheight.append(float(line.split()[3]))
                    if biasfcolnum==8:
                        time_hills.append(float(line.split()[0]))
                        gaussheight.append(float(line.split()[5]))
        gaussheight_kcal=np.array(gaussheight)/4.184
        time_hills_list.append(time_hills)
        gaussheightkcal_list.append(gaussheight_kcal)
        time_hills=[];gaussheight_kcal=[];gaussheight=[]

#READ COLVAR
time=[]
colvar_value=[]
colvar2_value=[]
biaspot_value=[]

colvar_value_deg_list=[]
colvar2_value_deg_list=[]
biaspot_value_kcal_list=[]
time_list=[]

for colvarfile in COLVARFILELIST:
    with open(colvarfile) as colvarf:
        for line in colvarf:
            if 'FIELDS' in line:
                biascolnum = [i for i, s in enumerate(line.split()) if '.bias' in s][0]
            if '#' not in line:
                try:
                    #1 CVs
                    if biascolnum==4:
                        CVnum=1
                        biaspot_value.append(float(line.split()[2]))
                        time.append(float(line.split()[0]))
                        colvar_value.append(float(line.split()[1]))
                    #2 CVs
                    elif biascolnum==5:
                        CVnum=2
                        biaspot_value.append(float(line.split()[3]))
                        time.append(float(line.split()[0]))
                        colvar_value.append(float(line.split()[1]))
                        colvar2_value.append(float(line.split()[2]))
                    else:
                        print("unknown format of COLVAR file. More than 2 CVs ??")
                        exit()
                except:
                    pass
    #convert to deg
    colvar_value_deg=np.array(colvar_value)*180/pi
    colvar2_value_deg=np.array(colvar2_value)*180/pi
    #Convert to kcal
    biaspot_value_kcal=np.array(biaspot_value)/4.184

    #New. For multiple COLVAR files we create lists of colvar_value_deg, colvar2_value_deg and biaspot_value_kcal
    colvar_value_deg_list.append(colvar_value_deg)
    colvar2_value_deg_list.append(colvar2_value_deg)
    biaspot_value_kcal_list.append(biaspot_value_kcal)
    time_list.append(time)
    time=[];biaspot_value_kcal=[];colvar2_value_deg=[];colvar_value_deg=[]
    biaspot_value=[];colvar2_value=[];colvar_value=[]
#READING fes.dat
#Reaction coordinates (radian if torsion)
rc=[]
rc2=[]
#Free energy (kJ/mol)
free_energy=[]

#Derivative of Free Energy vs. reaction-coordinate. Probably not useful
derivG=[]
derivG2=[]

#Reading file
##! FIELDS dihed1 dihed2 file.free der_dihed1 der_dihed2
with open("fes.dat") as fesfile:
    for line in fesfile:
        if '#' not in line and len(line.split()) > 0:
            if CVnum==1:
                rc.append(float(line.split()[0]))
                free_energy.append(float(line.split()[1]))
                derivG.append(float(line.split()[2]))
            else:
                rc.append(float(line.split()[0]))
                rc2.append(float(line.split()[1]))
                free_energy.append(float(line.split()[2]))
                derivG.append(float(line.split()[3]))
                derivG2.append(float(line.split()[4]))
#rc is in rad
#convert to deg
rc_deg=np.array(rc)*180/pi
rc2_deg=np.array(rc2)*180/pi
#Convert free energy from kJ/mol to kcal/mol
free_energy_kcal=np.array(free_energy)/4.184
Relfreeenergy_kcal=free_energy_kcal-min(free_energy_kcal)

###################
# Matplotlib part
###################
print("")
print("Now plotting via Matplotlib")

if CVnum==1:
    print("CVs:", 1)
    #Space between subplots
    plt.subplots_adjust(hspace=0.4)
    #Subplot 1: Free energy surface. From fes.dat via HILLS file (single-walker) or HILLS.X files (multiple-walker)
    plt.subplot(2, 2, 1)
    plt.gca().set_title('Free energy vs. CV')
    plt.xlabel('Torsion (°)')
    plt.ylabel('Energy (kcal/mol)')
    plt.xlim([-180,180])
    #plt.plot(rc_deg, free_energy_kcal, marker='o', linestyle='-', markerwidth is , linewidth=1, label='G (kcal/mol)')
    plt.plot(rc_deg, Relfreeenergy_kcal, marker='o', linestyle='-', linewidth=1, markersize=3, label='G (kcal/mol): {} K'.format(temperature))
    if PotCurve==True:
        plt.plot(potcurve_degs, potcurve_Relenergy_kcal, marker='o', linestyle='-', markersize=3, linewidth=1, label='E (kcal/mol): 0 K', color='orange')
    plt.legend(shadow=True, fontsize='xx-small', loc='upper left')
    #Subplot 2: CV vs. time. From COLVAR file/files.
    plt.subplot(2, 2, 2)
    plt.gca().set_title('CV vs. time')
    plt.xlabel('Time (ps)')
    plt.ylabel('Torsion (°)')
    #New: Using first timelist to get x-axis limit
    plt.xlim([0,max(time_list[0])+5])

    #New. For MW-MTD we have multiple trajectories. Time should be the same
    for num,(t,cv_deg) in enumerate(zip(time_list,colvar_value_deg_list)):
        plt.plot(t, cv_deg, marker='o', linestyle='-', linewidth=0.5, markersize=2, label='Walker'+str(num))
    #lg = plt.legend(shadow=True, fontsize='xx-small', bbox_to_anchor=(1.3, 1.0), loc='upper right')

    #Subplot 3: Bias potential from COLVAR
    plt.subplot(2, 2, 3)
    #plt.title.set_text('Bias potential')
    plt.gca().set_title('Bias potential')
    plt.xlabel('Torsion (°)')
    plt.xlim([-180,180])
    for num,(cv_deg,biaspot) in enumerate(zip(colvar_value_deg_list,biaspot_value_kcal_list)):
        plt.scatter(cv_deg, biaspot, marker='o', linestyle='-', s=3, linewidth=1, label='Walker'+str(num))
    #lg2 = plt.legend(shadow=True, fontsize='xx-small', bbox_to_anchor=(0.0, 0.0), loc='lower left')

    if WellTemp==True:
        #Subplot 4: Gaussian height from HILLS
        plt.subplot(2, 2, 4)
        plt.gca().set_title('G-height vs. time')
        plt.xlabel('Time (ps)')
        plt.xlim([0,max(time_hills_list[0])+5])
        for num,(th,gh) in enumerate(zip(time_hills_list,gaussheightkcal_list)):
            #plt.scatter(th, gh, marker='o', linestyle='-', s=3, linewidth=1, label='G height')
            plt.plot(th, gh, marker='o', linestyle='-', markersize=2, linewidth=0.5, label='Walker'+str(num))
        plt.legend(shadow=True, fontsize='xx-small', loc='lower right', bbox_to_anchor=(1.3, 0.0))

elif CVnum==2:
    print("CVs:", 2)

    def flatten(list):
        return [item for sublist in list for item in sublist]

    #2CV-MW plots will be too messy so combinining walker information
    colvar_value_deg_list_flat=flatten(colvar_value_deg_list)
    colvar2_value_deg_list_flat=flatten(colvar2_value_deg_list)
    biaspot_value_kcal_list_flat=flatten(biaspot_value_kcal_list)
    time_hills_flat=flatten(time_hills)
    time_flat=flatten(time_list)
    gaussheight_kcal_flat=flatten(gaussheight_kcal)

    #Space between subplots
    plt.subplots_adjust(hspace=0.4)
    plt.subplots_adjust(wspace=0.4)

    #Subplot 1: Free energy surface
    plt.subplot(2, 2, 1)
    plt.gca().set_title('Free energy vs. CV')
    plt.xlabel('Dihedral ({})'.format(dihed1atoms))
    plt.ylabel('Dihedral ({})'.format(dihed2atoms))
    plt.xlim([-180,180])
    plt.ylim([-180,180])
    cm = plt.cm.get_cmap('RdYlBu')
    colorscatter=plt.scatter(rc_deg, rc2_deg, c=Relfreeenergy_kcal, marker='o', linestyle='-', linewidth=1, cmap=cm)
    cbar = plt.colorbar(colorscatter)
    cbar.set_label('deltaG (kcal/mol)')

    #Subplot 2: CV vs. time
    plt.subplot(2, 2, 2)
    plt.gca().set_title('CV vs. time')
    plt.xlabel('Dihedral ({})'.format(dihed1atoms))
    plt.ylabel('Dihedral ({})'.format(dihed2atoms))
    #plt.xlim([0,max(time)+5])
    cm = plt.cm.get_cmap('RdYlBu')
    colorscatter=plt.scatter(colvar_value_deg_list_flat, colvar2_value_deg_list_flat, c=time_flat, marker='o', s=2, linestyle='-', linewidth=1, cmap=cm)
    cbar = plt.colorbar(colorscatter)
    cbar.set_label('Time (ps)')

    #Subplot 3: Bias potential
    plt.subplot(2, 2, 3)
    plt.gca().set_title('Bias potential')
    plt.xlabel('Dihedral ({})'.format(dihed1atoms))
    plt.ylabel('Dihedral ({})'.format(dihed2atoms))
    cm = plt.cm.get_cmap('RdYlBu')
    colorscatter2=plt.scatter(colvar_value_deg_list_flat, colvar2_value_deg_list_flat, c=biaspot_value_kcal_list_flat, marker='o', linestyle='-', linewidth=1, cmap=cm)
    cbar2 = plt.colorbar(colorscatter2)
    cbar2.set_label('Biaspot (kcal/mol)')
    #lg = plt.legend(fontsize='xxx-small', bbox_to_anchor=(1.05, 1.0), loc='lower left')

    #Subplot 4: Gaussian height
    plt.subplot(2, 2, 4)
    plt.gca().set_title('G-height vs. time')
    plt.xlabel('Time (ps)')
    plt.ylabel('Gaussian height (kcal/mol)')
    plt.xlim([0,max(time_hills_list[0])+5])
    plt.xlim([0,max(time_hills_list[0])+5])
    for num,(th,gh) in enumerate(zip(time_hills_list,gaussheightkcal_list)):
        plt.plot(th, gh, marker='o', linestyle='-', markersize=2, linewidth=0.5, label='Walker'+str(num))
    #plt.legend(shadow=True, fontsize='xx-small')
    plt.legend(fontsize='xxx-small', bbox_to_anchor=(1.3, 0.0), loc='lower right')

# loc='upper right', bbox_to_anchor=(0.5, 0.5)
#Saving figure
maxtime=int(max(time_list[0]))
plt.savefig("MTD_Plot-"+str(maxtime)+"ps"+".png",
            dpi=300,
            format='png')

 #           bbox_inches='tight')
#bbox_extra_artists = (lg),
plt.show()