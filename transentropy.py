#!/bin/env python3
import math
#Calculate translational entropy for reaction: C => A +B
#Calculated as in ORCA. Gives about the same values as G09 and results agree with N2H2 => N2 + H2 values in
#Entropy Explained: The Origin of Some Simple Trends J Chem Ed 2002 Paper by Watson and Eisenstein

########################
#INPUT
########################

#Input masses
A_mass=4.002603 #Ddhyd 122389.34
B_mass=2.011       #SH2
C_mass=A_mass+B_mass

#Pressure in atm
P=1
#Temperature
T=298.15

########################

#Constants
#R gas constant in kcal/molK
R=1.987E-3

#Conversion factor for formula. Taken from ORCA code. Confirmed to give same results as G09
factor=0.025607868


#Translation partition function and T*S_trans values for each species
qtrans_A=(factor*T**2.5*A_mass**1.5)/P
print(qtrans_A)
TS_trans_A=T*R*(math.log(qtrans_A)+2.5)
print(TS_trans_A)
qtrans_A=(factor*T**2.5*A_mass**1.5)/P
TS_trans_A=T*R*(math.log(qtrans_A)+2.5)

qtrans_B=(factor*T**2.5*B_mass**1.5)/P
TS_trans_B=T*R*(math.log(qtrans_B)+2.5)

qtrans_C=(factor*T**2.5*C_mass**1.5)/P
TS_trans_C=T*R*(math.log(qtrans_C)+2.5)

#Reaction trans entropy
delta_TS_trans=TS_trans_A+TS_trans_B-TS_trans_C

print("Printing translational entropy for reaction C => A + B   assuming P=", P,  "atm and T=", T, "298.15 K")
print("Using masses: A", A_mass, "B:", B_mass, "C:", C_mass)
print("")
print("TS_trans_A is", TS_trans_A, "kcal/mol")
print("TS_trans_B is", TS_trans_B, "kcal/mol")
print("TS_trans_B is", TS_trans_C, "kcal/mol")
print("")
print("delta_TS_trans is", delta_TS_trans, "kcal/mol")
