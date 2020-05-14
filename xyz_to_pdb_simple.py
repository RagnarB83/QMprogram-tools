import sys

#Very simple XYZ to PDB converter. No proper residue information information but in correct PDB Format.

try:
    filename=sys.argv[1]
except:
    print("Run:  python3 xyz_to_pdb_simple.pdb filename.xyz")

#Is variable an integer
def isint(s):
    try:
        int(s)
        return True
    except ValueError:
        return False
    
#Read XYZ file
def read_xyzfile(filename):
    #Will accept atom-numbers as well as symbols
    elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
            'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
            'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs',
            'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
            'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa',
            'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
    print("Reading coordinates from XYZfile {} ".format(filename))
    coords=[]
    elems=[]
    with open(filename) as f:
        for count,line in enumerate(f):
            if count == 0:
                if len(line.replace(' ','')) <10:
                    numatoms=int(line.split()[0])
                else:
                    print("Make sure XYZ file is in XMOL format with numatoms as 1st line and title as 2nd line")
                    print("Exiting...")
                    exit()
            if count > 1:
                if len(line.strip()) > 0:
                    if isint(line.split()[0]) is True:
                        elems.append(elements[int(line.split()[0])-1])
                    else:
                        elems.append(line.split()[0])
                    coords.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
    assert len(coords) == numatoms, "Number of coordinates does not match header line"
    assert len(coords) == len(elems), "Number of coordinates does not match elements."
    return elems,coords

#Write PDBfile (dummy version) for PyFrame
def write_pdbfile_dummy(elems,coords,name, atomlabels,residlabels):

    occ_column=1.00
    with open(name+'.pdb', 'w') as pfile:
        resnames = ['LIG'] * len(elems)
        #resnames=['QM', 'QM', 'QM', 'QM', 'QM', 'QM', 'QM', 'QM', 'QM', 'QM', 'QM', 'QM', 'QM', 'HOH', 'HOH','HOH']
        #resids=[1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2]
        #Example:
        #pfile.write("ATOM      1  N   SER A   2      65.342  32.035  32.324  1.00  0.00           N\n")
        for count,(el,c,resname,resid) in enumerate(zip(elems,coords, resnames, residlabels)):
            #print(count, el,c,resname)
            #Dummy resid for everything
            #resid=1
            #Using string format from: https://cupnet.net/pdb-format/
            line="{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(
                'ATOM', count+1, el, '', resname, '', resid, '',    c[0], c[1], c[2], 1.0, occ_column, el, '')
            pfile.write(line+'\n')
    print("Wrote PDB file:", name+'.pdb')


basename = filename.split('.')[0]

elems,coords=read_xyzfile(filename)
residlabels=[1 for i in elems]
write_pdbfile_dummy(elems,coords,basename,elems,residlabels)
