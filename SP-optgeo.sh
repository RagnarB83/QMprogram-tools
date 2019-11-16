#!/bin/ksh
# SP-optgeo. Ragnar Bjornsson 2013
# Script reads ORCA geometry optimization jobs, determines if they converged and creates new single-point inputfiles from optimized geometries.
# Script also checks if job is single-point job or surface scan and ignores if so.
echo
if [ "$1" == "" ]
then
echo " usage: SP-optgeo dir"
echo 'dir is directory with outputfiles. Do "SP-optgeo ."  for current dir or "SP-optgeo .." for parent dir'
exit 1
#echo "Reading ORCA output from current directory"
else
filereaddir=$1
echo "Reading ORCA output from: $filereaddir"
fi

ls $filereaddir | grep .out > blax457
while read line
do
	sline=`echo $line | cut -d'.' --complement -f2-`
	hurrtest=`grep 'HURRAY' $filereaddir/$line`
	scantest=`grep 'Relaxed Surface Scan' $filereaddir/$line`
        sptest=`grep 'Single Point Calculation' $filereaddir/$line`


if [[ $sptest == *Single* ]]
then
echo "$line : single-point calculation"
:
else
if [[ $scantest == *Relaxed* ]]
then
echo "$line : surface scan"
:
else
	if [[ $hurrtest == *HURRAY* ]]
	then
	echo "$line : ----------FINISHED---------"
	cp $filereaddir/$sline.inp $sline-XXX.inp
	sed -i '/xyz/,/*/d' $sline-XXX.inp
	sed -i 's/Opt/ /I' $sline-XXX.inp

	numatom=`grep -m 1 'Number of atoms                         ....' $filereaddir/$line | awk '{print $5}'`
	natom5=$(( numatom + 5 ))
	charge=`grep -m 1 'Total Charge           Charge          ....' $filereaddir/$line | awk '{print $5}'`
	mult=`grep -m 1 ' Multiplicity           Mult            ....' $filereaddir/$line | awk '{print $4}'`

        echo "# Optimized geometry from $line " >> $sline-XXX.inp
	echo "*xyz $charge $mult" >> $sline-XXX.inp
	grep -A$natom5 'FINAL ENERGY EVALUATION' $filereaddir/$line | tail -n +7 >blax4k05
	cat $sline-XXX.inp blax4k05 > $sline-optgeo-SP.inp
	echo "*" >> $sline-optgeo-SP.inp

	rm $sline-XXX.inp
	else
	echo "$line  : optimization not converged"
	fi
fi
fi
done <"blax457"
rm blax4k05 blax457

echo
echo "New ORCA input files:"
ls -1 *optgeo-SP.inp


