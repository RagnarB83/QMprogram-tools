#! /bin/bash
#
# Originally by J. Song.
# Modified and converted to bash by R. Bjornsson
# Orbital-range version
normal=`tput sgr0`
bold=`tput bold`
red=`tput setaf 1`
yellow=`tput setaf 3`
green=`tput setaf 2`
cyan=`tput setaf 6`
magenta=`tput setaf 5`
underline=`tput smul`

orcadir=/opt/orca_4.2.0
echo "${cyan}Using orca_plot present in $orcadir${normal}"
INPUT_FILE=$1
spin=$2

if [[ "$1" == "" ]]
then
  echo "${green}usage: getorbitals-range.sh gbwfile [alpha | beta | both] range-begin range-end ..."
  echo "Example: getorbitals-range.sh file.gbw alpha 5 9    . Will print alpha orbitals 5-9 ${normal}"
  echo
  exit
fi

if [[ "$3" == "" ]]
then
echo "${red}Please give orbital range: e.g. getorbitals-range.sh file.gbw alpha 5 9 ${normal}"
exit
else
range1=$3
fi

if [[ "$4" == "" ]]
then
range2=$range1
else
range2=$4
fi

if [[ "$spin" == "alpha" ]]
then
echo "Printing alpha orbitals..."
elif [[ "$spin" == "beta" ]]
then
echo "Printing beta orbitals..."
elif [[ "$spin" == "both" ]]
then
echo "Printing both alpha and beta orbitals..."
else
echo "${red}Second argument should be alpha, beta or both. Example: getorbitals-range.sh file.gbw alpha 5 9 ${normal}"
exit
fi

echo "Orbital range: $range1  to $range2 "
SPLOTFILE="SPLOT_FILE"
if [[ -e $SPLOTFILE ]]
then
 rm $SPLOTFILE
fi

# grid
echo 4 >> $SPLOTFILE
echo 80 >> $SPLOTFILE
# Gaussian cube file print
echo 5 >> $SPLOTFILE
echo 7 >> $SPLOTFILE
#

for orb in $(seq $range1 $range2)
do
   if [[ $spin == alpha ]]
   then
       echo 3 >> $SPLOTFILE
       echo 0 >> $SPLOTFILE
       echo 2 >> $SPLOTFILE
       echo $orb >> $SPLOTFILE
       echo 10 >> $SPLOTFILE
   fi
   if [[ $spin == beta ]]
   then
       echo 3 >> $SPLOTFILE
       echo 1 >> $SPLOTFILE
       echo 2 >> $SPLOTFILE
       echo $orb >> $SPLOTFILE
       echo 10 >> $SPLOTFILE
   fi
   if [[ $spin == both ]]
   then
       echo 3 >> $SPLOTFILE
       echo 0 >> $SPLOTFILE
       echo 2 >> $SPLOTFILE
       echo $orb >> $SPLOTFILE
       echo 10 >> $SPLOTFILE

       echo 3 >> $SPLOTFILE
       echo 1 >> $SPLOTFILE
       echo 2 >> $SPLOTFILE
       echo $orb >> $SPLOTFILE
       echo 10 >> $SPLOTFILE

   fi
done

echo 11 >> $SPLOTFILE

#run
$orcadir/orca_plot $INPUT_FILE -i < $SPLOTFILE >$SPLOTFILE.log
rm $SPLOTFILE
echo
echo "Finished plotting."
