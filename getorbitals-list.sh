#! /bin/bash
#
# Originally by J. Song.
# Modified by R. Bjornsson and converted to bash
# Orbital list version
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
var="$@"
orbs=$(echo "${var#* }");orbs=$(echo "${orbs#* }")
echo "Orbitals : $orbs"
if [[ "$1" == "" ]]
then
   echo " usage: getorbitals-list.sh gbwfile [alpha | beta | both] n1  n2  n3 ..."
   echo
   exit
fi

if [[ "$spin" == "alpha" ]]
then
echo "Doing: alpha orbitals..."
elif [[ "$spin" == "beta" ]]
then
echo "Doing: beta orbitals..."
elif [[ "$spin" == "both" ]]
then
echo "Doing: both alpha and beta orbitals..."
else
echo "${red}Second argument should be alpha, beta or both. Example: getorbitals-list.sh file.gbw alpha 5 9 ${normal}"
exit
fi


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

for orb in $orbs
do
echo "Doing orbital:  $orb"
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
$orcadir/orca_plot $INPUT_FILE -i < $SPLOTFILE > $SPLOTFILE.log
rm $SPLOTFILE
echo
echo "Finished plotting."

