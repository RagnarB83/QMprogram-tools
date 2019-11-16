#!/bin/bash

echo " Script to add an ORCA block to all ORCA inputfiles in dir (assumes .inp extension) after the ! line."
echo " usage: addblocktoinput <read in lines>"
echo ""
echo "Type or paste in the lines (e.g. ORCA block) you want added."
echo "Hit Return key twice when done:"
unset tmp
while :
do
 read line
 [[ $line == "" ]] && tmp="${tmp:0:$((${#tmp}-1))}" && break
 #tmp="$tmp"$line$'\n'
 echo $line >> temporca54blockfile
done

for i in $(ls *.inp )
do
sed -i "/!/r temporca54blockfile" $i
echo "Added block to inputfile $i"
done

rm temporca54blockfile
