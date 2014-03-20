#!/bin/sh

## This script submits the single run shell script to the cluster via the Sun Grid Engine (qsub)
Nprocs=72
bcy=0
bcz=0
channelType=0
contactType=0
mul='-1D0'
mur='0D0'
muT='0D0'
a0='1D0'
Npx=10 #30
Npy=10 #8
Npz=8 #8

randSeed=1
corrLength='0D0'
includePoisson=0
iteMax=20
abPhi='0D0'


for disorderStrength in '0D0' #'0D0' '1D0' '2D0'
do
 for corrLength in '0D0' #'1D-1' '2D-1'
 do
  if [ $contactType -eq 0 ]
  then
    p2="BigSI"
  elif [ $contactType -eq 1 ]
  then
    p2="BigP"
  elif [ $contactType -eq 2 ]
  then
    p2="SmallP"
  fi
  p1="${randSeed}"

  RESULTSDIR="TR${p1}${p2}_Mat${channelType}_a0${a0}_${Npx}x${Npy}x${Npz}_bc${bcy}${bcz}_mu${mul}_${mur}_muT${muT}_Phi${abPhi}_Disorder${disorderStrength}x${corrLength}_Poisson${includePoisson}"
  PROGRAM="Single_Run.sh $Nprocs $Npx $Npy $Npz $bcy $bcz $mul $mur $muT $channelType $contactType $a0 $abPhi $disorderStrength $corrLength $randSeed $includePoisson $iteMax $RESULTSDIR"
  echo "qsub -pe orte $Nprocs $PROGRAM" | sh -x
  sleep 1
 done
done
exit 0

