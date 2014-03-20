#!/bin/sh

if [ $# == 19 ] 
then
## Embedded Grid Engine 'qsub' arguments to simplify things
#$ -S /bin/bash
# Execute from the current working directory
#$ -cwd
# export current user environment
#$ -V
#make sure machines are chosen exclusively
#$ -l excl=true
#specify the cluster nodes to be used
#$ -q qttg.hp@qttg1,qttg.hp@qttg2,qttg.hp@qttg3,qttg.hp@qttg4,qttg.hp@qttg5,qttg.hp@qttg6,qttg.hp@qttg7,qttg.hp@qttg8,qttg.hp@qttg9,qttg.hp@qttg10,qttg.hp@qttg11,qttg.hp@qttg12,qttg.hp@qttg13,qttg.hp@qttg14,qttg.hp@qttg15,qttg.hp@qttg16,qttg.hp@qttg17,qttg.hp@qttg18,qttg.hp@qttg19,qttg.hp@qttg20

   PROGRAM="nonEquilibriumGreensFunction.exe $*"
   RESULTSDIR=${19}
   Nprocs=${1}

   if [ -d $RESULTSDIR ]
      then
         echo "$RESULTSDIR already exists"
      else
         mkdir -m 755 $RESULTSDIR
   fi
   touch $RESULTSDIR/{messagelog.txt,scalars.txt,Transmissions.txt,Greens.txt,ERCD.txt,ORCD.txt,Poissons.txt,SigExp.txt}
   cat /dev/null > $RESULTSDIR/messagelog.txt
   cat /dev/null > $RESULTSDIR/scalars.txt
   cat /dev/null > $RESULTSDIR/Transmissions.txt
   cat /dev/null > $RESULTSDIR/Greens.txt
   cat /dev/null > $RESULTSDIR/ERCD.txt
   cat /dev/null > $RESULTSDIR/ORCD.txt
   cat /dev/null > $RESULTSDIR/Poissons.txt
   cat /dev/null > $RESULTSDIR/SigExp.txt

   date >& $RESULTSDIR/messagelog.txt

   if [ $Nprocs == 1 ]
   then
     ./$PROGRAM >> $RESULTSDIR/messagelog.txt 2>&1
   else
     mpirun --mca mpi_paffinity_alone 1 -np $Nprocs $PROGRAM >> $RESULTSDIR/messagelog.txt 2>&1
   fi
   date >> $RESULTSDIR/messagelog.txt 2>&1

   rm *.sh.{e,o,p}*
   rm core.*
else
   echo "Not enough input parameters: only $# given!">>errormessage.log
fi

exit 0 
