#!/bin/bash

FILESET=$1
DATASET=$2
BOOK=$3
CATALOG=$4
ISDATA=$5
USEGEN=$6
FSRMODE=$7
NEVENTS=$8
SKIPHLTFAIL=$9

OUTPUTDIR=/data/blue/ksung/EWKAna/bacon/$DATASET
SCRAMDIR=/home/ksung/releases/CMSSW_5_2_3/src
RUNDIR=$SCRAMDIR/EWKAna/Ntupler/macros
WORKDIR=`pwd`
export KRB5CCNAME=/home/ksung/.krb5/ticket

cp $RUNDIR/runNtupler.C ./
cp $RUNDIR/rootlogon.C ./
cp $RUNDIR/START50_V15_*_AK5PF.txt ./
cp $RUNDIR/Subdet*BDTG.weights.xml ./

cd $SCRAMDIR
SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`
cd $WORKDIR

echo `hostname`
echo "Running runNtupler.C $FILESET $DATASET $BOOK $CATALOG $ISDATA $USEGEN $FSRMODE $NEVENTS $SKIPHLTFAIL"
echo "Working in $WORKDIR"
ls -l
echo "   Using $CMSSW_BASE"

root -b -l -q ./rootlogon.C runNtupler.C+\(\"${FILESET}\",\"${DATASET}\",\"${BOOK}\",\"${CATALOG}\",${ISDATA},${USEGEN},${FSRMODE},${NEVENTS},${SKIPHLTFAIL}\)

mkdir -p $OUTPUTDIR
mv ${DATASET}_${FILESET}_ntuple.root* $OUTPUTDIR/

rm -vf runNtupler*
rm -vf rootlogon*
rm -vf START50_V15_*_AK5PF.txt
rm -vf Subdet*BDTG.weights.xml

exit 0
