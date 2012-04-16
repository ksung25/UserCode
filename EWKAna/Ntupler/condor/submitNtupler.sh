#!/bin/bash
#===================================================================================================
# Submit a set of jobs to run over a given dataset, splitting the jobs according to the filesets.
#
#===================================================================================================

# Read the arguments
echo "Starting data processing with arguments: $*"
DATASET=$1
BOOK=$2
CATALOG=$3
ISDATA=$4
USEGEN=$5
FSRMODE=$6
NEVENTS=$7
SKIPHLTFAIL=$8

datestring=`date +-%g%m%d-%H%M%S`

DATADIR=`tail -1  ${CATALOG}/${BOOK}/${DATASET}/Filesets | cut -d' ' -f2`

# Prepare environment
echo " "
echo "Process dataset: $DATASET  of book: $BOOK"
echo "     in catalog: $CATALOG"
echo " "
WORKDIR=`pwd`

# Create the directory for the results
OUTPUTDIR=${WORKDIR}/condor${datestring}
mkdir -p ${OUTPUTDIR}/res

SCRIPT=$WORKDIR/runNtupler.sh

# Looping through each fileset and submitting the condor jobs
LIST=${CATALOG}/${BOOK}/${DATASET}/Filesets

for FILESET in `cat $LIST | cut -d' ' -f1 `
do
  condor=${OUTPUTDIR}/submit_${DATASET}_${FILESET}${datestring}.condor
  
cat > $condor <<EOF
Universe                = vanilla
Notification            = Error
Executable              = $SCRIPT
Arguments               = $FILESET $DATASET $BOOK $CATALOG $ISDATA $USEGEN $FSRMODE $NEVENTS $SKIPHLTFAIL
Input                   = /dev/null
Initialdir              = $OUTPUTDIR
Output                  = ${OUTPUTDIR}/res/${DATASET}_${FILESET}.out
Error                   = ${OUTPUTDIR}/res/${DATASET}_${FILESET}.err
Log                     = /tmp/condor${datestring}_${FILESET}.log
GetEnv                  = True
Rank                    = Mips
Requirements            = (Arch == "X86_64") && (OpSys == "LINUX") && (Disk >= DiskUsage) && ((Memory * 1024) >= ImageSize) && (HasFileTransfer) && (HAS_SSSE3)
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
+AccountingGroup = "research.ksung"
Queue
EOF
  
  # Finally, submit the jobs
  echo "Submitting job for ${DATASET} ${FILESET}..."  
  condor_submit $condor

done

exit 0
