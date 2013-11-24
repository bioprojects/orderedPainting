#! /bin/bash
#$ -S /bin/bash

##$ -e /dev/null
##$ -o /dev/null

. lib/env_func.bashrc

usage () {
  echo "qsub -cwd -t 1:NUM_HAP -S /bin/bash `basename $0` -r prefix_linked.rec -n prefix_linked.neest -l prefix_ordered_hap.list "
  exit 1
}

while getopts r:n:l: OPTION
do
  case $OPTION in
    r)  if [ ! -z "${OPTARG}" ];then RECOMB_FILE=${OPTARG} ;else usage ;fi
        ;;
    n)  if [ ! -z "${OPTARG}" ];then NeFILE=${OPTARG} ;else usage ;fi
        ;;
    l)  if [ ! -z "${OPTARG}" ];then HAP_LIST=${OPTARG} ;else usage ;fi
        ;;
    \?) usage ;;
  esac
done

if [ $# -lt 2 ] ; then
  usage
fi

EXE_PAINT=./lib/chromopainter

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

globalparams=`cat $NeFILE` # $NeFILE contains the commands to tell ChromoPainter about both Ne and global mutation rate, e.g. "-n 10000 -M 0.01".

#
# chromopainter linkage model 
#   -a 1 1 (recipient is always the 1st haplotype)
#
#   by specifying a recombination map file and globalparams
#

declare -a array

array=($(cat ${HAP_LIST}))

phasefile=""
if [ "${QUEUE_TYPE}" == "SGE" ]; then
  phasefile=${array["SGE_TASK_ID"-1]}
elif [ "${QUEUE_TYPE}" == "LSF" ]; then
  phasefile=${array["LSB_JOBINDEX"-1]}
else
  echo_fail "unknown QUEUE_TYPE: ${QUEUE_TYPE}"
fi

OUT_PREFIX=`echo ${phasefile} | perl -pe 's/\.hap$//g'`

IND_RECIP=`head -2 ${phasefile} | tail -1`

CMD="${EXE_PAINT} $globalparams -a ${IND_RECIP} ${IND_RECIP} -j -g $phasefile -r $RECOMB_FILE -o ${OUT_PREFIX}"
if [ "${SELF_COPY}" == "TRUE" ]; then
  CMD="${CMD} -c"
fi
# -b: print-out zipped file  with suffix '.copyprobsperlocus.out' 
CMD="${CMD} -b"
if [ "${VERBOSE}" == "FALSE" ]; then
  CMD=${CMD}" > /dev/null"
fi
echo ${CMD}
eval ${CMD}

if [ $? -ne 0 ]; then 
  echo_fail "Error: ${CMD} "
else
  echo "$phasefile was painted, so it is removed here."
  # remove the input .hap file created by step2 (to save disk space)
  if [ -f "$phasefile" ]; then
    CMD="/bin/rm -f $phasefile"
    echo ${CMD}
    eval ${CMD}
  fi
fi

