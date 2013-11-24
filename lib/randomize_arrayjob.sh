#! /bin/bash
#$ -S /bin/bash

#$ -e /dev/null
#$ -o /dev/null

. lib/env_func.bashrc

usage () {
  echo "qsub -cwd -t 1:NUM_HAP -S /bin/bash `basename $0` "
  echo " -h file.hap "
  echo " -p outprefix "
  echo " -l strainName.list "
  echo " -o strainName_dispOrder.list "
  echo " -t 1,...,10 (th ordering to be randoized, which is also used to run from the reverse) "
  echo " -s 1 (this value + counter of random ordering => seed of random number generator) "

  exit 1
}

while getopts h:p:l:o:t:s: OPTION
do
  case $OPTION in
    h)  if [ ! -z "${OPTARG}" ];then PHASEFILE=${OPTARG} ;else usage ;fi
        ;;
    p)  if [ ! -z "${OPTARG}" ];then OUT_PREFIX_BASE=${OPTARG} ;else usage ;fi
        ;;
    l)  if [ ! -z "${OPTARG}" ];then HAP_LIST=${OPTARG} ;else usage ;fi
        ;;
    o)  if [ ! -z "${OPTARG}" ];then HAP_LIST_OUTDISP=${OPTARG} ;else usage ;fi
        ;;
    t)  if [ ! -z "${OPTARG}" ];then i_forward_reverse=${OPTARG} ;else usage ;fi
        ;;
    s)  if [ ! -z "${OPTARG}" ];then SEED=${OPTARG} ;else usage ;fi
        ;;
    \?) usage ;;
  esac
done

if [ $# -lt 2 ] ; then
  usage
fi

i_recipient=""
if [ "${QUEUE_TYPE}" == "SGE" ]; then
  i_recipient=${SGE_TASK_ID}
elif [ "${QUEUE_TYPE}" == "LSF" ]; then
  i_recipient=${LSB_JOBINDEX}
else
  echo_fail "unknonw QUEUE_TYPE: ${QUEUE_TYPE}"
fi
let i_recipient=${i_recipient}-1 # 1-indexed => 0-indexed

#
# execute
#
EXE_RANDOMIZE=lib/randomize/rd # already checked by the parent shell

CMD=""
CMD=${CMD}"${EXE_RANDOMIZE} "
CMD=${CMD}" -h ${PHASEFILE}"
CMD=${CMD}" -p ${OUT_PREFIX_BASE}"
CMD=${CMD}" -l ${HAP_LIST}"
CMD=${CMD}" -o ${HAP_LIST_OUTDISP}"
CMD=${CMD}" -t ${i_forward_reverse}"
CMD=${CMD}" -i ${i_recipient}"
CMD=${CMD}" -s ${SEED}" 

echo ${CMD}
eval ${CMD}
