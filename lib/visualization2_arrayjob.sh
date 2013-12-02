#! /bin/bash
#$ -S /bin/bash

. lib/env_func.bashrc

#$ -e /dev/null
#$ -o /dev/null

usage () {
  echo "/bin/bash $0 "
  echo " -a min "
  echo " -b max" 
  echo " -l target_pos_matrix.list" 
  exit 1
}

usage () {
  echo "$0   "
  exit 1
}

echo_fail(){
  echo "`date '+%Y%m%d_%H%M%S'` $*"
  exit 1
}

VISUAL_LIST=""
while getopts a:b:l: OPTION
do
  case $OPTION in
    a)  if [ ! -z "${OPTARG}" ];then MIN=${OPTARG} ;else usage ;fi
        ;;
    b)  if [ ! -z "${OPTARG}" ];then MAX=${OPTARG} ;else usage ;fi
        ;;
    l)  if [ ! -z "${OPTARG}" ];then VISUAL_LIST=${OPTARG} ;else usage ;fi
        ;;
    \?) usage ;;
  esac
done

if [ $# -lt 6 ] ; then
  usage
fi

#
# env
#
LIB_DIR=lib # created by setup.sh

R_MAIN2=./${LIB_DIR}/visualization2.R
R_LIB=./${LIB_DIR}/plotHeatmap.R

#
# visualization of each pos by using arrayjob
#

declare -a array

array=($(cat ${VISUAL_LIST}))

POS_MATRIXFILE=""

if [ "${QUEUE_TYPE}" == "SGE" ]; then
  POS_MATRIXFILE=${array["SGE_TASK_ID"-1]}
elif [ "${QUEUE_TYPE}" == "LSF" ]; then
  POS_MATRIXFILE=${array["LSB_JOBINDEX"-1]}
else
  echo_fail "unknown QUEUE_TYPE: ${QUEUE_TYPE}"
fi

CMD="R --vanilla --quiet < ${R_MAIN2} --args ${R_LIB} ${POS_MATRIXFILE} ${MIN} ${MAX} "
if [ "${VISUAL_LIST}" != "" ]; then
  CMD=${CMD}" ${VISUAL_LIST} "
fi
CMD=${CMD}" > /dev/null 2>&1"
echo ${CMD}
eval ${CMD}


