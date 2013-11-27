#! /bin/bash
#$ -S /bin/bash

#$ -e /dev/null
#$ -o /dev/null

. lib/env_func.bashrc

usage () {
  echo "qsub -cwd -t 1:NUM_GZ -S /bin/bash `basename $0` -l each_ordering_dir/target_gz.list "
  exit 1
}

while getopts l: OPTION
do
  case $OPTION in
    l)  if [ ! -z "${OPTARG}" ];then GZ_LIST=${OPTARG} ;else usage ;fi
        ;;
    \?) usage ;;
  esac
done

if [ $# -lt 2 ] ; then
  usage
fi

EXE_SPLITN=lib/splitN/sp

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# decompress *copyprobsperlocus.out.gz
#
# sort by position (ascending)
#
# save as *copyprobsperlocus.out files
#
#   in order to make the following sorting much faster by using "-m" 
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

declare -a array

array=($(cat ${GZ_LIST}))

gzfile=""

if [ "${QUEUE_TYPE}" == "SGE" ]; then
  gzfile=${array["SGE_TASK_ID"-1]}
elif [ "${QUEUE_TYPE}" == "LSF" ]; then
  gzfile=${array["LSB_JOBINDEX"-1]}
else
  echo_fail "unknown QUEUE_TYPE: ${QUEUE_TYPE}"
fi

each_copyprobsperlocus=`echo ${gzfile} | perl -pe 's/\.gz//g'`

#
# decompress
# sort (required for msort within each ordering across recipient haplotypes, 
#       maybe not needed if chromopainter is modified in the future)
# and split 
#
date +%Y%m%d_%T

CMD="gzip -dc ${gzfile} | grep -v '^pos' | grep -v '^HAP' | sort -n | ${EXE_SPLITN} -p ${each_copyprobsperlocus}_"
#CMD="gzip -dc ${gzfile} | grep -v '^pos' | grep -v '^HAP' | sort -n | split -50000 - ${each_copyprobsperlocus}_"
echo ${CMD}
eval ${CMD}
if [ $? -ne 0 ]; then 
  echo_fail "Error: ${CMD} "
fi

#
# rm ${gzfile}
#
date +%Y%m%d_%T

CMD="/bin/rm -f $gzfile";
echo ${CMD}
eval ${CMD}
if [ $? -ne 0 ]; then 
  echo_fail "Error: ${CMD} "
fi

