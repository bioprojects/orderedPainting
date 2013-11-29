#! /bin/bash
. lib/env_func.bashrc

usage () {

  if [ "${QUEUE_TYPE}" == "SGE" ]; then
    echo "qsub -cwd  <<< '"
  elif [ "${QUEUE_TYPE}" == "LSF" ]; then
    echo "bsub -o logfile.txt <<< '"
  else
    echo_fail "unknown QUEUE_TYPE: ${QUEUE_TYPE}"
  fi

  echo "/bin/bash $0 "
  echo "  -g file.hap "
  echo "  -l strainName.list (individual name in the hap file)" 
  echo " [-n 20 (num. of dirs where uncompressed tmp files are processed simultaneously: default=20) ]"
  echo " [-m pos2missingInd.txt (pos[tab]missing_individual_name]"
  echo " [-o strainName.list (individual name in an order for output: default is the file specified by -l) ]"
  
  echo " [-t 10 (num. of orderings and the reverse, default=10) ]"
  echo " [-s 1  (seed of random number generator: default=1) ]"

  echo "' "
  exit 1
}

################################################################################################################

#
# env
#
LIB_DIR=lib # created by setup.sh
LOG_DIR=log # created by setup.sh

EXE_PAINT=./${LIB_DIR}/chromopainter # also used in SH_PAINT_QSUB

PL_CHECK_LINE_LEN=./${LIB_DIR}/check_lineLen.pl
PL_MAKE_RECMAP=./${LIB_DIR}/makeuniformrecfile.pl
PL_ESTIMATE_Ne=./${LIB_DIR}/neaverage.pl

EXE_RANDOMIZE=./${LIB_DIR}/randomize/rd
SH_RANDOMIZE=./${LIB_DIR}/randomize_arrayjob.sh
SH_PAINT_QSUB=./${LIB_DIR}/chromopainter_linkage_orderings_arrayjob.sh
PL_SITE_BY_SITE=./${LIB_DIR}/create_examine_site_by_site_matrices.pl
SH_DECOMPRESS_SORT_SPLIT_EACH_ORDERING=./${LIB_DIR}/decompress_sort_split_gz_arrayjob.sh
PL_MSORT_EACH_ORDERING=./${LIB_DIR}/msort_each_ordering.pl
EXE_SORT=./${LIB_DIR}/sort 
EXE_SPLITN=./${LIB_DIR}/splitN/sp
EXE_PP=./${LIB_DIR}/postprocess/pp

R_MAIN1=./${LIB_DIR}/visualization1.R
R_LIB_HEATMAP=./${LIB_DIR}/plotHeatmap.R

R_MAIN2=./${LIB_DIR}/visualization2.R
SH_R_MAIN2=./${LIB_DIR}/visualization2_arrayjob.sh

R_CHECK_MISSING_STAT=./${LIB_DIR}/check_missingCnt_distStat.R


# to check existence below
arr_executable_files=(
  $EXE_PAINT
  $PL_CHECK_LINE_LEN
  $PL_MAKE_RECMAP
  $PL_ESTIMATE_Ne
  $EXE_RANDOMIZE
  $SH_RANDOMIZE
  $EXE_PAINT
  $PL_SITE_BY_SITE
  $SH_DECOMPRESS_SORT_SPLIT_EACH_ORDERING
  $PL_MSORT_EACH_ORDERING
  $EXE_SORT
  $EXE_SPLITN
  $EXE_PP
)

arr_R_files=(
  $R_MAIN1
  $R_MAIN2
  $R_LIB_HEATMAP
  $R_CHECK_MISSING_STAT
)

#
# constant
#
NUM_EM=30 # 0-indexed. 10 is sometimes not enough for convergence
NUM_SPLITN=9

#
# rule
#
GZ_CAT_COPYPROB_EACH_DIR=copyprobsperlocus.cat.gz

#
# vars
#
OUT_DIR=""


################################################################################################################

#
# func 
#   used only in this script
#   other functions are defined in lib/env_func.bashrc
#

move_log_files() {
  declare -a arr_logs
  arr_logs=`find . -maxdepth 1 -name $1\*$2\* -print`
  #if ls $1*$2* &> /dev/null; then
  for each_log in ${arr_logs[@]}
  do
    if [ -f "${each_log}" ]; then
      /bin/mv ${each_log} ${LOG_DIR}/ 
      echo "log file ${each_log} was moved to ${LOG_DIR}/"
    fi
  done
}

returnQSUB_CMD() {
  QSUB_CMD=""
  if [ "${QUEUE_TYPE}" == "SGE" -o "${QUEUE_TYPE}" == "UGE" ]; then
    QSUB_CMD="${QSUB_COMMON} -o $1.log -e $1.log -N $1"
  elif [ "${QUEUE_TYPE}" == "LSF" ]; then
    QSUB_CMD="${QSUB_COMMON} -o $1.log -e $1.log -J $1"
  fi
  
  if test "$2" = "" ; then
    QSUB_CMD=${QSUB_CMD}
  else
    if [ "${QUEUE_TYPE}" == "SGE" -o "${QUEUE_TYPE}" == "UGE" ]; then
      QSUB_CMD="${QSUB_CMD} -t $2:$3"
    elif [ "${QUEUE_TYPE}" == "LSF" ]; then
      QSUB_CMD="${QSUB_CMD}[$2-$3]"
    fi
  fi
  echo ${QSUB_CMD}
}

submit_calcAveDist_ordering() {
  CMD=`returnQSUB_CMD ${STAMP}`
  CMD=${CMD}" <<< '"
  CMD=${CMD}" perl ${PL_SITE_BY_SITE}"
  CMD=${CMD}" -g ${PHASEFILE} "
  CMD=${CMD}" -d ${ORDER_DIR_LIST} "
  CMD=${CMD}" -l ${ORDER_STRAIN_LIST} "
  #if [ "${CONTRAST_MAX}" -gt 0 ]; then
  #  CMD=${CMD}" -c ${CONTRAST_MAX} " 
  #fi
  CMD=${CMD}" -s ${HAP_LIST_OUTDISP} "
  #CMD=${CMD}" -n ${i_ordering}"
  CMD=${CMD}" -n $1"
  if [ "${MISSING_POS_IND_FILE}" != "" ]; then
    CMD=${CMD}" -m ${MISSING_POS_IND_FILE}"
  fi
  if [ "${CONSTRAINT_FILE}" != "" ]; then
    CMD=${CMD}" -c ${CONSTRAINT_FILE}"
  fi
  CMD=${CMD}"'"

  echo ${CMD}
  eval ${CMD}
  if [ $? -ne 0 ]; then 
    echo_fail "Execution error: ${CMD} (step${STEP})"
  fi
}

wait_until_finish() {
  while :
  do
    QSTAT_CMD=qstat
    END_CHECK=`${QSTAT_CMD} | grep -e $1 -e 'not ready yet' | wc -l`
    if [ "${END_CHECK}" -eq 0 ]; then
      break
    fi
    sleep 10
  done
}

submit_msort_each_ordering() {
  i_dir=0
  while [ "$i_dir" -lt ${#arr_dirs_for_msort[@]} ]; do
    date +%Y%m%d_%T

    ARRAY_S=1
    ARRAY_E=9

    CMD=`returnQSUB_CMD ${STAMP} ${ARRAY_S} ${ARRAY_E}`
    CMD=${CMD}" ${PL_MSORT_EACH_ORDERING} -d ${arr_dirs_for_msort[$i_dir]} -g ${PHASEFILE} "
    echo ${CMD}
    eval ${CMD}
    if [ $? -ne 0 ]; then 
      echo_fail "Execution error: ${CMD} (step${STEP}) "
    fi
    
    #
    # update the array by removing the submitted dir
    #
    unset arr_dirs_for_msort[$i_dir]
    arr_dirs_for_msort=(${arr_dirs_for_msort[@]})
  done
  # this loop automatically ends after submitting msort jobs for the decompressed dirs
}

get_stamp() {
  DATE_N=`date +%N`
  TIME_DIGIT=`echo "$(date +%H%S)$(printf '%02d' $(expr $DATE_N / 10000000))"` # hour(2)sec(2)milisec(2)
  #if ls "s$1_${TIME_DIGIT}"* &> /dev/null; then
  while ls "s$1_${TIME_DIGIT}"* &> /dev/null
  do
    sleep 1
    DATE_N=`date +%N`
    TIME_DIGIT=`echo "$(date +%H%S)$(printf '%02d' $(expr $DATE_N / 10000000))"` # hour(2)sec(2)milisec(2)
  done
  #fi
  STAMP="s$1_${TIME_DIGIT}"
}

disp_punctuate() {
  echo "*************** STEP$1, log=$2.log ****************************** "
  date +%Y%m%d_%T
}

################################################################################################################

#
# args
#
SEED=1
TYPE_NUM_ORDERING=10
VERBOSE=FALSE
CONTRAST_MAX=-9999
MAX_PARALLEL_DECOMPRESS=20 

MISSING_POS_IND_FILE=""
CONSTRAINT_FILE=""
while getopts g:t:s:l:o:n:m:c:v OPTION
do
  case $OPTION in
    g)  if [ ! -z "${OPTARG}" ];then PHASEFILE=${OPTARG} ;else usage ;fi
        ;;

    l)  if [ ! -z "${OPTARG}" ];then 
          HAP_LIST=${OPTARG}
          HAP_LIST_OUTDISP=${OPTARG}
        else 
          usage
        fi
        ;;

    o)  if [ ! -z "${OPTARG}" ];then HAP_LIST_OUTDISP=${OPTARG} ;else usage ;fi
        ;;
    t)  if [ ! -z "${OPTARG}" ];then TYPE_NUM_ORDERING=${OPTARG} ;else usage ;fi
        ;;
    s)  if [ ! -z "${OPTARG}" ];then SEED=${OPTARG} ;else usage ;fi
        ;;

    n)  if [ ! -z "${OPTARG}" ];then MAX_PARALLEL_DECOMPRESS=${OPTARG} ;else usage ;fi
        ;;
    m)  if [ ! -z "${OPTARG}" ];then MISSING_POS_IND_FILE=${OPTARG} ;else usage ;fi
        ;;
    c)  if [ ! -z "${OPTARG}" ];then CONSTRAINT_FILE=${OPTARG} ;else usage ;fi
        ;;

    v)  VERBOSE=TRUE
        ;;
    \?) usage ;;
  esac
done

#
# output files (for error check)
#

# EACH_DIR
OUTF_SITE_DISTSCORE=site_distScore.txt

# COMBINED_RES_DIR
OUTF_SITE_STATS=results_siteStats.txt.gz
OUTF_SUMMARY_POS=results_siteStats_summary.pos.txt
OUTF_SUMMARY_TXT=sum_site_minus_average.summary.txt.gz
OUTF_SUMMARY_RANGE=sum_site_minus_average.summary.range.txt

PNG_HIST=results_siteStats_hist.png
PNG_ALONG_SEQ=results_siteStats_along_seq.png

#
# check env
#
if [ "${QUEUE_TYPE}" == "SGE" -o "${QUEUE_TYPE}" == "UGE" ]; then
  CHECK=`which qsub`
  if [ "${CHECK}" == "" ]; then
    echo_fail "Error: qsub is not available"
  fi
elif [ "${QUEUE_TYPE}" == "LSF" ]; then
  CHECK=`which bsub`
  if [ "${CHECK}" == "" ]; then
    echo_fail "Error: bsub is not available"
  fi
fi

if [ ! -d "${LOG_DIR}" ]; then
  echo_fail "Error: ${LOG_DIR} doesn't exist.  Please execute setup.sh first."
fi

if [ ! -d "${LIB_DIR}" ]; then
  echo_fail "Error: ${LIB_DIR} doesn't exist.  Please execute setup.sh first."
fi

CHECK=`which R`
if [ "${CHECK}" == "" ]; then
  echo_fail "Error: R is not installed in PATH"
fi

for aa in ${arr_executable_files[@]}
do
  if [ ! -x "${aa}" ]; then
    echo_fail "Environment error: ${aa} doesn't exist or is not executable.  Please execute setup.sh first."
  fi
done

for aa in ${arr_R_files[@]}
do
  if [ ! -f "${aa}" ]; then
    echo_fail "Environment error: ${aa} doesn't exist"
  fi
done

#
# check args
#
if [ $# -lt 2 ] ; then
  usage
fi

if [ ! -f "${PHASEFILE}" ]; then
  echo_fail "Error: ${PHASEFILE} doesn't exist"
fi

WC_UQ_HAP_LEN=`${PL_CHECK_LINE_LEN} -f ${PHASEFILE} | tail -n +5 | uniq | wc -l`
if [ "${WC_UQ_HAP_LEN}" -gt 1 ]; then
  echo_fail "Error: haplotype sequences with different length are found in ${PHASEFILE}"
fi

NUM_IND=`head -2 ${PHASEFILE} | tail -1`

if [ ! -f "${HAP_LIST}" ]; then
  echo_fail "Error: ${HAP_LIST} doesn't exist"
fi

if [ "${MISSING_POS_IND_FILE}" != "" ]; then
  if [ ! -f "${MISSING_POS_IND_FILE}" ]; then
    echo_fail "Error: ${MISSING_POS_IND_FILE} doesn't exist"
  fi

  MAX_MISSING=`awk '{print $1}' ${MISSING_POS_IND_FILE} | uniq -c |sort -n -r |head -1 | awk '{print $1}'`
  if [ "${MAX_MISSING}" -ge "${NUM_IND}" ]; then
    echo_fail "Error: too many missng data, with maximum=${MAX_MISSING} "
  fi
fi

#
# check ${HAP_LIST}
#
CHECK_CR=`head -1 "${HAP_LIST}" | od -c | grep "\r"`
if [ "${CHECK_CR}" != "" ]; then
  perl -i -pe 's/\r//g' ${HAP_LIST}
fi

CHECK_PIPE=`cat "${HAP_LIST}" | grep '|'`
if [ "${CHECK_PIPE}" != "" ]; then
  perl -i -pe 's/\|/_/g' ${HAP_LIST}
fi

WC_HAP_LIST=`wc -l ${HAP_LIST} | awk '{print $1}'`
if [ "${WC_HAP_LIST}" -ne "${NUM_IND}" ]; then
  echo_fail "Error: The number of rows of ${HAP_LIST} must be ${NUM_IND}, but ${WC_HAP_LIST}"
fi

# 
# check ${HAP_LIST_OUTDISP}
#
if [ "${HAP_LIST}" != "${HAP_LIST_OUTDISP}" ]; then
  CHECK_CR=`head -1 "${HAP_LIST_OUTDISP}" | od -c | grep "\r"`
  if [ "${CHECK_CR}" != "" ]; then
    perl -i -pe 's/\r//g' ${HAP_LIST_OUTDISP}
  fi

  CHECK_PIPE=`cat "${HAP_LIST_OUTDISP}" | grep '|'`
  if [ "${CHECK_PIPE}" != "" ]; then
    perl -i -pe 's/\|/_/g' ${HAP_LIST_OUTDISP}
  fi

  WC_HAP_LISTDISP=`wc -l ${HAP_LIST_OUTDISP} | awk '{print $1}'`
  if [ "${WC_HAP_LISTDISP}" -ne "${NUM_IND}" ]; then
    echo_fail "Error: The number of rows of ${HAP_LIST_OUTDISP} must be ${NUM_IND}, but ${WC_HAP_LISTDISP}"
  fi

  if [ ! -s "${HAP_LIST_OUTDISP}.names" ]; then
    awk '{print $1}' ${HAP_LIST}         | sort > ${HAP_LIST}.names
    awk '{print $1}' ${HAP_LIST_OUTDISP} | sort > ${HAP_LIST_OUTDISP}.names
    DIFF=`diff ${HAP_LIST}.names ${HAP_LIST_OUTDISP}.names`
    if [ "${DIFF}" != "" ]; then
      echo_fail "Error: difference of names between ${HAP_LIST} and ${HAP_LIST_OUTDISP}: ${DIFF}"
    fi
  fi
fi

#
# check ${MISSING_POS_IND_FILE}
#
if [ "${MISSING_POS_IND_FILE}" != "" ]; then
  CHECK_CR=`head -1 "${MISSING_POS_IND_FILE}" | od -c | grep "\r"`
  if [ "${CHECK_CR}" != "" ]; then
    perl -i -pe 's/\r//g' ${MISSING_POS_IND_FILE}
  fi

  if [ ! -s "${MISSING_POS_IND_FILE}.names" ]; then
    awk '{print $2}' ${MISSING_POS_IND_FILE} | sort -u > ${MISSING_POS_IND_FILE}.names

    if [ ! -s "${HAP_LIST}.names" ]; then
      awk '{print $1}' ${HAP_LIST}  | sort > ${HAP_LIST}.names
    fi

    DIFF=`comm -23 ${MISSING_POS_IND_FILE}.names ${HAP_LIST}.names`
    if [ "${DIFF}" != "" ]; then
      echo_fail "Error: some strain names in ${MISSING_POS_IND_FILE} are not found in ${HAP_LIST}"
    fi
  fi
fi

#
# prepare
#
OUT_PREFIX_BASE=`echo ${PHASEFILE} | perl -pe 's/^.*\///g' | perl -pe "s/\.hap//g"`

TGZ_ORDER_PAINTINGS=${OUT_PREFIX_BASE}_orderedS${SEED}_both_paintings.tgz

# output files of STEP2
ORDER_HAP_LIST=${OUT_PREFIX_BASE}_orderedS${SEED}_hap.list
ORDER_DIR_LIST=${OUT_PREFIX_BASE}_orderedS${SEED}_rnd_1_${TYPE_NUM_ORDERING}_dirs.list
ORDER_STRAIN_LIST=${OUT_PREFIX_BASE}_orderedS${SEED}_rnd_1_${TYPE_NUM_ORDERING}_dirs_strainOrder.list

#declare -a arr_STAMP

cwd=`dirname $0` 
cd $cwd

################################################################################################################

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# create (uniform) recombination map
# estimation of Ne
#   according to http://paintmychromosomes.com/
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
STEP=1

get_stamp ${STEP}
#arr_STAMP=("${arr_STAMP[@]}" "${STAMP}")
disp_punctuate ${STEP} ${STAMP}

NUM_IND=`head -2 ${PHASEFILE} | tail -1 | perl -pe 's/\n//g'`

OUT_DIR_linked1_est=${OUT_PREFIX_BASE}_linked_out1_respartsEM
if [ ! -d "${OUT_DIR_linked1_est}" ]; then
  mkdir ${OUT_DIR_linked1_est}
fi

#
# (uniform) recombination map
#
RECOMB_FNANE=${OUT_PREFIX_BASE}_linked.rec
if [ -s "${RECOMB_FNANE}" ]; then
  echo "${RECOMB_FNANE} already exists.  Skipped."
else
  perl ${PL_MAKE_RECMAP} ${PHASEFILE} ${RECOMB_FNANE}
fi

#
# EM inference of Ne (recombination rate) on a subset of data 
#
globalparams=""
N_e=""
N_e_FNAME=${OUT_PREFIX_BASE}_linked.neest

if [ -s "${N_e_FNAME}" ]; then
  echo "${N_e_FNAME} already exists.  Skipped."
else

  #
  # check format of hap file 
  #
  CMD="${EXE_PAINT} -a 1 1 -j -g $PHASEFILE -r $RECOMB_FNANE -o ${OUT_DIR_linked1_est}/formatcheck > /dev/null 2>&1 &"
  echo ${CMD}
  eval ${CMD}
  PID=$!
  if [ $? -ne 0 ]; then 
    echo_fail "Error: format of $PHASEFILE is wrong. "
  else
    kill ${PID}
    wait ${PID} > /dev/null 2>&1
    if ls ${OUT_DIR_linked1_est}/formatcheck* &> /dev/null; then
      /bin/rm -f ${OUT_DIR_linked1_est}/formatcheck*
    fi
  fi

  #
  # start
  #
  INCREMENT=1
  if [ "${NUM_IND}" -gt 30 ]; then
    INCREMENT=`expr $NUM_IND / 30`
    let INCREMENT=${INCREMENT}+1
  fi

  for ind in `seq 1 ${INCREMENT} ${NUM_IND}`
  do
    out_prefix_ind=${OUT_PREFIX_BASE}_${ind}

    CMD=`returnQSUB_CMD ${STAMP}`
    CMD=${CMD}" <<< '"
    CMD=${CMD}"${EXE_PAINT} -i ${NUM_EM} -in -n 1 -a $ind $ind -j -g $PHASEFILE -r $RECOMB_FNANE -o ${OUT_DIR_linked1_est}/${out_prefix_ind}"
    CMD=${CMD}"'"
    echo ${CMD}
    eval ${CMD}
    if [ $? -ne 0 ]; then 
      echo_fail "Execution error: ${CMD} (step${STEP}_1) "
    fi
  done

  OUT_DIR=${OUT_DIR_linked1_est}
  wait_until_finish "${STAMP}"

  # remove unnecessary files
  ls ${OUT_DIR_linked1_est}/${OUT_PREFIX_BASE}*.out | grep -v EMprobs | xargs rm
  ls ${OUT_DIR_linked1_est}/${OUT_PREFIX_BASE}*.copyprobsperlocus.out.gz | xargs rm

  # summarize *.EMprobs.out files and estimate Ne
  let WC_CONVERGED=${NUM_EM}+2
  for aa in `wc -l ${OUT_DIR_linked1_est}/${OUT_PREFIX_BASE}*.EMprobs.out | grep -v " ${WC_CONVERGED} " | grep -v total | awk '{print $2}'`
  do
    if [ -f ${aa} ]; then
      /bin/mv ${aa} ${aa}.notconverged
    fi
  done
  CMD="${PL_ESTIMATE_Ne} -o ${N_e_FNAME} ${OUT_DIR_linked1_est}/${OUT_PREFIX_BASE}*.EMprobs.out "
  echo "executing ${PL_ESTIMATE_Ne} -o ${N_e_FNAME} ${OUT_DIR_linked1_est}/${OUT_PREFIX_BASE}*.EMprobs.out ... "
  #echo ${CMD}
  eval ${CMD}
  if [ $? -ne 0 ]; then 
    echo_fail "Execution error: ${CMD} (step${STEP}_2) "
  fi
fi


globalparams=`cat $N_e_FNAME` # $N_e_FNAME contains the commands to tell ChromoPainter about both Ne and global mutation rate, e.g. "-n 10000 -M 0.01".
N_e=`cat $N_e_FNAME | awk '{print $2}' `

if [ "${globalparams}" == "" ]; then 
  echo_fail "Error (step${STEP}): global params infered by EM are empty"
fi

move_log_files "${STAMP}"

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# prepare ordered *.hap files 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
STEP=2

get_stamp ${STEP}
#arr_STAMP=("${arr_STAMP[@]}" "${STAMP}")
disp_punctuate ${STEP} ${STAMP}

arr_s2_target_ordering=()

#
# check unfinished directories
#
i_ordering=1
while [ "${i_ordering}" -le "${TYPE_NUM_ORDERING}"  ]
do
  echo "checking ordering (forward and reverse) ${i_ordering} ..."
  EACH_DIR_PREFIX=$(printf %s_orderedS%s_rnd%02d ${OUT_PREFIX_BASE} ${SEED} ${i_ordering})

  FINISHED_FLAG=FALSE

  #
  # whether preparation of hap files (STEP2) finished or not
  #
  if ls ${EACH_DIR_PREFIX}_forward/*.hap &> /dev/null; then
    CHECK=`ls ${EACH_DIR_PREFIX}_???????/*.hap | wc -l`
    let CHECK=${CHECK}+2
    
    let CORRECT=${NUM_IND}+${NUM_IND}
    if [ "${CHECK}" == "${CORRECT}" ]; then
      FINISHED_FLAG=TRUE
    fi
  #
  # whether painting (STEP3) finished or not
  #
  elif ls ${EACH_DIR_PREFIX}_forward/*.copyprobsperlocus.out.gz &> /dev/null; then
    CHECK=`ls ${EACH_DIR_PREFIX}_???????/*.copyprobsperlocus.out.gz | wc -l`
    let CHECK=${CHECK}+2

    let CORRECT=${NUM_IND}+${NUM_IND}
    if [ "${CHECK}" == "${CORRECT}" ]; then
      FINISHED_FLAG=TRUE
    fi
  #
  # whether ${GZ_CAT_COPYPROB_EACH_DIR}.?? (STEP4) exist or not
  #
  elif   ls ${EACH_DIR_PREFIX}_forward/${GZ_CAT_COPYPROB_EACH_DIR}.?? &> /dev/null; then
    CHECK=`ls ${EACH_DIR_PREFIX}_???????/${GZ_CAT_COPYPROB_EACH_DIR}.?? | wc -l`
    
    let CORRECT=${NUM_SPLITN}+${NUM_SPLITN}
    if [ "${CHECK}" == "${CORRECT}" ]; then
      FINISHED_FLAG=TRUE
    fi
  #
  # incomplete decompressing (STEP4) 
  #
  elif ls ${EACH_DIR_PREFIX}_forward/*.copyprobsperlocus.out_?? &> /dev/null; then
    CHECK=`ls ${EACH_DIR_PREFIX}_???????/*.copyprobsperlocus.out_?? | wc -l`
    
    let CORRECT=${NUM_IND}-1
    CORRECT=`expr ${CORRECT} \* ${NUM_SPLITN}`
    let CORRECT=${CORRECT}+${CORRECT}
    if [ "${CHECK}" == "${CORRECT}" ]; then
      FINISHED_FLAG=TRUE
    fi
    
    for aa in `ls ${EACH_DIR_PREFIX}_???????/*.copyprobsperlocus.out_??`
    do
      if [ ! -s "${aa}" ]; then
        FINISHED_FLAG=FALSE
        break
      fi
    done
  fi


  #
  # store unfinished orderings
  #
  if [ "${FINISHED_FLAG}" == "FALSE" ]; then
    arr_s2_target_ordering+=(${i_ordering})
  fi
  
  let i_ordering=${i_ordering}+1
done


#
# execute
#
for i_s2_target_ordering in ${arr_s2_target_ordering[@]}
do
  #
  # prepare CMD
  #
  CMD=""
  
  ARRAY_S=1
  ARRAY_E=`wc -l ${HAP_LIST} | awk '{print $1}'`

#    # arrayjob for i_recipient (can be slower if the num. of available cores is limited)
##    CMD=`returnQSUB_CMD ${STAMP} ${ARRAY_S} ${ARRAY_E}`
#    CMD=${CMD}" ${SH_RANDOMIZE}"
#    CMD=${CMD}" -h ${PHASEFILE}"
#    CMD=${CMD}" -p ${OUT_PREFIX_BASE}"
#    CMD=${CMD}" -l ${HAP_LIST}"
#    CMD=${CMD}" -o ${HAP_LIST_OUTDISP}"
#    CMD=${CMD}" -t ${i_s2_target_ordering}"
#    CMD=${CMD}" -s ${SEED}" 

  # qsub for each_ordering
  QSUB_FILE=${STAMP}_${OUT_PREFIX_BASE}_orderedS${SEED}_${i_s2_target_ordering}.sh
/bin/cat  > ${QSUB_FILE} << EOF
#! /bin/bash
#$ -S /bin/bash
EOF
  for i_recipient in `seq ${ARRAY_S} ${ARRAY_E}`
  do
/bin/cat >> ${QSUB_FILE} << EOF
${SH_RANDOMIZE} \
  -h ${PHASEFILE} \
  -p ${OUT_PREFIX_BASE} \
  -l ${HAP_LIST} \
  -o ${HAP_LIST_OUTDISP} \
  -t ${i_s2_target_ordering} \
  -i ${i_recipient} \
  -s ${SEED} 
EOF
  done

  chmod 755 ${QSUB_FILE}
  CMD=`returnQSUB_CMD ${STAMP} `
  CMD=${CMD}" ./${QSUB_FILE}" # " <<< /bin/bash " doesn't work in LSF

  #
  # submit
  #
  echo ${CMD}
  QSUB_MSG=`${CMD}`
  if [ $? -ne 0 ]; then 
    echo_fail "Execution error: ${CMD} (step${STEP}) "
  fi

done

wait_until_finish "${STAMP}"

if [ ls ${STAMP}_${OUT_PREFIX_BASE}_orderedS${SEED}_*.sh &> /dev/null; then
  CMD="/bin/rm -f ${STAMP}_${OUT_PREFIX_BASE}_orderedS${SEED}_*.sh"
  eval ${CMD}
  if [ $? -ne 0 ]; then 
    echo_fail "Error: ${CMD}  "
  fi
fi


#
# create ${ORDER_HAP_LIST} as a preparation for the next step (painting as arrayjobs)
#
echo "listing up all .hap files to be processed into ${ORDER_HAP_LIST} ..."
/bin/cat /dev/null > ${ORDER_HAP_LIST}

i_ordering=1
while [ "${i_ordering}" -le "${TYPE_NUM_ORDERING}"  ]
do
  EACH_DIR_PREFIX=$(printf %s_orderedS%s_rnd%02d ${OUT_PREFIX_BASE} ${SEED} ${i_ordering})

  if ls ${EACH_DIR_PREFIX}_*/*.hap &> /dev/null; then
    CMD="ls ${EACH_DIR_PREFIX}_*/*.hap >> ${ORDER_HAP_LIST}"
    #echo ${CMD}
    eval ${CMD}
    if [ $? -ne 0 ]; then 
      echo_fail "Error: ${CMD}  "
    fi
  fi

  let i_ordering=${i_ordering}+1
done

wc -l ${ORDER_HAP_LIST}

#
# two list files
#
if [ ! -s "${ORDER_DIR_LIST}" ]; then
  CMD="find ./ -maxdepth 1 -type d -name ${OUT_PREFIX_BASE}_orderedS${SEED}_rnd\* | grep -v results | perl -pe 's/^\.\///g' | sort > ${ORDER_DIR_LIST}"
  #echo ${CMD}
  eval ${CMD}
  if [ $? -ne 0 ]; then 
    echo_fail "Error: ${ORDER_DIR_LIST} "
  fi
  echo "${ORDER_DIR_LIST} was created"
fi

if [ ! -s "${ORDER_STRAIN_LIST}" ]; then
  CMD="/bin/cat ${ORDER_DIR_LIST} | perl -pe 's/\n/.strainOrder\n/g' > ${ORDER_STRAIN_LIST}"
  #echo ${CMD}
  eval ${CMD}
  if [ $? -ne 0 ]; then 
    echo_fail "Error: ${ORDER_STRAIN_LIST} "
  fi
  echo "${ORDER_STRAIN_LIST} was created"
fi


#
# cleaning tmp files
#
if [ -f "${HAP_LIST}.names" ]; then
  /bin/rm -f ${HAP_LIST}.names
fi

if [ -f "${HAP_LIST_OUTDISP}.names" ]; then
  /bin/rm -f ${HAP_LIST_OUTDISP}.names
fi


move_log_files "${STAMP}"


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# execute chromopainter for each ordered recipient haplotype
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
STEP=3

get_stamp ${STEP}
#arr_STAMP=("${arr_STAMP[@]}" "${STAMP}")
disp_punctuate ${STEP} ${STAMP}

TARGET_HAP_FNANE="target_hap.list"

#NUM_HAP=`head -2 ${PHASEFILE} | tail -1`
#NUM_SITE=`head -3 ${PHASEFILE} | tail -1`
#let NUM_ROW_COPYPROB=${NUM_SITE}+2

while read EACH_DIR
do
  #
  # prepare $TARGET_HAP_LIST, $NUM_TARGET_HAP
  #
  NUM_TARGET_HAP=0
  TARGET_HAP_LIST="${EACH_DIR}/${TARGET_HAP_FNANE}"
  /bin/cat /dev/null > ${TARGET_HAP_LIST}
  echo "preparing ${TARGET_HAP_LIST} ... "

  CHECK_GZ_CAT_COPYPROB_EACH_DIR=0
  if [ -f "${EACH_DIR}/${GZ_CAT_COPYPROB_EACH_DIR}" ]; then
    CHECK_GZ_CAT_COPYPROB_EACH_DIR=`gzip -dc ${EACH_DIR}/${GZ_CAT_COPYPROB_EACH_DIR} | head | wc -l`
  fi
  # skip this ordered directory if there is ${GZ_CAT_COPYPROB_EACH_DIR} created by the next step
  if [ ${CHECK_GZ_CAT_COPYPROB_EACH_DIR} -gt 0 ]; then
    echo "  painting in ${EACH_DIR} is skipped because there is already ${GZ_CAT_COPYPROB_EACH_DIR}"
  else
    for EACH_HAP in `grep ${EACH_DIR} ${ORDER_HAP_LIST}`
    do
      EACH_COPYPROB_GZ=`echo ${EACH_HAP} | perl -pe 's/\.hap$/.copyprobsperlocus.out.gz/g'`

      target_flag=0
      if [ ! -s "${EACH_COPYPROB_GZ}" ]; then
        target_flag=1
      else
        CHECK_HEAD=`gzip -dc ${EACH_COPYPROB_GZ} | head | wc -l`
        if [ $? -ne 0 -o "${CHECK_HEAD}" -eq 0 ]; then
          target_flag=1
          
          /bin/rm -f ${EACH_COPYPROB_GZ}
          echo "incomplete ${EACH_COPYPROB_GZ} was removed"
        fi
      fi
      # otherwise, execute painting of unfinished hap files in this ordered directory 
      if [ "${target_flag}" -eq 1 ]; then
        echo ${EACH_HAP} >> ${TARGET_HAP_LIST}
        let NUM_TARGET_HAP=${NUM_TARGET_HAP}+1
      else
        echo "${EACH_COPYPROB_GZ} already exists and is not empty.  Skipped."
      fi
    done
  fi

  #
  # execute paining of the target hap files
  #
  if [ "${NUM_TARGET_HAP}" -eq 0 ]; then
    echo "${EACH_DIR} was skipped because there is no hap file to be painted."
  else
    ARRAY_S=1
    ARRAY_E=${NUM_TARGET_HAP}
    
    CMD=`returnQSUB_CMD ${STAMP} ${ARRAY_S} ${ARRAY_E}`
    #CMD=${CMD}" -t 1:${NUM_TARGET_HAP} "
    CMD=${CMD}" ${SH_PAINT_QSUB}"
    CMD=${CMD}"  -r ${RECOMB_FNANE}"
    CMD=${CMD}"  -n ${N_e_FNAME}"
    CMD=${CMD}"  -l ${TARGET_HAP_LIST}"

    echo ${CMD}
    QSUB_MSG=`${CMD}`
    if [ $? -ne 0 ]; then 
      echo_fail "Execution error: ${CMD} (step${STEP}) "
    fi
    #QSUB_ID=`echo ${QSUB_MSG} | perl -pe 's/ \(.*$//g' | perl -pe 's/^.* //g' | perl -pe 's/\..*$//g'`
  fi

done < ${ORDER_DIR_LIST}

wait_until_finish "${STAMP}"


echo "checking whether there is any incomplete output file (.copyprobsperlocus.out.gz) of the step${STEP} ..."

i_failed=0
while read EACH_DIR
do
  echo "${EACH_DIR}"
  i_HAP=0
  TARGET_HAP_LIST="${EACH_DIR}/${TARGET_HAP_FNANE}"
  while read EACH_HAP
  do
    let i_HAP=${i_HAP}+1
    EACH_COPYPROB_GZ=`echo ${EACH_HAP} | perl -pe 's/\.hap/.copyprobsperlocus.out.gz/g'`
    
    if [ -f "${EACH_COPYPROB_GZ}" ]; then
      CHECK_HEAD=`gzip -dc ${EACH_COPYPROB_GZ} | head | wc -l`
      if [ $? -ne 0 -o "${CHECK_HEAD}" -eq 0 ]; then
        echo "painting of ${EACH_HAP} failed, because ${EACH_COPYPROB_GZ} is an incomplete file"
        let i_failed=${i_failed}+1
        
        /bin/rm -f ${EACH_COPYPROB_GZ}
        echo "incomplete ${EACH_COPYPROB_GZ} was removed"
      fi
    fi

  done < "${TARGET_HAP_LIST}"
done < ${ORDER_DIR_LIST}

if [ "${i_failed}" -gt 0 ]; then
  echo_fail "There are ${i_failed} failed jobs.  Please execute this program again (already finished jobs will be skipped)."
fi

move_log_files "${STAMP}"



################################################################################################################

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Formatting output files of ChromoPainter
#
#   for each ordering,
#     split
#     cat to ${GZ_CAT_COPYPROB_EACH_DIR}
#
#     it can require a large temporary disk in each ordering
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
STEP=4

get_stamp ${STEP}
#arr_STAMP=("${arr_STAMP[@]}" "${STAMP}")
disp_punctuate ${STEP} ${STAMP}

TARGET_GZ_FNANE="target_gz.list"

declare -a arr_dirs_for_msort

while read EACH_DIR
do
  #
  # prepare $TARGET_HAP_LIST, $NUM_TARGET_HAP
  #
  TARGET_GZ_LIST="${EACH_DIR}/${TARGET_GZ_FNANE}"
  /bin/cat /dev/null > "${TARGET_GZ_LIST}"
  echo "preparing ${TARGET_GZ_LIST} ... "

  if ls ${EACH_DIR}/*copyprobsperlocus.out.gz &> /dev/null; then
    CMD="ls ${EACH_DIR}/*copyprobsperlocus.out.gz > ${TARGET_GZ_LIST}"
    #echo ${CMD}
    eval ${CMD}
    if [ $? -ne 0 ]; then 
      echo_fail "Execution error: ${CMD} (step${STEP}) "
    fi
  fi

  NUM_TARGET_GZ=`wc -l ${TARGET_GZ_LIST} | awk '{print $1}'`

  CHECK_GZ_CAT_COPYPROB_EACH_DIR=0
  if [ -f "${EACH_DIR}/${GZ_CAT_COPYPROB_EACH_DIR}" ]; then
    CHECK_GZ_CAT_COPYPROB_EACH_DIR=`gzip -dc ${EACH_DIR}/${GZ_CAT_COPYPROB_EACH_DIR} | head | wc -l`
  fi
  
  #
  # skip this ordered directory if there is already ${GZ_CAT_COPYPROB_EACH_DIR}.??
  if ls ${EACH_DIR}/${GZ_CAT_COPYPROB_EACH_DIR}.?? &> /dev/null; then
    CHECK=`ls ${EACH_DIR}/${GZ_CAT_COPYPROB_EACH_DIR}.?? | wc -l`

    if [ "${CHECK}" == "${NUM_SPLITN}" ]; then
      echo "msort in ${EACH_DIR} is skipped because there are already ${NUM_SPLITN} sorted files"
    fi
  else
    arr_dirs_for_msort=("${arr_dirs_for_msort[@]}" "${EACH_DIR}")
    
    #
    # if there are .copyprobsperlocus.out.gz files to be decompressed,
    #
    # decompress
    # split
    #
    ARRAY_S=1
    ARRAY_E=${NUM_TARGET_GZ}

    if [ "${NUM_TARGET_GZ}" -gt 0 ]; then
      CMD=`returnQSUB_CMD ${STAMP} ${ARRAY_S} ${ARRAY_E}`
      #CMD=${CMD}" -t 1:${NUM_TARGET_GZ} "
      CMD=${CMD}" ${SH_DECOMPRESS_SORT_SPLIT_EACH_ORDERING}"
      CMD=${CMD}"  -l ${TARGET_GZ_LIST}"

      echo ${CMD}
      QSUB_MSG=`${CMD}`
      if [ $? -ne 0 ]; then 
        echo_fail "Execution error: ${CMD} (step${STEP}) "
      fi
    else
      echo "${SH_DECOMPRESS_SORT_SPLIT_EACH_ORDERING} was not executed for ${EACH_DIR} (no *.gz file) "
    fi

    while :
    do
      #
      # when the number of dirs to be processed > ${MAX_PARALLEL_DECOMPRESS}
      #
      if [ "${#arr_dirs_for_msort[@]}" -ge "${MAX_PARALLEL_DECOMPRESS}" ]; then
        #
        # wait decompression (for dirs submitted above)
        #
        wait_until_finish "${STAMP}"

        #
        # then submit msort for each decompressed dir
        #
        submit_msort_each_ordering "${STAMP}"
        #
        # wait until the submitted msort jobs are finished
        #
        wait_until_finish "${STAMP}"
      else 
      #
      # otherwise, proceed to qsub of the next ordering
      #
        break
      fi

      sleep 10
    done

  fi
done < ${ORDER_DIR_LIST}

#
# wait decompression
#   (used only if ${MAX_PARALLEL_DECOMPRESS} was not used above)
#
wait_until_finish "${STAMP}"



#
# submit msort for dirs recorded in arr_dirs_for_msort
#    if ${MAX_PARALLEL_DECOMPRESS} was not used above, 
#      all dirs are submitted here
#
#    if ${MAX_PARALLEL_DECOMPRESS} was used above, 
#      no dir is submitted here 
#      because arr_dirs_for_msort is already empty
#
submit_msort_each_ordering "${STAMP}"
#
# wait until the submitted msort jobs are finished
#
wait_until_finish "${STAMP}"



#
# cat ${GZ_CAT_COPYPROB_EACH_DIR}.?? in each dir into ${GZ_CAT_COPYPROB_EACH_DIR} 
#   for calculating the average in the postprocessing below
#
while read EACH_DIR
do
  if [ ! -s "${EACH_DIR}/${GZ_CAT_COPYPROB_EACH_DIR}" ]; then
    if ls ${EACH_DIR}/${GZ_CAT_COPYPROB_EACH_DIR}.?? &> /dev/null; then
      CMD=`returnQSUB_CMD ${STAMP}`
      CMD=${CMD}" <<< '"
      CMD=${CMD}" /bin/cat ${EACH_DIR}/${GZ_CAT_COPYPROB_EACH_DIR}.?? > ${EACH_DIR}/${GZ_CAT_COPYPROB_EACH_DIR}"
      CMD=${CMD}"'"

      echo ${CMD}
      eval ${CMD}
      if [ $? -ne 0 ]; then 
        echo_fail "Execution error: ${CMD} (step${STEP})"
      fi
    fi
  fi
done < ${ORDER_DIR_LIST}

wait_until_finish "${STAMP}"


#
# check ${GZ_CAT_COPYPROB_EACH_DIR} in each dir
#
while read EACH_DIR
do
  if [ ! -s "${EACH_DIR}/${GZ_CAT_COPYPROB_EACH_DIR}" ]; then
    echo_fail "Error: ${EACH_DIR}/${GZ_CAT_COPYPROB_EACH_DIR} doesn't exist or empty"
  else
    CHECK_HEAD=`gzip -dc ${EACH_DIR}/${GZ_CAT_COPYPROB_EACH_DIR} | head | wc -l`
    if [ $? -ne 0 -o "${CHECK_HEAD}" -eq 0 ]; then
      echo_fail "Error: ${EACH_DIR}/${GZ_CAT_COPYPROB_EACH_DIR} is an incomplete file"
    fi
  fi

  #
  # temporary files, if they remain for some reason
  #
  if ls ${EACH_DIR}/sort?????? &> /dev/null; then # old (non-arrayjob) version
    /bin/rm -f ${EACH_DIR}/sort??????
  fi

  # to be removed in chromopainter_linkage_orderings_arrayjob.sh
  if ls ${EACH_DIR}/*.hap &> /dev/null; then
    /bin/rm -f ${EACH_DIR}/*.hap
  fi

  # to be removed in mort_each_ordering.pl
  if ls ${EACH_DIR}/*copyprobsperlocus.out_?? &> /dev/null; then # arrayjob version
    /bin/rm -f ${EACH_DIR}/*copyprobsperlocus.out_??
  fi

done < ${ORDER_DIR_LIST}

move_log_files "${STAMP}"



#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# postprocessing (l1,l2)
#   calculate average, and distance to the average for each ordering
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
STEP=5

get_stamp ${STEP}
#arr_STAMP=("${arr_STAMP[@]}" "${STAMP}")
disp_punctuate ${STEP} ${STAMP}

i_submitted=0

i_ordering=1
while read EACH_DIR
do
  if [ ! -s "${EACH_DIR}/${OUTF_SITE_DISTSCORE}" ]; then # doesn't exist or empty
    submit_calcAveDist_ordering "${i_ordering}"
    let i_submitted=${i_submitted}+1
  else
    echo "step${STEP} of ${EACH_DIR} was skipped, because ${OUTF_SITE_DISTSCORE} already exists there";
  fi
  let i_ordering=${i_ordering}+1
done < ${ORDER_DIR_LIST}


if [ "${i_submitted}" -gt 0 ]; then

  wait_until_finish "${STAMP}"

  #
  # check the output file
  # 
  while read EACH_DIR
  do
    if [ ! -s "${EACH_DIR}/${OUTF_SITE_DISTSCORE}" ]; then
      echo_fail "Error (step${STEP}): ${EACH_DIR}/${OUTF_SITE_DISTSCORE} doesn't exist or empty"
    fi
  done < ${ORDER_DIR_LIST}
  
fi

move_log_files "${STAMP}"


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# postprocessing 3 (l3)
#   combine results of all orderings and produce final results
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
STEP=6

get_stamp ${STEP}
#arr_STAMP=("${arr_STAMP[@]}" "${STAMP}")
disp_punctuate ${STEP} ${STAMP}

COMBINED_RES_DIR=`echo ${ORDER_DIR_LIST} | perl -pe 's/\.list/_results/g'`

SKIP_FLAG=0
if [ -s "${COMBINED_RES_DIR}/${OUTF_SITE_STATS}" ]; then
  if [ -s "${COMBINED_RES_DIR}/${OUTF_SUMMARY_POS}" ]; then
    if [ -s "${COMBINED_RES_DIR}/${OUTF_SUMMARY_TXT}" ]; then
      if [ -s "${COMBINED_RES_DIR}/${OUTF_SUMMARY_RANGE}" ]; then
        DIRS_OK_FLAG=1
        for VISUALIZE_TYPE_DIR in `find ${COMBINED_RES_DIR} -type d -name visualize\*`
        do
          if [ ! -d "${VISUALIZE_TYPE_DIR}" ]; then
            DIRS_OK_FLAG=0
          fi
        done

        if [ "${DIRS_OK_FLAG}" -eq 1 ]; then
          SKIP_FLAG=1
        fi
      fi
    fi
  fi
fi

if [ "${SKIP_FLAG}" -eq 0 ]; then

  CMD=`returnQSUB_CMD ${STAMP}` 
  CMD=${CMD}" ${MAX_MEMORY} " # ${MAX_MEMORY} is used only here
  CMD=${CMD}" <<< '"
  CMD=${CMD}"perl ${PL_SITE_BY_SITE}"
  CMD=${CMD}" -g ${PHASEFILE} "
  CMD=${CMD}" -d ${ORDER_DIR_LIST} "
  CMD=${CMD}" -l ${ORDER_STRAIN_LIST} "
  #if [ "${CONTRAST_MAX}" -gt 0 ]; then
  #  CMD=${CMD}" -c ${CONTRAST_MAX} " 
  #fi
  CMD=${CMD}" -s ${HAP_LIST_OUTDISP} "
  CMD=${CMD}" -r "   # only one difference from the previous step
  CMD=${CMD}"'"
  
  echo ${CMD}
  eval ${CMD}
  if [ $? -ne 0 ]; then 
    echo_fail "Execution error: ${CMD} (step${STEP}) "
  fi
  
  wait_until_finish "${STAMP}"


  if [ ! -s "${COMBINED_RES_DIR}/${OUTF_SITE_STATS}" ]; then
    echo_fail "Error (step${STEP}): ${COMBINED_RES_DIR}/${OUTF_SITE_STATS} doesn't exist or empty "
  fi

  if [ ! -s "${COMBINED_RES_DIR}/${OUTF_SUMMARY_POS}" ]; then
    echo_fail "Error (step${STEP}): ${COMBINED_RES_DIR}/${OUTF_SUMMARY_POS} doesn't exist or empty "
  fi

  if [ `gzip -dc "${COMBINED_RES_DIR}/${OUTF_SUMMARY_TXT}" | wc -l` -le 1 ]; then
    echo_fail "Error (step${STEP}): ${COMBINED_RES_DIR}/${OUTF_SUMMARY_TXT} is empty "
  fi

  if [ ! -s "${COMBINED_RES_DIR}/${OUTF_SUMMARY_RANGE}" ]; then
    echo_fail "Error (step${STEP}): ${COMBINED_RES_DIR}/${OUTF_SUMMARY_RANGE} doesn't exist or empty "
  fi


  echo "The step${STEP} normally finished."

else
  echo "step${STEP} was skipped because all output files are in ${COMBINED_RES_DIR}"
fi

move_log_files "${STAMP}"


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# visualization 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
STEP=7

get_stamp ${STEP}
#arr_STAMP=("${arr_STAMP[@]}" "${STAMP}")
disp_punctuate ${STEP} ${STAMP}

#
# basic plot (histgram & along the sequemce) of all sites
#
if [ ! -s "${COMBINED_RES_DIR}/${PNG_HIST}" -o ! -s "${COMBINED_RES_DIR}/${PNG_ALONG_SEQ}" ]; then
  CMD="R --vanilla --quiet < ${R_MAIN1} --args ${R_LIB_HEATMAP} ${COMBINED_RES_DIR}/${OUTF_SITE_STATS} > /dev/null 2>&1"
  echo ${CMD}
  eval ${CMD}
  if [ $? -ne 0 ]; then 
    echo_fail "Execution error: ${CMD} (step${STEP}, ${R_MAIN1}) "
  fi
else 
  echo "${R_MAIN1} was skipped because there are already ${COMBINED_RES_DIR}/${PNG_HIST} and ${COMBINED_RES_DIR}/${PNG_HIST}"
fi 

#
# if imputed data are specified,
# plot relation between missing count per site and the distance statistic
#
if [ "${MISSING_POS_IND_FILE}" != "" ]; then
  CMD="R --vanilla --quiet < ${R_CHECK_MISSING_STAT} --args ${MISSING_POS_IND_FILE}  ${COMBINED_RES_DIR}/${OUTF_SITE_STATS} > /dev/null 2>&1"
  echo ${CMD}
  eval ${CMD}
  if [ $? -ne 0 ]; then 
    echo_fail "Execution error: ${CMD} (step${STEP}, ${R_CHECK_MISSING_STAT}) "
  fi
fi

#
# heatmaps of summary sites (not executed by default)
#
MIN=`awk '{print $1}' "${COMBINED_RES_DIR}/${OUTF_SUMMARY_RANGE}"`
MAX=`awk '{print $2}' "${COMBINED_RES_DIR}/${OUTF_SUMMARY_RANGE}"`

echo "If you would like to visualize representative sites, please execute the following commands"

for VISUALIZE_TYPE_DIR in `find ${COMBINED_RES_DIR} -type d -name visualize\*`
do
  if [ -d "${VISUALIZE_TYPE_DIR}" ]; then
    ls ${VISUALIZE_TYPE_DIR}/* > ${VISUALIZE_TYPE_DIR}.list

    NUM_TARGET_POS=`wc -l ${VISUALIZE_TYPE_DIR}.list | awk '{print $1}'`

    ARRAY_S=1
    ARRAY_E=${NUM_TARGET_POS}

    CMD=`returnQSUB_CMD ${STAMP} ${ARRAY_S} ${ARRAY_E}`
    #CMD=${CMD}" -t 1:${NUM_TARGET_POS} "
    CMD=${CMD}" ${SH_R_MAIN2}"
    CMD=${CMD}"  -a ${MIN}"
    CMD=${CMD}"  -b ${MAX}"
    CMD=${CMD}"  -l ${VISUALIZE_TYPE_DIR}.list"
 
    echo ${CMD}
    #QSUB_MSG=`${CMD}`
    if [ $? -ne 0 ]; then 
      echo_fail "Execution error: ${CMD} (step${STEP}, ${SH_R_MAIN2}) "
    fi
    #QSUB_ID=`echo ${QSUB_MSG} | perl -pe 's/ \(.*$//g' | perl -pe 's/^.* //g' | perl -pe 's/\..*$//g'`
  else
    echo_fail "Error: ${VISUALIZE_TYPE_DIR} doesn't exist"
  fi
done

#wait_until_finish "${STAMP}"
#move_log_files "${STAMP}"

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

echo "************************************************************** "

date +%Y%m%d_%T
echo "Done. Please look at output files in ${COMBINED_RES_DIR}"

