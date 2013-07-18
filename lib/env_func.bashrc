#########################################################################
# definition of grid engine
#########################################################################

QUEUE_TYPE="SGE"
QSUB_COMMON="qsub -cwd "
#QSUB_COMMON="qsub -cwd -l s_vmem=8G -l mem_req=8"

#QUEUE_TYPE="LSF"
#QSUB_COMMON="bsub "
#QSUB_COMMON="bsub -M 8000000 "


#########################################################################
# common functions
#########################################################################

echo_fail(){
  echo "`date '+%Y%m%d_%H%S'` $*" >&2
  exit 1
}

